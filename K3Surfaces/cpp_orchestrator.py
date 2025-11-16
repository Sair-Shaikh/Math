import os
import re
import uuid
import subprocess
from pathlib import Path
import pandas as pd
from concurrent.futures import ProcessPoolExecutor, as_completed
import time

BASE_PATH = Path("Dataset")
COEFFS_DIR = Path("CppLib/Coeffs")
# CPP_FILE = Path("CppLib/count_cubic_batched.cpp")
CPP_FILE = Path("CppLib/count_quad_batched.cpp")
CSV_FILE = BASE_PATH / "point_counts_quad_2.csv"
CONST_HEADER = Path("CppLib/constants.h")
DATA_FILES = [ 
                Path("Dataset/CppCoeffs/quad_coeffs_table3.txt"),
            #    Path("Dataset/CppCoeffs/cube_coeffs_table_1.txt"), 
            #    Path("Dataset/CppCoeffs/cube_coeffs_table_2.txt"), 
            #    Path("Dataset/CppCoeffs/cube_coeffs_table_3.txt"), 
            #    Path("Dataset/CppCoeffs/cube_coeffs_table_4.txt"), 
            #    Path("Dataset/CppCoeffs/cube_coeffs_table_5.txt"), 
            ]

# MISSING_KEYS = ["7-960699","7-1171368", "7-78040","5-115860","7-1351482","7-258285","7-119622",
#                 "7-420257","7-798571","7-240306","5-794630","7-420220","7-438238","7-618572",
#                 "7-1349484","7-240287","7-79843","7-1351485","7-1349505","7-258291","7-258289",
#                 "7-101626","7-798575","7-240312","5-614635","7-101630","7-60047","7-1351491",
#                 "7-1351511","7-960709","7-101631","7-61850","5-115870","7-600540","7-1349512",
#                 "7-780540","7-978712","7-1171341","7-60049","7-60051","7-61851","7-79851","7-618541",
#                 "5-614798","7-438230","6-300799","6-299447","7-1490780","7-843888","6-297780",
#                 "6-299655","6-300325","6-301358","6-300680","6-297619","6-305232","6-297987",
#                 "7-1496123","6-302989","6-304923","6-305123","6-300784","6-300327","6-301342", 
#                 "7-1171366", "7-61876","5-434693","7-258283","7-438255","7-600567","7-1349499","7-1351499"  
#                 ]


key_pattern = re.compile(r"^Key:\s*(.+)")
corrections_pattern = re.compile(r"^Corrections:\s*\[(.*?)\]")

BATCH_SIZE = 100

def stream_chunks():
    for file in DATA_FILES:
        key = None
        cpp_lines = []
        in_cpp_block = False

        with open(file, "r") as f:
            for line in f:

                key_match = key_pattern.match(line)
                corr_match = corrections_pattern.match(line)

                if key_match:
                    key = key_match.group(1).strip()
                    in_cpp_block = True
                    continue

                elif corr_match:
                    # End of chunk
                    corrections = [int(x) for x in re.findall(r"-?\d+", corr_match.group(1))]
                    in_cpp_block = False

                    if key and cpp_lines:
                        if key.startswith("5"):
                            yield (key, cpp_lines, corrections)
                        key, cpp_lines = None, []
                        cpp_lines = []
                    continue

                elif in_cpp_block:
                    if "ABC" not in line and line.strip():
                        line = line.strip()
                        if not line.endswith("\\"): 
                            if not line.endswith(";"): # First chunk ends here
                                line += "; \\[SPLIT]"
                            else:
                                line += " \\"
                        cpp_lines.append(line)
                else:
                    raise NotImplementedError

def stream_batched_chunks():
    BATCH_SIZE = 100

    results = {}
    count = 0
    out_lines = []
    out_lines_temp = []
    out_lines.append("// Auto-generated batched header\n")
    out_lines.append("#define BATCHED_ABC \\\n")
    out_lines_temp.append("#define BATCHED_ABC2 \\\n")
    for key, macro, corr in stream_chunks():
        body = "\n\t".join(macro)
        j = count % BATCH_SIZE
        results[j] = (key, corr)

        body_mod = re.sub(r'\bA\s*=', f"As[{j}] =", body)
        body_mod = re.sub(r'\bB\s*=', f"Bs[{j}] =", body_mod)
        body_mod = re.sub(r'\bC\s*=', f"Cs[{j}] =", body_mod)
        # body_mod = re.sub(r'\bD\s*=', f"Ds[{j}] =", body_mod)

        first, second, _ = body_mod.split(sep="[SPLIT]")
        first, second = first.strip(), second.strip()

        out_lines.append(f"        /* Key: {key} */ \\\n        {first}\n")
        out_lines_temp.append(f"        /* Key: {key} */ \\\n        {second}\n")

        count += 1
        if not (count % BATCH_SIZE):
            # Write to file
            file_text = "".join(out_lines) + "\n\n\n" + "".join(out_lines_temp)
            yield (results, file_text)
            
            # Reset
            out_lines = []
            out_lines_temp = []
            out_lines.append("// Auto-generated batched header\n")
            out_lines.append("#define BATCHED_ABC \\\n")
            out_lines_temp.append("#define BATCHED_ABC2 \\\n")
            results = {}
    
    # Final batch, if there is any
    # Write to file
    if len(out_lines) > 2: 
        file_text = "".join(out_lines) + "\n\n\n" + "".join(out_lines_temp)
        yield (results, file_text)

def process_chunk(chunk):
    results, text = chunk

    # Write header with coefficients
    id = uuid.uuid4()
    header_path = COEFFS_DIR / f"coeffs_{id}.h"
    header_path.write_text(text)

    executable = COEFFS_DIR / f"exec_{id}"

    try: 
        compile_cmd = ["g++","-std=c++17",f"-include{header_path}", str(CPP_FILE),"-o", str(executable),"-O2", f"-DBATCH_SIZE={BATCH_SIZE}", "-DEXT_COEFFS=1"]
        subprocess.run(compile_cmd, check=True)

        # Run executable and capture output
        run_result = subprocess.run(
            [str(executable)],
            capture_output=True,
            text=True,
            check=True,
        )
        output_str = run_result.stdout.strip()
        try:
            numbers = [int(n) for n in output_str.split()]
            curr_results = { k : [numbers[i*BATCH_SIZE + j]+corrs[i] for i in range(11)] for j, (k, corrs) in results.items()} 

        except ValueError:
            raise ValueError(f"Result could not be parsed correctly: {output_str[:50]}")

    finally: 
        # Clean up generated files
        pass
        if header_path.exists():
            header_path.unlink()
        if executable.exists():
            executable.unlink()
    return curr_results

def process_chunks_in_batches(chunk_generator, max_workers=None):
    if max_workers is None:
        max_workers = os.cpu_count()

    all_results = {}
    progress = 0
    batch = []

    start_time = time.time()
    with ProcessPoolExecutor(max_workers=max_workers) as executor:
        futures = set()

        breaker = False
        for chunk in chunk_generator:
            # Submit new task
            futures.add(executor.submit(process_chunk, chunk))
            batch.append(chunk)

            # If batch is full, wait for some futures to complete
            if len(futures) >= max_workers:
                done, futures = wait_some(futures)
                for results in done:
                    all_results = {**all_results, **results}

                    progress += 1
                    if progress % 10 == 0:
                        print(f"Progress: {progress*100}")
                        # breaker = True

                    if progress % 100 == 0:
                        df = pd.DataFrame.from_dict(all_results, orient="index")
                        df.to_csv(CSV_FILE,  mode="a", index=True)
                        all_results = {}
                        print(f"Total Time Taken: {time.time()-start_time}")

            if breaker: 
                break
                        
        # Process remaining futures
        while futures:
            done, futures = wait_some(futures)
            for results in done:
                all_results = {**all_results, **results}
                progress += 1
                df = pd.DataFrame.from_dict(all_results, orient="index")
                df.to_csv(CSV_FILE, mode="a", index=True)
                all_results = {}
        if all_results: 
            for results in done:
                all_results = {**all_results, **results}
                progress += 1
                df = pd.DataFrame.from_dict(all_results, orient="index")
                df.to_csv(CSV_FILE, mode="a", index=True)
                all_results = {}


    return all_results

def wait_some(futures):
    """Wait for at least one future to finish, return (done_results, remaining_futures)"""
    done_futures, remaining_futures = set(), set()
    for future in as_completed(futures):
        done_futures.add(future)
        remaining_futures = futures - done_futures
        done_futures
        break  # only wait for one to complete
    results = [f.result() for f in done_futures]
    return results, remaining_futures

if __name__ == "__main__":

    progress = 0
    process_chunks_in_batches(stream_batched_chunks(), max_workers=8)