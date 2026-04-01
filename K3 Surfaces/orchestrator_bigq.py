import os
import re
import uuid
import subprocess
from pathlib import Path
import pandas as pd
from concurrent.futures import ProcessPoolExecutor, as_completed
import time
import argparse

key_pattern = re.compile(r"^Key:\s*(.+)")
corrections_pattern = re.compile(r"^Corrections:\s*\[(.*?)\]")


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

                    assert len(corrections) == 1

                    if key and cpp_lines:
                        if key in SELECTED_KEYS: 
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
        body_mod = re.sub(r'\bD\s*=', f"Ds[{j}] =", body_mod)

        first, second, _ = body_mod.split(sep="[SPLIT]")
        first, second = first.strip(), second.strip()

        out_lines.append(f"        /* Key: {key} */ \\\n        {first}\n")
        out_lines_temp.append(f"        /* Key: {key} */ \\\n        {second}\n")

        count += 1
        if not (count % BATCH_SIZE):
            # Write to file
            file_text = "".join(out_lines).rstrip("\\\n") + "\n\n\n" + "".join(out_lines_temp).rstrip("\\\n") + "\n\n"
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
        file_text = "".join(out_lines).rstrip("\\\n") + "\n\n\n" + "".join(out_lines_temp).rstrip("\\\n") + "\n\n"
        yield (results, file_text)

def process_chunk(chunk):
    results, text = chunk
    
    # Write header with coefficients
    id = uuid.uuid4()
    header_path = COEFFS_DIR / f"coeffs_{id}.h"
    header_path.write_text(text)

    executable = COEFFS_DIR / f"exec_{id}"

    try: 
        compile_cmd = ["g++","-std=c++17", f"-include{header_path}", str(CPP_FILE),"-o", str(executable),"-O2", f"-DBATCH_SIZE={BATCH_SIZE}", "-DEXT_COEFFS=1", f"-DN={N}"]
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
            numbers = [int(i) for i in output_str.split()]
            assert len(numbers) == min(BATCH_SIZE, len(results))
            curr_results = { k : [numbers[j]+corrs[0]] for j, (k, corrs) in results.items()} 

        except ValueError:
            raise ValueError(f"Result could not be parsed correctly: {output_str[:50]}")

    finally: 
        # Clean up generated files
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
                        print(f"Progress: {progress*BATCH_SIZE}")
                        # breaker = True

                    if progress % 1 == 0:
                        df = pd.DataFrame.from_dict(all_results, orient="index")
                        df.to_csv(CSV_FILE,  mode="a", index=True)
                        all_results = {}
                        print(f"Total Time Taken: {time.time()-start_time}")
                        breaker = True

                        
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
    results = []
    for f in done_futures:
        try: 
            results.append(f.result())
        except Exception as e:
            print(e)
    return results, remaining_futures

def parse_args():
    parser = argparse.ArgumentParser(
        description="Run the job with configurable workers, N, and batch size."
    )

    parser.add_argument(
        "--max_workers",
        type=int,
        required=True,
        help="Number of workers to use in the pool"
    )

    parser.add_argument(
        "--N",
        type=int,
        required=True,
        help="The value of N to process"
    )

    parser.add_argument(
        "--batch_size",
        type=int,
        required=True,
        help="Number of items per batch"
    )

    return parser.parse_args()



BASE_PATH = Path("Dataset/zeta_functions")
COEFFS_DIR = Path("CppLib/Coeffs")
CONST_HEADER = Path("CppLib/constants.h")

args = parse_args()
MAX_WORKERS = args.max_workers
BATCH_SIZE = args.batch_size
N = args.N


# For surfaces without a line
# CPP_FILE = Path("CppLib/count_cubic_bigq.cpp")
# CSV_FILE = BASE_PATH / f"point_counts_cubic_bigq_{N}.csv"
# DATA_FILES = [f"Dataset/cpp_coeffs/cubic_coeffs_bigq_{N}.txt" ]

# For surfaces with a line
CPP_FILE = Path("CppLib/count_quad_bigq.cpp")
CSV_FILE = BASE_PATH / f"point_counts_quad_bigq_{N}_test.csv"
DATA_FILES = [f"Dataset/cpp_coeffs/quad_coeffs_bigq_{N}.txt" ]


SELECTED_KEYS = set(map(str.strip, list(open(BASE_PATH/ f"zeta_fails_{N}.txt").readlines())))
print(f"Computing values for N={N} and {len(SELECTED_KEYS)} surfaces")

if os.path.exists(CSV_FILE):
    SELECTED_KEYS -= set([line.split(",")[0] for line in open(CSV_FILE).readlines()])

print(f"Computing values for N={N} and {len(SELECTED_KEYS)} surfaces")


if __name__ == "__main__":
    process_chunks_in_batches(stream_batched_chunks(), max_workers=MAX_WORKERS)