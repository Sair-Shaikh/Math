import os
import re
import uuid
import subprocess
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor, as_completed
import pandas as pd

# === Paths ===
BASE_DIR = Path("K3Surfaces/CppLib")
DATA_FILE = Path("K3Surfaces/Dataset/container_quad_coeffs_table.txt")
HEADER_PATH = Path("K3Surfaces/CppLib/coeffs.h")
CPP_FILE = Path("K3Surfaces/CppLib/point_count_containers.cpp")
EXECUTABLE = Path("K3Surfaces/CppLib/point_count_containers.out")
FQ_HEADER = Path("K3Surfaces/CppLib/Fq.h")
CONST_HEADER = Path("K3Surfaces/CppLib/constants.h")
FQ_GCH = FQ_HEADER.with_suffix(".gch")
CONST_GCH = CONST_HEADER.with_suffix(".gch")

# === Regex patterns for parsing ===
key_pattern = re.compile(r"^Key:\s*(.+)")
corrections_pattern = re.compile(r"^Corrections:\s*\[(.*?)\]")

# === Stream Chunks ===
def stream_chunks():
    results = {}

    with open(DATA_FILE, "r") as f:
        key = None
        cpp_lines = []
        in_cpp_block = False

        for line in f:

            key_match = key_pattern.match(line)
            corr_match = corrections_pattern.match(line)

            if key_match:
                key = key_match.group(1).strip()
                in_cpp_block = True
                continue

            elif corr_match:
                # End of chunk
                corrections = [int(x) for x in re.findall(r"\d+", corr_match.group(1))]
                in_cpp_block = False

                if key and cpp_lines:
                    yield (key, cpp_lines, corrections)
                    key, cpp_lines = None, []
                    cpp_lines = []
                continue

            elif in_cpp_block:
                cpp_lines.append(line)
            else:
                raise NotImplementedError

def process_chunk(chunk):
    key, cpp_header, corrections = chunk
    uncorrected_counts = compile_and_run(cpp_header)
    result = {k : v+corrections[k-1] for k, v in uncorrected_counts.items()}
    return key, result

def process_chunks_in_batches(chunk_generator, max_workers=None, batch_size=None):
    if max_workers is None:
        max_workers = os.cpu_count()
    if batch_size is None:
        batch_size = max_workers  # submit one batch per core by default

    results = {}
    progress = 0
    batch = []

    with ProcessPoolExecutor(max_workers=max_workers) as executor:
        futures = set()

        breaker = False
        for chunk in chunk_generator:
            # Submit new task
            futures.add(executor.submit(process_chunk, chunk))
            batch.append(chunk)

            # If batch is full, wait for some futures to complete
            if len(futures) >= batch_size:
                done, futures = wait_some(futures)
                for key, result in done:
                    results[key] = result
                    progress += 1
                    if progress % 100 == 0:
                        print(f"Progress: {progress}")
                        breaker = True

                    if progress % 50 == 0:
                        df = pd.DataFrame.from_dict(results, orient="index")
                        df.to_csv(BASE_DIR / "partial_progress.csv",  mode="a")
                        results = {}

                batch.clear()

            if breaker: 
                break
                        
        # Process remaining futures
        while futures:
            done, futures = wait_some(futures)
            for key, result in done:
                results[key] = result
                progress += 1
                df = pd.DataFrame.from_dict(results, orient="index")
                df.to_csv(BASE_DIR / "partial_progress.csv", mode="a")
                results = {}


    return results

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

# === For each chunk: write, compile, run, collect, clean ===
def compile_and_run(cpp_lines):

    # Write coeffs.h
    id = uuid.uuid4()
    header_path = BASE_DIR / f"coeffs_{id}.h"
    header_path.write_text("".join(cpp_lines))

    executable = BASE_DIR / f"exec_{id}"

    results = {}

    try: 
        for i in range(8,12):
            # FQ_HEADERN = BASE_DIR / f"Fq_N{i}.h"
            # compile_cmd = ["g++","-std=c++17", f"-include{FQ_HEADERN}",f"-include{header_path}", str(CPP_FILE),"-o", str(executable),"-O2"]
            compile_cmd = ["g++","-std=c++17",f"-include{header_path}", str(CPP_FILE),"-o", str(executable),"-O2",f"-DN={i}"]
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
                value = int(output_str)
            except ValueError:
                print(f"Warning: Non-integer output: N={i}: {output_str}")
                value = None
            results[i] = value
    finally: 
        # Clean up generated files
        if header_path.exists():
            header_path.unlink()
        if executable.exists():
            executable.unlink()

    return results

# === Run ===
if __name__ == "__main__":
    # for i in range(8,11): 
    #     # Compile constants.h and Fq.h with N plugged in
    #     FQ_GCH = BASE_DIR / f"Fq_N{i}.h.gch"
    #     subprocess.run(["g++","-std=c++17","-O2",f"-DN={i}","-x", "c++-header",str(FQ_HEADER),"-o", FQ_GCH], check=True)
    
    progress = 0
    results = process_chunks_in_batches(stream_chunks(), max_workers=os.cpu_count(), batch_size=8)

