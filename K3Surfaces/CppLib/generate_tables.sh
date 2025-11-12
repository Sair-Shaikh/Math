
#!/bin/bash
set -e  # exit on any error

for N in {1..11}; do
    echo "=============================="
    echo "Building and running for N = $N"
    echo "=============================="

    # Compile
    g++ -std=c++17 -O2 generate_Fq_tables.cpp -DN=$N -o generate_tables

    # Run
    ./generate_tables
    echo "Finished N = $N"
    echo
done