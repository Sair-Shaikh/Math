#!/bin/bash

SRC="CppLib/generate_bigFq_tables.cpp"

for i in {1..20}; do
    echo "----------------------------------------"
    echo "Compiling for N=$i ..."

    exe="bigFq_$i"

    # Compile with macros ARM and N=i
    g++ -O3 -std=c++17 -DN=$i "$SRC" -o "$exe"
    if [ $? -ne 0 ]; then
        echo "Compilation failed for N=$i"
        exit 1
    fi

    echo "Running $exe ..."
    ./"$exe"

    echo "Deleting $exe ..."
    rm "$exe"
done