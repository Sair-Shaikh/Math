
for N in {18..20}; do
    echo "Running N=$N"
    uv run orchestrator_fast.py \
        --batch_size 30 \
        --max_workers 6 \
        --N $N
done