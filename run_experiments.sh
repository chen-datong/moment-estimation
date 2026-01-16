#!/bin/bash

# 固定参数
N=15
MAX_DIM=500
S=1000000
H=1.0

# 扫描的 depth 列表（你可以自己加）
DEPTHS=(20 30)

# 扫描的 sampleNum 列表
PS=(0.0 0.05 0.10 0.15 0.25 0.30)

# 每个条件重复次数
REPEATS=1000

# 输出目录
OUTDIR="estimate_third"
mkdir -p $OUTDIR

echo "Starting experiments..."
echo ""

exp_id=0

for DEPTH in "${DEPTHS[@]}"
do
    echo "=== Running depth = $DEPTH ==="

    for P in "${PS[@]}"
    do
        echo "  p = $P"

        for ((i=1; i<=REPEATS; i++))
        do
            ./main $N $DEPTH $S $MAX_DIM $P $H $exp_id 
            ((exp_id++))

            if (( i % 500 == 0 )); then
                echo "    progress: $i / $REPEATS"
            fi
        done
    done
done

echo ""
echo "ALL EXPERIMENTS FINISHED!"
echo "Total experiments = $exp_id"
