#!/bin/bash

# 固定参数
N=15
P=0.2
MAX_DIM=500

# 扫描的 depth 列表（你可以自己加）
DEPTHS=(15 20 25 30 35 40)

# 扫描的 sampleNum 列表
SAMPLES=(1000000)

# 每个条件重复次数
REPEATS=10000

# 输出目录
OUTDIR="estimate_third"
mkdir -p $OUTDIR

echo "Starting experiments..."
echo ""

exp_id=0

for DEPTH in "${DEPTHS[@]}"
do
    echo "=== Running depth = $DEPTH ==="

    for S in "${SAMPLES[@]}"
    do
        echo "  sampleNum = $S"

        for ((i=1; i<=REPEATS; i++))
        do
            ./main $N $DEPTH $S $MAX_DIM $P $exp_id 
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
