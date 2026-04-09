#!/bin/bash

# --- Automatic Python Detector ---
if command -v python3 &>/dev/null; then
    PY_CMD="python3"
elif command -v python &>/dev/null; then
    PY_CMD="python"
else
    echo "Error: Python not found! You may need to run 'module load python' or 'module load anaconda'."
    exit 1
fi
# ---------------------------------

# ==============================================================================
# MODE 1: GENERATE TESTCASES (--gen)
# ==============================================================================
if [ "$1" == "--gen" ]; then
    echo "=========================================================="
    echo " Generating HARD Testcases and Checker Scripts using $PY_CMD"
    echo "=========================================================="
    
    mkdir -p testcases

    # 1. The Ultra-Hard Graph Generator
    cat << 'EOF' > generate_tests.py
import random
import os

random.seed(42)

configs = [
    ("01_tight_s",   300,  0.55,    1,  20,   5),
    ("02_tight_m",   400,  0.55,    1,  20,   5),
    ("03_dense_s",   300,  0.65,    1,  20,   6),
    ("04_dense_m",   400,  0.65,    1,  20,   6),
    ("05_ta_clone",  500,  0.60,    1,  20,   9),
    ("06_hard_500",  500,  0.65,    1,  20,   7),
    ("07_hard_600",  600,  0.60,    1,  20,   8),
    ("08_stress_a",  600,  0.65,    1,  20,   7),
    ("09_stress_b",  700,  0.60,    1,  20,   8),
    ("10_ultra",     500,  0.70,    1,  15,   7),
]

os.makedirs("testcases", exist_ok=True) # Works in Python 3.2+

for name, N, density, clo, chi, bmul in configs:
    filepath = "testcases/{}.txt".format(name)
    if os.path.exists(filepath):
        print("Skipping {} (already exists)".format(name))
        continue

    avg_cost = (clo + chi) / 2.0
    B = int(bmul * avg_cost)

    print("Generating {}: N={}, density={}, cost=[{},{}], B={} (~{} nodes fit)...".format(
        name, N, density, clo, chi, B, bmul))

    edges = []
    for i in range(N):
        for j in range(i + 1, N):
            if random.random() < density:
                edges.append((i, j))

    with open(filepath, "w") as f:
        f.write("{} {} {}\n".format(N, len(edges), B))
        for i in range(N):
            f.write("{} {}\n".format(random.randint(1, 100), random.randint(clo, chi)))
        for u, v in edges:
            f.write("{} {}\n".format(u, v))

    print("  -> {} edges written.".format(len(edges)))

print("Done generating hard testcases.")
EOF

    # 2. The Universal Verifier Script
    cat << 'EOF' > checker.py
import sys

def check(input_file, output_file, expected_profit):
    with open(input_file, 'r') as f:
        lines = [l.strip() for l in f.readlines() if l.strip()]
    N, E, B = map(int, lines[0].split())
    profits, costs = [], []
    for i in range(1, N + 1):
        p, c = map(int, lines[i].split())
        profits.append(p), costs.append(c)
    
    edges = set()
    for i in range(N + 1, N + 1 + E):
        u, v = map(int, lines[i].split())
        edges.add((u, v)), edges.add((v, u))

    with open(output_file, 'r') as f:
        out_lines = [l.strip() for l in f.readlines() if l.strip()]
    if not out_lines:
        print("FAIL: Output file is empty")
        sys.exit(1)
        
    reported_profit = int(out_lines[0])
    clique = list(map(int, out_lines[1].split())) if len(out_lines) > 1 and out_lines[1] else []

    if any(v < 0 or v >= N for v in clique):
        print("FAIL: Invalid vertex IDs")
        sys.exit(1)
    if not all(clique[i] < clique[i+1] for i in range(len(clique)-1)):
        print("FAIL: Vertices not sorted in ascending order")
        sys.exit(1)
    for i in range(len(clique)):
        for j in range(i + 1, len(clique)):
            if (clique[i], clique[j]) not in edges:
                print("FAIL: Not a clique! Missing edge ({}, {})".format(clique[i], clique[j]))
                sys.exit(1)
                
    actual_cost = sum(costs[v] for v in clique)
    if actual_cost > B:
        print("FAIL: Budget exceeded ({} > {})".format(actual_cost, B))
        sys.exit(1)
        
    actual_profit = sum(profits[v] for v in clique)
    if actual_profit != reported_profit:
        print("FAIL: Profit math mismatch (Printed: {}, Actual: {})".format(reported_profit, actual_profit))
        sys.exit(1)
        
    if actual_profit != expected_profit:
        print("FAIL: Suboptimal. Expected Pmax={}, Got={}".format(expected_profit, actual_profit))
        sys.exit(1)

    print("PASS")
    sys.exit(0)

check(sys.argv[1], sys.argv[2], int(sys.argv[3]))
EOF

    if [ -f "tainput.txt" ]; then
        cp tainput.txt testcases/00_tainput.txt
        echo "Copied tainput.txt to testcases/"
    fi

    # Use the detected python command
    $PY_CMD generate_tests.py
    exit 0
fi

# ==============================================================================
# PREREQUISITE CHECK
# ==============================================================================
if [ ! -d "testcases" ]; then
    echo "Error: 'testcases' directory not found. Please run './run_tests.sh --gen' first."
    exit 1
fi

LOG_FILE="baseline_logs.csv"

# ==============================================================================
# MODE 2: BASELINE RUN (--baseline)
# ==============================================================================
if [ "$1" == "--baseline" ]; then
    echo "=========================================================="
    echo " Building Baseline Logs (Sequential C++)"
    echo "=========================================================="
    
    module purge
    module load compiler/gcc/9.1.0
    
    g++ -O3 -std=c++17 seq.cpp -o seq_bin
    if [ $? -ne 0 ]; then echo "Sequential compile failed!"; exit 1; fi

    echo "TC_Name,Seq_Time_Sec,Optimal_Profit" > $LOG_FILE
    
    for tc in testcases/*.txt; do
        tc_name=$(basename "$tc")
        echo -n "Running $tc_name ... "
        
        start_time=$(date +%s.%N)
        ./seq_bin $tc tmp_out.txt
        end_time=$(date +%s.%N)
        
        seq_time=$(echo "$end_time - $start_time" | bc | awk '{printf "%.4f", $0}')
        profit=$(head -n 1 tmp_out.txt)
        
        echo "Done! (Time: ${seq_time}s, Profit: ${profit})"
        echo "$tc_name,$seq_time,$profit" >> $LOG_FILE
    done
    rm tmp_out.txt
    echo "Baseline saved to $LOG_FILE"
    exit 0
fi

# ==============================================================================
# MODE 3: MPI RUN (No Flags)
# ==============================================================================
if [ ! -f "$LOG_FILE" ]; then
    echo "Error: $LOG_FILE not found! Please run './run_tests.sh --baseline' first."
    exit 1
fi

if [ ! -f "checker.py" ]; then
    echo "Error: checker.py not found! Please run './run_tests.sh --gen' first."
    exit 1
fi

if [ "$1" != "" ]; then
    echo "Usage: ./run_tests.sh [--gen | --baseline]"
    exit 1
fi

echo "=========================================================="
echo " Compiling and Running MPI Version"
echo "=========================================================="

module purge
module load compiler/gcc/9.1.0
module load mpi/openmpi/4.1/gnu/mpivars

mpicxx -O3 -std=c++17 main.cpp -o mpi_bin
if [ $? -ne 0 ]; then echo "MPI compile failed!"; exit 1; fi

echo ""
printf "%-20s | %-5s | %-10s | %-10s | %-8s | %-8s\n" "Testcase" "Ranks" "Seq Time" "MPI Time" "Speedup" "Status"
printf "%s\n" "---------------------------------------------------------------------------------"

export UCX_LOG_LEVEL=error

for tc in testcases/*.txt; do
    tc_name=$(basename "$tc")
    
    baseline_info=$(grep "^${tc_name}," $LOG_FILE)
    seq_time=$(echo "$baseline_info" | cut -d',' -f2)
    optimal_profit=$(echo "$baseline_info" | cut -d',' -f3)
    
    for ranks in 2 4 8 12 16; do
        start_time=$(date +%s.%N)
        mpirun -np $ranks ./mpi_bin $tc mpi_out.txt 2>/dev/null
        end_time=$(date +%s.%N)
        
        mpi_time=$(echo "$end_time - $start_time" | bc | awk '{printf "%.4f", $0}')
        
        if (( $(echo "$mpi_time > 0" | bc -l) )); then
            speedup=$(echo "scale=2; $seq_time / $mpi_time" | bc)
        else
            speedup="INF"
        fi
        
        # Use detected python command for verification
        checker_out=$($PY_CMD checker.py $tc mpi_out.txt $optimal_profit)
        
        if [[ "$checker_out" == "PASS" ]]; then
            status="\e[32mPASS\e[0m"
        else
            status="\e[31mFAIL\e[0m"
        fi
        
        printf "%-20s | %-5s | %-10ss | %-10ss | %-8sx | %b\n" "$tc_name" "$ranks" "$seq_time" "$mpi_time" "$speedup" "$status"
        
        if [[ "$checker_out" != "PASS" ]]; then
            echo -e "   -> Reason: \e[31m$checker_out\e[0m"
        fi
    done
    printf "%s\n" "---------------------------------------------------------------------------------"
done

rm -f mpi_out.txt