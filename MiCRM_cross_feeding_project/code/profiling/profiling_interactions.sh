#!/bin/bash

for i in 16; #2 4 8 16 32 64
    do

        python3 -m cProfile -s cumtime -o prof_results/profiling_comms_NO_scs$i get_interactions_profiling_NO_scs.py --size $i

        wait

    done
