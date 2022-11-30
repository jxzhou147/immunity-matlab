#! /bin/bash
for file in `ls cmil10_cmccl2_1_il10low`
    do
        sed -i '12s/2 1 0/1 1 0/' cmil10_cmccl2_1_il10low/$file
    done