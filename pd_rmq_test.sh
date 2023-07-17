#!/bin/bash

# Compile the program
g++ main.cpp -o ads_programm -std=c++20 -O3 -g

if [ $? -ne 0 ]; then
    echo "Compilation failed. Exiting..."
    exit 1
fi

# Run the program three times
for ((k=1; k<=3; k++))
do
    ./ads_programm pd "predecessor_example_$k.txt" "pred_$k.txt"
done

for ((k=1; k<=2; k++))
do
    ./ads_programm rmq "rmq_example_$k.txt" "rmq_$k.txt"
done


# Check for differences using diff
differences_found=0

for ((k=1; k<=3; k++))
do
    if ! diff -q "test/pd_example_$k.txt" "pred_$k.txt" >/dev/null; then
        echo "Differences found between pd_example_$k.txt and pred_$k.txt"
        differences_found=1
    fi
done

for ((k=1; k<=2; k++))
do
    if ! diff -q "test/out$k.txt" "rmq_$k.txt" >/dev/null; then
        echo "Differences found between out$k.txt and rmq_$k.txt"
        differences_found=1
    fi
done

# Report if differences were found
if [ $differences_found -eq 1 ]; then
    echo "Differences were found."
else
    echo "No differences were found."
fi
