#!/bin/bash

# Run the commands and generate diff files
run_ads_programm() {
    # Remove any existing output files
    rm -f diff*.txt
    
    differences_found=false
    
    for ((i=1; i<=20; i++))
    do
        echo "Iteration $i"
        
        python3 gen.py > input.txt

        # Run the first command
        ./ads_programm rmq input.txt sparse.txt
        
        # Run the second command
        ./ads_programm rmq input.txt naive.txt
        
        # Generate diff file between sparse.txt and naive.txt
        diff_output=$(diff sparse.txt naive.txt)
        
        if [[ -n "$diff_output" ]]; then
            differences_found=true
            echo "Difference found in iteration $i"
        fi
        
        echo "$diff_output" > "diff${i}.txt"
    done
    
    if [[ $differences_found = true ]]; then
        echo "Differences found in at least one iteration."
    else
        echo "No differences found in any iteration."
    fi
}

# Call the function to execute the commands
run_ads_programm