#!/bin/bash

# Loop over all folders scenario_1 ... scenario_18
for i in {1..18}; do
    echo "Running scenario_$i/run.sh..."
    # Enter the folder
    cd scenario_$i
    # Execute the script
    ./run.sh
    # Return to the main folder
    cd ..
    echo "Scenario_$i completed."
done

echo "All scenarios have been executed."
