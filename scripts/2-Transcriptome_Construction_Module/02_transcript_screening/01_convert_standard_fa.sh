#!/bin/bash

if [ $# -ne 2 ]; then
    echo "Usage: $0 <input_file> <output_file>"
    exit 1
fi

input_file=$1
output_file=$2

echo -n "" > $output_file

sequence_id=""
sequence=""

while IFS= read -r line
do
    if [[ $line == \>* ]]; then
        if [[ -n $sequence ]]; then
            echo -e "$sequence_id\n$sequence" >> $output_file
        fi

        sequence_id=$line
        sequence=""
    else
        sequence+=$line
    fi
done < "$input_file"

if [[ -n $sequence ]]; then
    echo -e "$sequence_id\n$sequence" >> $output_file
fi