#!/bin/bash

# Output file
outfile="sq_input_files.txt"

# Empty or create the output file
> "$outfile"

# Get sorted lists of each file type
class_files=($(ls -d "$PWD/"*classification.txt 2>/dev/null | sort))
junction_files=($(ls -d "$PWD/"*junctions.txt 2>/dev/null | sort))
tpm_files=($(ls -d "$PWD/"*TPM.tsv 2>/dev/null | sort))

# Find the maximum number of rows among the three lists
max_len=$(( ${#class_files[@]} > ${#junction_files[@]} ? ${#class_files[@]} : ${#junction_files[@]} ))
max_len=$(( $max_len > ${#tpm_files[@]} ? $max_len : ${#tpm_files[@]} ))

# Write rows
for ((i=0; i<max_len; i++)); do
    printf "%s\t%s\t%s\n" \
        "${class_files[i]:-}" \
        "${junction_files[i]:-}" \
        "${tpm_files[i]:-}" >> "$outfile"
done
