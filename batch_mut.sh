#!/bin/bash

search_dir = `ls *.fasta`
for file in $search_dir
do
    python pdcli.py --mode dna -s $file -m sub -t C -r A -p 25 "${file}-primer.fasta"
done