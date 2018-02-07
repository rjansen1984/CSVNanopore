# CSVNanopore
Reads the output CSV files from the EPI2ME 16S workflow and calculates the species and genus occurance based on the NCBI taxonomy IDs.
To run the script: python csv_nanopore.py
The script asks for the location of the CSV files and which rank to use (1: phylum ----> 2: class ----> 3: order ----> 4: family ----> 5: genus ----> 6: species). Two new files will be created. One with the rank names and ocurrence count and one with the rank names occurrence count and percentages.
