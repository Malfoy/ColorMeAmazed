# ColorMeAmazed

## Description
ColorMeAmazed is a tool able to find shared sequences from a large amount of FASTA file (zipped or not)

## Usage


```
./cma -p . -c 10 -o output.fa
```

Here ColorMeAmazed will 

-Read all fasta files from the current directory "." indicated by the option "-p"

-Output the sequences shared by more than 10 files indicated by the option "-c"

-Write the shared sequences in output.fa indicated by the option "-o" 


## Implementation details

ColorMeAmazed will consider any file containing ".fa" or ".fasta" in its filename as a FASTA file

Each line will be considered a single sequence regarless of its length

Thus seq1 and seq2 will not be matched
```
>seq1
ATC
>seq2
ATCG
```
## Algorithm overview

### First step

ColorMeAmazed will insert all sequences from the files list in a Count Min Sketch (from https://github.com/barrust/count-min-sketch)

### Second step
ColorMeAmazed will index all sequence seen more than C times in the files (according to the Count Min SKetch wich can overestimate the count)

Each sequence is associated with the list of filename containing it

### Third step
ColorMeAmazed will output all sequences seen more than C times in the files (without overestimation) along with the list of filename containing it

### C parameter

The C parameter is a critical  as the shared sequences are stored in memory in an expensive manner.

A low C parameter can index a large amount of sequences and need a very high amount of memory
## Sample output
```
>./file1.fa ./file2.da
CAGCTAGCTGATCGTACGTGCATGACTGACTGATCGTCA
```

Mean that CAGCTAGCTGATCGTACGTGCATGACTGACTGATCGTCA is found in file1.fa file2.fa

## Potential future features (aka todo list)
Multithreading

Handle FASTQ/truncated FASTA

Handle reverse-complement

Use file of file input instead of folder

