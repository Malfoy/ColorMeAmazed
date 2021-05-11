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

## Algorithm overview

### First step

ColorMeAmazed will insert all sequences from the files list in a Count Min Sketch (from https://github.com/barrust/count-min-sketch)

### Second step
ColorMeAmazed will index all sequence seen more than C times in the files (according to the Count Min SKetch wich can overestimate the count)

Each sequence is associated with the list of filename containing it

### Third step
ColorMeAmazed will output all sequences seen more than C times in the files (without overestimation) along with the list of filename containing it

## Sample output
```
>./file1.fa ./file2.da
CAGCTAGCTGATCGTACGTGCATGACTGACTGATCGTCA
```

Mean that CAGCTAGCTGATCGTACGTGCATGACTGACTGATCGTCA is found in file1.fa file2.fa
