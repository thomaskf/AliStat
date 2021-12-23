# AliStat
Quantifying the Completeness of Multiple Sequence Alignments used in Phylogenetic and Phylogenomic Research

## Installation
The software is written in C++, and has been tested under the Linux and MacOS platforms. You need to have C++
compiler installed in your computer to compile the source code. The compilation steps are shown as follows:
```
$ tar -zxvf AliStat_1.xx.tar.gz
$ cd AliStat_1.xx
$ make
```
Then an executable file named *alistat* will appear.

## Metrics
Consider a multiple sequence alignment with *m* sequences and *n* sites. Given this alignment, we may compute:

- C<sub>a</sub> - Completeness score for the alignment
```
Ca = Number of unambiguous characters in the alignment / (m x n)
```

- C<sub>r</sub> - Completeness score for individual sequences
```
Cr = Number of unambiguous characters in a row (sequence) of the alignment / n
```

- C<sub>c</sub> - Completeness score for individual sites
```
Cc = Number of unambiguous characters in a column (site) of the alignment / m
```

- C<sub>ij</sub> - Completeness score for pairs of sequences
```
Cij = Number of columns with homologous pairs of unambiguous characters in sequences i and j / n
By definition, Cii = 1.
```

- I<sub>ij</sub> - Incompleteness score for pairs of sequences
```
Iij = 1 - Cij
```

- P<sub>ij</sub> - P-distance for pairs of sequences.
```
Let Wij be the set of columns with unambiguous homologous characters of sequences i and j.
Pij = (Number of sites with differences between sequences i and j across the columns in Wij) / |Wij|
Pij is undefined if |Wij| = 0, and Pii = 0.
```

## Usage of AliStat

Syntax:
```
./alistat <alignment file> <data type> [other options]
./alistat -h
```

```
  <alignment file> : Multiple alignment file in FASTA format
  
  <data type>      : 1  - Single nucleotides (SN);
                     2  - Di-nucleotides (DN);
                     3  - Codons (CD);
                     4  - 10-state genotype data (10GT);
                     5  - 14-state genotype data (14GT);
                     6  - Amino acids (AA);
                     7  - Mixture of nucleotides and amino acids (NA)
                          (User has to specify the data type for each partition
                           inside the partition file)
```

```
other options:

  -b               : Report the brief summary of figures to the screen
                     Output format:
                     File,#seqs,#sites,Ca,Cr_max,Cr_min,Cc_max,Cc_min,Cij_max,Cij_min
                     (this option cannot work with the options: -o,-t,-r,-m,-i,-d)
  -c <coding_type> : 0  - A, C, G, T [default];
                     1  - C, T, R (i.e. C, T, AG);
                     2  - A, G, Y (i.e. A, G, CT);
                     3  - A, T, S (i.e. A, T, CG);
                     4  - C, G, W (i.e. C, G, AT);
                     5  - A, C, K (i.e. A, C, GT);
                     6  - G, T, M (i.e. G, T, AC);
                     7  - K, M    (i.e. GT,   AC);
                     8  - S, W    (i.e. GC,   AT);
                     9  - R, Y    (i.e. AG,   CT);
                     10 - A, B    (i.e. A,   CGT);
                     11 - C, D    (i.e. C,   AGT);
                     12 - G, H    (i.e. G,   ACT);
                     13 - T, V    (i.e. T,   ACG);
                     (this option is only valid for <data type> = 1)
  -o <FILE>        : Prefix for output files
                     (default: <alignment file> w/o .ext)
  -n <FILE>        : Only consider sequences with names listed in FILE
  -p <FILE>        : Specify the partitions
                     For <data type> = 1 - 6, partition file format:
                       "<partition name 1>=<start pos>-<end pos>, ..."
                       "<partition name 2>= ..."
                       Example:
                           part1=1-50,60-100
                           part2=101-200
                       (enumeration starts with 1)
                     For <data type> = 7, partition file format:
                       "<SN/DN/CD/10GT/14GT/AA>, <partition name 1>=<start pos>-<end pos>, ..."
                       "<SN/DN/CD/10GT/14GT/AA>, <partition name 2>= ..."
                       Example:
                           SN,part1=1-50,60-100
                           AA,part2=101-200
                           CD,part3=201-231
                     (this option cannot be used with '-s' at the same time)
  -s <n1,n2>       : Sliding window analysis: window size = n1; step size = n2
                     (this option cannot be used with '-p' at the same time)
  -t <n3,n4,...>   : Only output the tables n3, n4, ...
                     1 - C scores for individual sequences (Cr)
                     2 - C scores for individual sites (Cc)
                     3 - Distribution of C scores for individual sites (Cc)
                     4 - Matrix with C scores for pairs of sequences (Cij)
                     5 - Matrix with incompleteness scores for pairs of
                         sequences (Iij = 1 - Cij)
                     6 - Table with C score and incompleteness scores for pairs
                         of sequences (Cij & Iij)
                     (default: the program does not output any tables)
                     If "-t" option is used but no <n3,n4,...>, then the program
                     outputs all tables
  -r <row|col|both>: Reorder the rows/columns (or both) of the alignment
                     according to the Cr/Cc scores
                     All the tables are displayed according to the reordered
                     alignment
                     To output the reordered alignment, please also use the
                     option -m
  -m <n5>          : Mask the alignment; 0 <= n5 <= 1
                     Output (1) the alignment with columns Cc >= n5 in the file
                     'Mask.fst', (2) the alignment with columns Cc < n5 in the
                     file 'Disc.fst', and (3) the alignment with an extra row in
                     the first line to indicate whether the column is masked in
                     the file 'Stat.fst'.
                     (Special case: if no <n5>, whole alignment is outputted in
                      the file 'Mask.fst'; if <n5> is 0, the alignment with
                      columns Cc > 0 is outputted in the file 'Mask.fst')
  -i <1|2|3>       : Generate heat map image for Cij scores of sequence pairs
                     1 - Triangular heat map
                     2 - Rectangular heat map
                     3 - Both
                     (if no number, then both triangular & rectangular heat map
                      files are outputted)
  -d               : Report the p distances between sequences
                     in the file with extension '.p-dist.csv'
                     (default: disabled)
                     Note: computation of p distances may take long time for
                          large number of sequences
  -u               : Color scheme of the heatmaps
                     1 - Default color scheme, suitable for color-blind persons
                     2 - Another color scheme
  -h               : This help page
```

Please refer to the enclosed user manual for detailed information.

## Data Types and Alphabets

Seven types of sequence data are considered: single nucleotides, di-nucleotides, codons, 10-state genotypes, 14-state genotypes, amino acids, and a mix of single nucleotides and amino acids.

The single nucleotides, di-nucleotides, codons, and mixture of nucleotides and amino acids use the nomenclature outlined above.

The 10- and 14-state genotype data are experimental in the sense that they have not yet been recognised by the wider scientific community as sound sources of genetic information pertaining to genomes of diploid and triploid species. Simulation-based research is currently under way to examine the usefulness of these types of data. How these types of data are generated may be revealed in due course.

When the sequences are regarded as

- single nucleotides, the data represent a 4-state alphabet (i.e., it comprises four letters: A, C, G, and T).
- di-nucleotides, the data represent a 16-state alphabet (i.e., it comprises 16 pairs of letters: AA, …, TT).
- codons, the data represent a 64-state alphabet (i.e., it comprises 64 triplets of letters: AAA, …, TTT).
- amino acids, the data represent a 20-state alphabet (i.e., it comprises 20 letters: A, …, Y).
- 10-state genotype data, the data represent a 10-state alphabet: A, C, G, T, R, Y, K, M, S, and W.
- 14-state genotype data, the data represent a 14-state alphabet: A, C, G, T, R, Y, K, M, S, W, B, D, H, and V.

## Contacts

- Dr. Lars Jermin (Email: lars.jermiin@anu.edu.au)
- Dr. Thomas Wong (Email: thomas.wong@anu.edu.au)

## Citation

Thomas K F Wong, Subha Kalyaanamoorthy, Karen Meusemann, David K Yeates, Bernhard Misof, Lars S Jermiin, **A minimum reporting standard for multiple sequence alignments**, *NAR Genomics and Bioinformatics*, Volume 2, Issue 2, June 2020, lqaa024

DOI: https://doi.org/10.1093/nargab/lqaa024

