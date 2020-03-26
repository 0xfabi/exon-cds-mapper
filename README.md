# Leaderless transcription in eukaryotes
![Overview](./overview.jpg)

Code is written in Python 3.7.4. Required python packages are:
* pandas==1.0.1

# Task 1
For each gene in a given `gff3` file validate if first exon and first CDS have same starting position.
If they have the same starting position, they are returned as list in an output file.
Keep in mind that each gene is located on upper or lower strand of the DNA sequence. For the lower strand the stop position is the wanted start position.

## corresponding:
code: src/map_first_exon_cds.py

input file: gff file (version 3) is provided by course organizers

example entry:
```text
##gff-version   3
Chr1	Araport11	gene	3631	5899	.	+	.	ID=AT1G01010;Name=AT1G01010;Note=NAC domain containing protein 1;symbol=NAC001;some more explanantion ...additional information
```

output file: output_data/matching_exon_cds.tsv

output file is tab seperated and contains following columns:
```text
GeneIdentifier	Start	Strand
```

# Task 2
Validate if mapped reads overlap position of first exon in a way they start before the first exon position.

## corresponding:
code: TODO
file: bam file is provided by course organizers

example entry:
```text
HWI-1KL168:41:H7GWJADXX:2:2113:20686:32469	0	0	3630	255	129M1S	-1	-1	129	AAATTATTAGATATACCAAACCAGAGAAAACAAATACATAATCGGAGAAATACAGATTACAGAGAGCGAGAGAGATCGACGGCGAAGCTCTTTACCCGGAAACCATTGAAATCGGACGGTTTAGTGAAAN	array('B', [30, 31, 31, 26, 33, 33, 35, 35, 27, 34, 32, 39, 37, 40, 33, 38, 38, 38, 39, 37, 36, 37, 35, 37, 37, 36, 39, 34, 38, 36, 37, 32, 32, 35, 35, 38, 39, 39, 31, 29, 30, 30, 27, 35, 31, 30, 37, 39, 40, 40, 23, 30, 19, 34, 32, 34, 37, 38, 31, 34, 19, 28, 28, 31, 38, 36, 32, 36, 32, 27, 28, 24, 32, 28, 32, 26, 34, 33, 33, 29, 33, 25, 33, 33, 29, 27, 27, 27, 29, 32, 25, 32, 34, 34, 34, 34, 34, 34, 24, 20, 27, 20, 30, 32, 33, 34, 34, 34, 18, 31, 32, 29, 31, 30, 33, 27, 33, 33, 33, 20, 27, 33, 34, 34, 19, 31, 29, 25, 18, 29])	[('NH', 1), ('HI', 1), ('AS', 127), ('nM', 0)]

```
