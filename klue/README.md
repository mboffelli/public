# klue: k-mer based local uniqueness exploration.

klue is a bioinformatics tool designed for the analysis of genomic data through k-mer processing. It excels in constructing and interpreting colored de Bruijn graphs, a powerful representation where nodes symbolize k-mers (subsequences of length k) of DNA sequences, and edges connect k-mers that overlap by k−1 nucleotides. The tool's specialty lies in accurately extracting contigs —contiguous sequences composed of unitigs— from genomic data stored in FASTA/FASTQ files.

A unitig is a non-branching path within the de Bruijn graph. It represents a sequence of overlapping k-mers that can be confidently assembled together, offering a continuous and unambiguous stretch of the genome. klue then assembles these unitigs into larger contigs, effectively reconstructing significant portions of the genome. The process extends to refining these contigs based on specific criteria, ensuring their biological relevance and accuracy. This careful assembly and refinement make klue a valuable tool in genomic research, particularly for uncovering complex genomic structures and variations.

## Key components

- `KmerIndex.cpp`: Implements k-mer indexing for graph construction and traversal.
- `main.cpp`: Provides the command-line interface, orchestrating the tool's functions and interactions.
- `ExpressionParser.cpp`: Parses set operation expressions for complex k-mer set operations.

## Getting started with klue

To use klue, follow the steps below:

### Step 1: Build a colored compact de Bruijn Graph (ccdBG) with Bifrost

Start by building a colored k-mer compact de Bruijn graph over long-read genome assemblies (FASTA/FASTQ files).

### Step 2: Extract unique contigs with klue

Once the graph is constructed, use klue to extract contigs that are unique to each strain or sample (default), extract contigs that are shared between `n-1` samples but not all `n`. klue can also extract contigs that are unique to a specified subset or combination of sets of strains or samples. klue selects contigs of length `2k-1`; this specific length ensures that each k-mer within a specific contig overlaps a putative single nucleotide polymorphism (SNP), which is vital for the subsequent mapping and identification process.

For example, using klue with `n=4` strains or samples, there are 2<sup>n</sup> - 1 = 15 distinct areas representing all unique intersections and combinations of the four sets. klue  can be used to extract contigs that are unique to each of these areas, allowing identification of sequences that are specific to individual strains, as well as those that are shared between any combination of 2, 3, or all 4 of the strains.

<img src="https://github.com/mboffelli/projects/blob/main/images/venn-4.png" width="400" height="400">

### Step 3: Map reads and assign strains

Use tools like kallisto to index the contigs and map your sequencing reads for strain assignment.

## Example workflow

1. Input a minimum of two FASTA/FASTQ files to klue to distinguish unique contigs to each file.
2. Align these unique contigs to a reference genome using methods such as kallisto or other suitable tools.
3. Analyze klue's classification to determine the identity or differences between the samples. The classification is be based on which strain has the most unique molecular identifiers (UMIs) assigned to it.

## Results

The results from KLUE (a FASTA file containing extracted contigs) will allow you to classify the cell identity based on the predominant strain indicated by the UMIs.

For example, when tasked with distinguishing and classifying cells from eight distinct and genetically varied inbred mouse strains, klue performed with an impressive precision that was on par with pairwise comparisons. The plots below highlight klue's robustness, showcasing that it maintains high accuracy in identifying unique contigs across all eight strains, surpassing a 10<sup>2</sup> UMI threshold. Such efficacy in multi-strain classification is not common in genomic analysis tools, marking a significant step forward in the field.

klue is able to not only extract unique contigs from a complex mixture of multiple strains but also to classify them with high reliability, thus proving the tool's advanced capability in managing and interpreting intricate genomic datasets.

<img src="https://github.com/mboffelli/public/blob/main/klue/results/c57bl6j.png" width="400" height="400">
<img src="https://github.com/mboffelli/public/blob/main/klue/results/129s1svimj.png" width="400" height="400">
<img src="https://github.com/mboffelli/public/blob/main/klue/results/nzohlltj.png" width="400" height="400">
<img src="https://github.com/mboffelli/public/blob/main/klue/results/casteij.png" width="400" height="400">
<img src="https://github.com/mboffelli/public/blob/main/klue/results/pwkphj.png" width="400" height="400">
<img src="https://github.com/mboffelli/public/blob/main/klue/results/aj.png" width="400" height="400">
<img src="https://github.com/mboffelli/public/blob/main/klue/results/nodshiltj.png" width="400" height="400">
<img src="https://github.com/mboffelli/public/blob/main/klue/results/wsbeij.png" width="400" height="400">



## Features

- **Graph Construction**: Builds colored de Bruijn graphs for unique and shared sequence identification.
- **Contig Reconstruction**: Produces extended contig sequences from reads by joining contiguous unitigs/nucleotides.
- **Command-Line Parsing**: Customizes operations via command-line arguments.

## Usage
- **Distinguish k-mers belonging exclusively to one color**: Extract unique contigs from input FASTA/FASTQ files:

  ```klue distinguish -k 31 -o output.fasta input1.fasta input2.fasta```
- **Distinguish k-mers not shared between all colors**: Extract unique contigs that are not found in all input files.
  
  ```klue distinguish --all-but-one -k 31 -o output.fasta input1.fasta input2.fasta```
- **Distinguish k-mers shared between every combination of subsets**:  Extract unique contigs that are found in each combination of input files (e.g., only A, only B, only C, only A and B, only A, B, and C, etc.)

  ```klue distinguish --combinations -k 31 -o output.fasta input1.fasta input2.fasta```
- **Distinguish k-mers shared between a combination of subsets**:  Extract contigs based on complex set operations (e.g., shared between A and B but not in C)
  
  ```klue distinguish -s "(AIB)\C" -k 31 -o output.fasta input1.fasta input2.fasta```
- **Extend extracted contigs**: Perform a forward depth-first search traversal to extend unique contigs from input FASTA/FASTQ files to get longer sequences.
  
  ```klue --extend -k 31 -o output.fasta input1.fasta input2.fasta```

