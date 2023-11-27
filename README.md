# Mayuko Boffelli
### Research repository
Welcome! This repository houses two primary projects: `klue` and `RNA length bias`. Both projects reflect my commitment to addressing challenges in computational genomics methods; the technology `klue` is still in development, but preliminary results can be found in https://github/mboffelli/projects/tree/main/klue/results.

## klue
Genomic research often hinges on using a "reference genome" to analyze gene expression and identify structural variants. However, the limitations of reference genomes, such as their incompleteness and variation across individuals and species, can hinder the discovery of critical genetic differences and unique transcripts crucial for cell identity. 

To address these limitations, we have developed the method k-mer based local uniqueness exploration (klue). klue stands out for its ability to detect de novo DNA and RNA sequences that are distinct between samples or divergent from reference genomes, without relying on pre-existing structures or references. klue addresses the challenges posed by incomplete reference genomes and the genetic diversity between individuals, strains, and species. By assembling sequences that are exclusive to a sample from a particular experimental group, klue offers a sensitive solution to uncover genomic irregularities, especially those undetectable by traditional reference-based approaches.

Implemented in C++ and utilizing Bifrost for constructing an index for k-mers in a colored and compacted de Bruijn graph, klue offers a sensitive solution to uncover genomic irregularities.

**Note**: Since klue is still in development, I have omitted source files in `/src` that are integral to klue's implementation, but to which I did not contribute. My contributions to the project may be found in `/src/KmerIndex.cpp` and `/src/ExpressionParser.cpp`.

#### More about klue
[https://github.com/mboffelli/public/tree/main/klue#readme]

## RNA-seq length bias in transcriptomics
The "RNA length bias" project investigates the impact of transcript length on gene expression measurements in RNA sequencing (RNA-seq). My analysis quantifies the dependence of RNA sequencing technologies on the length of these RNA molecules, both unspliced and spliced, across various human and mouse samples. This kind of analysis is crucial in understanding how different sequencing technologies might bias toward or against certain RNA lengths, which can influence data interpretation and comparisons across different methods and studies.

#### More about length bias
[https://github.com/mboffelli/public/tree/main/length-bias#readme]

## Acknowledgements
This repository reflects collaborative efforts and individual research undertakings. Key collaborators include MD/PhD student Delaney K. Sullivan, with whom I'm working on klue.
