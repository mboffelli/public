# RNA Length Bias Analysis in Transciptomics

This collection of Jupyter notebooks provides a comprehensive analysis of various mouse and human datasets, examining unspliced and spliced RNA. We explore metrics like CPM vs. length, TPM vs. length, dropout rates, among others.

Microarrays, 10x, bulk RNA-seq, CEL-Seq, SPLiT-seq, Smart-seq3, are different technologies used in transcriptomics, the study of RNA molecules within a cell. Each of these methods has unique features and applications, particularly in how they sequence and analyze RNA. Understanding RNA length-dependence across these technologies is crucial for several reasons, including data normalization, comparison of datasets, and improving overall data quality and interpretation.

**Microarrays**: Microarrays are a gene expression analysis technique where DNA is affixed to a solid surface like a glass slide, allowing for the examination of thousands of genes simultaneously. They measure gene expression through cDNA hybridization to DNA spots on the array, quantified by signal intensity. However, this method can introduce biases related to RNA length, hybridization efficiency, and background noise.

**10x Genomics**: 10x is  amethod to perform large-scale single-cell RNA (scRNA) sequencing. Analyzing RNA length-dependence is crucial due to its high-throughput nature and the potential introduction of biases during the capture and amplification of RNA from a large number of cells.

**Bulk RNA-seq**: Bulk RNA-seq is a method of sequencing where total RNA is extracted from a sample containing many cells and sequenced en masse, without distinguishing individual cell contributions. It provides a holistic view of gene expression, making RNA length-dependence analysis crucial for overall expression profiling.

**Celseq**: Celseq (Cell Expression by Linear amplification and Sequencing) is another technique for scRNA sequencing. It uses linear amplification rather than exponential, which can reduce certain biases associated with PCR amplification. Investigating RNA length-dependence in Celseq can provide insights into how linear amplification impacts the representation of different length transcripts compared to other methods.

**Split-seq**: Split-seq (Split-Pool Ligation-based Transcriptome sequencing) is also a method for scRNA sequencing. It labels the RNA from individual cells with unique molecular identifiers (UMIs) during a series of split-pool barcoding steps. This process allows for the sequencing of RNA from thousands of cells in a single experiment. Analyzing RNA length-dependence in split-seq can help in understanding how the barcoding and amplification steps affect different RNA molecules, especially considering the potential biases introduced during these steps.

**Smart-seq3**: Smart-seq3 is another scRNA sequencing tool, offering higher sensitivity and detection of genes per cell compared to its predecessors. Studying RNA length-dependence in smart-seq3 data is important because the method aims to capture full-length transcripts. Analysis can reveal how well it performs across a range of transcript lengths.


## Notebooks Overview

1. `microarrays.ipynb` - Examination of [specific dataset] using microarrays.
2. `splitseq.ipynb` - Analysis of [specific dataset] using split-seq method.
3. `smartseq3.ipynb` - Detailed exploration of smart-seq3 data from [specific dataset].
4. `10x.ipynb` - Insights from 10x Genomics data analysis with [specific dataset].
5. `bulk.ipynb` - Bulk RNA-seq data examination using [specific dataset].
6. `celseq.ipynb` - Celseq method applied to [specific dataset].

Each notebook, with its histograms and scatterplots, illuminates various aspects of RNA sequencing data. Additionally, your analysis quantifies length dependence with microarray data, an important aspect for comprehensive understanding.

The conclusion that all methods and technologies exhibit some level of length bias is significant. It not only underscores the inherent challenges in RNA sequencing but also highlights the necessity for continual methodological improvements and careful data interpretation.
