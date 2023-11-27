# RNA Length Bias Analysis in Transciptomics

This collection of Jupyter notebooks provides a comprehensive analysis of various mouse and human datasets, examining unspliced and spliced RNA. We explore metrics like CPM vs. length, TPM vs. length, dropout rates, among others.

10x, bulk RNA-seq, CEL-Seq, SPLiT-seq, Smart-seq3, are different technologies used in transcriptomics, the study of RNA molecules within a cell. Each of these methods has unique features and applications, particularly in how they sequence and analyze RNA. Understanding RNA length-dependence across these technologies is crucial for several reasons, including data normalization, comparison of datasets, and improving overall data quality and interpretation.

**10x Genomics**: 10x is  amethod to perform large-scale single-cell RNA (scRNA) sequencing. Analyzing RNA length-dependence is crucial due to its high-throughput nature and the potential introduction of biases during the capture and amplification of RNA from a large number of cells.

**Bulk RNA-seq**: Bulk RNA-seq is a method of sequencing where total RNA is extracted from a sample containing many cells and sequenced en masse, without distinguishing individual cell contributions. It provides an overall picture of gene expression in the sample. Analyzing RNA length-dependence in bulk RNA-seq is important for understanding overall gene expression profiles and how they might be influenced by transcript length.

**Celseq**: Celseq (Cell Expression by Linear amplification and Sequencing) is another technique for scRNA sequencing. It uses linear amplification rather than exponential, which can reduce certain biases associated with PCR amplification. Investigating RNA length-dependence in Celseq can provide insights into how linear amplification impacts the representation of different length transcripts compared to other methods.

**Split-seq**: Split-seq (Split-Pool Ligation-based Transcriptome sequencing) is also a method for scRNA sequencing. It labels the RNA from individual cells with unique molecular identifiers (UMIs) during a series of split-pool barcoding steps. This process allows for the sequencing of RNA from thousands of cells in a single experiment. Analyzing RNA length-dependence in Split-seq can help in understanding how the barcoding and amplification steps affect different RNA molecules, especially considering the potential biases introduced during these steps.

**Smart-seq3**: Smart-seq3 is another scRNA sequencing tool, offering higher sensitivity and detection of genes per cell compared to its predecessors. Studying RNA length-dependence in Smart-seq3 data is important because the method aims to capture full-length transcripts. Analysis can reveal how well it performs across a range of transcript lengths.

## Why Conduct Length-Dependence Analysis
Conducting an analysis on RNA length-dependence in transcriptomics across these technologies is vital for several reasons:

## Notebooks Overview

1. `splitseq.ipynb` - Analysis of [specific dataset] using Split-seq method.
2. `smartseq3.ipynb` - Detailed exploration of Smart-seq3 data.
3. `10x.ipynb` - Insights from 10x Genomics data analysis.
4. `bulk.ipynb` - Bulk RNA-seq data examination.
5. `celseq.ipynb` - Celseq method applied to [specific context].

Each notebook contains histograms and scatterplots for the mentioned metrics, providing valuable insights into RNA sequencing data.
