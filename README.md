# HiCat
HiCat: A Semi-Supervised Approach for Cell Type Annotation

![Graphical Abstract](graphicAbstract.png)

This is the repo for the demo code of HiCat (Hybrid Cell Annotation using Transformative embeddings), a novel
semi-supervised pipeline for annotating cell types from single-cell RNA sequencing data. HiCat fuses
the strengths of supervised learning for known cell types with unsupervised learning to identify novel
types. This hybrid approach incorporates both reference and query genomic data for feature engineer-
ing, enhancing the embedding learning process, increasing the effective sample size for unsupervised
techniques, and improving the transferability of the supervised model trained on reference data when
applied to query datasets.