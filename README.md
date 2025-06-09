# HiCat: A Semi-Supervised Approach for Cell Type Annotation

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![arXiv](https://img.shields.io/badge/arXiv-2412.06805-b31b1b.svg)](https://arxiv.org/abs/2412.06805)

This repository contains the official demo code for **HiCat (Hybrid Cell Annotation using Transformative embeddings)**, a novel semi-supervised pipeline for annotating cell types in single-cell RNA sequencing (scRNA-seq) data.

Accurate cell type annotation is crucial for understanding cellular heterogeneity. However, supervised methods struggle to identify novel cell types, while unsupervised approaches often suffer from cluster impurity. HiCat was developed to overcome these limitations by integrating both approaches. It leverages reference (labeled) and query (unlabeled) data to enhance annotation accuracy for known cell types while simultaneously improving the discovery and differentiation of novel ones.

Our comprehensive benchmarking shows that HiCat consistently outperforms existing methods, particularly in its ability to distinguish and identify multiple novel cell types, including rare populations.

![HiCat Workflow](Figure1.png)
*Fig 1. The HiCat workflow, from data input and feature engineering to classification and final prediction.*

---

## 📖 Table of Contents

- [Key Features](#-key-features)
- [Workflow](#-workflow)
- [Performance Highlights](#-performance-highlights)
- [Installation](#-installation)
- [Usage](#-usage)
- [Data](#-data)
- [Citation](#-citation)
- [License](#-license)
- [Contact](#-contact)

---

## ✨ Key Features

HiCat introduces several key innovations to address the challenges of cell type annotation:

* **Multi-Resolution Feature Engineering**: HiCat creates a unique 53-dimensional feature space by combining batch-corrected principal components, UMAP embeddings, and unsupervised cluster identities. This provides a rich, condensed representation of cellular identity.
* **Fusion of Supervised & Unsupervised Signals**: It intelligently fuses supervised predictions from a CatBoost classifier with unsupervised clustering information from DBSCAN. This strategy ensures high accuracy for known cell types while enabling the identification of novel populations where the classifier has low confidence.
* **Discovery of Multiple Novel Cell Types**: A core strength of HiCat is its explicit capability to distinguish between multiple, different unseen cell types within the query data—a feature largely unaddressed by previous methods.
* **High Performance on Rare Cell Types**: The pipeline demonstrates a superior ability to identify unknown cell types, including rare populations with as few as 20 cells.

---

## ⚙️ Workflow

The HiCat pipeline consists of six main steps:

1.  **Batch Effect Removal**: Batch effects between reference and query datasets are corrected using Harmony on the top 50 principal components.
2.  **Dimensionality Reduction**: The 50-dimensional harmonized embedding is further reduced to 2 dimensions using UMAP to capture key local and global patterns.
3.  **Unsupervised Clustering**: DBSCAN is applied to the 2D UMAP embedding to identify dense clusters of cells, which can represent novel cell type candidates.
4.  **Feature Space Concatenation**: The outputs from the previous steps (50 PCs + 2 UMAP dimensions + 1 DBSCAN cluster label) are concatenated to create a 53-dimensional multi-resolution feature space.
5.  **Supervised Classification**: A CatBoost classifier is trained on the 53-dimensional features from the reference set. This trained model is then used to predict cell types for the query set.
6.  **Resolving Discrepancies**: A novel strategy resolves conflicts between the supervised (CatBoost) prediction and the unsupervised (DBSCAN) cluster label. The CatBoost prediction is used when confidence is high; otherwise, the DBSCAN label is assigned, enabling the robust identification of unseen cell types.

---

## 📊 Performance Highlights

In comprehensive evaluations using 10 public genomic datasets, HiCat demonstrated superior performance:

* **Seen Cell Types**: Achieved the highest overall accuracy compared to nine other well-established methods when all query cell types were present in the reference set.
* **Unseen Cell Types**: Remained the most resilient and accurate method when query sets included one, two, or three unseen cell types.
* **Distinguishing Novel Types**: Successfully identified and distinguished four major unseen cell types in a human pancreas dataset, assigning them distinct labels with over 90% accuracy.

---

## 💻 Installation

To set up the environment for running the demo, we recommend using `conda`.

1.  **Clone the repository:**
    ```bash
    git clone [https://github.com/changbiHub/HiCat.git](https://github.com/changbiHub/HiCat.git)
    cd HiCat
    ```

2.  **Create and activate a conda environment:**
    ```bash
    conda create --name hicat-env python=3.9
    conda activate hicat-env
    ```

3.  **Install Python dependencies:**
    The required Python packages are listed in `requirements.txt`.
    ```bash
    pip install -r requirements.txt
    ```

4.  **Install R dependencies:**
    Launch an R session and run the following command:
    ```R
    install.packages(c("Seurat", "harmony", "dplyr", "ggplot2", "cowplot"))
    ```

---

## 🚀 Usage

This demo showcases the core functionality of the HiCat pipeline.

1.  **Download the data:**
    Download the dataset pair (`train.rds` and `test.rds`) from the link in the [Data](#-data) section and place them in the `data/` directory.

2.  **Run the pipeline:**
    ```bash
    python script/run.py
    ```

3.  **Inspect the results:**
    Output files will be saved in the `outputs/` folder. The `walkThough.ipynb` notebook can be used to explore the results.

---

## 💾 Data

The data for Experiment 3 of the paper can be downloaded from Google Drive:

[**Download Dataset**](https://drive.google.com/drive/folders/1gjBLkGrORXKwmiRdb860viwZK3mOIrcn?usp=sharing)

---

## 📜 Citation

If you use HiCat in your research, please cite our preprint:

```bibtex
@misc{bi2024hicat,
      title={HiCat: A Semi-Supervised Approach for Cell Type Annotation}, 
      author={Chang Bi and Kailun Bai and Xuekui Zhang},
      year={2024},
      eprint={2412.06805},
      archivePrefix={arXiv},
      primaryClass={q-bio.GN}
}
