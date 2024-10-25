# HiCat: A Semi-Supervised Approach for Cell Type Annotation

![Graphical Abstract](graphicAbstract.png)

This is the repo for the demo code of HiCat (Hybrid Cell Annotation using Transformative embeddings), a novel semi-supervised pipeline for annotating cell types from single-cell RNA sequencing data. HiCat fuses the strengths of supervised learning for known cell types with unsupervised learning to identify novel types. This hybrid approach incorporates both reference and query genomic data for feature engineering, enhancing the embedding learning process, increasing the effective sample size for unsupervised techniques, and improving the transferability of the supervised model trained on reference data when applied to query datasets.

---
A complete package is working in progress. To run this demo, first download the dataset pair(train.rds and test.rds) and put them under the data folder, then run 

```
python script/run.py
```

the results would be generated in the outputs folder. Once the result is generated, walkThough.ipynb can be used to inspect the results.

Data for experiment 3 of the paper can be download from:
[Google Drive](https://drive.google.com/drive/folders/1gjBLkGrORXKwmiRdb860viwZK3mOIrcn?usp=sharing)
