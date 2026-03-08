# Differential Gene Expression Analysis of Blood Samples from Patients with Chronic Spontaneous Urticaria

## Project Overview
This project investigates the molecular mechanisms underlying **Chronic Spontaneous Urticaria (CSU)** through **differential gene expression analysis** of blood samples. The analysis identifies genes that are significantly upregulated or downregulated in CSU patients compared to healthy controls and explores their biological significance through functional enrichment analysis.

The study utilizes publicly available microarray data from the **Gene Expression Omnibus (GEO)** database to analyze gene expression patterns associated with CSU.

---

## Dataset

- **Accession ID:** GSE72541  
- **Source:** Gene Expression Omnibus (GEO), NCBI  
- **Platform:** Agilent SurePrint G3 Human Gene Expression v2 8×60K Microarray (GPL16699)  

**Samples**
- 20 patients with Chronic Spontaneous Urticaria (CSU)  
- 10 healthy control individuals  

**Original study:**  
Giménez-Arnau et al. (2017)

---

## Objectives

The main objectives of this analysis are:

- Identify **Differentially Expressed Genes (DEGs)** between CSU patients and healthy controls
- Explore **biological processes and pathways** associated with CSU
- Investigate the **molecular mechanisms** contributing to CSU pathophysiology

---

## Methods

### 1. Data Retrieval
The dataset was retrieved from the **GEO database** using the `GEOquery` package in R.

### 2. Data Preprocessing
- Log2 transformation of gene expression values
- Normalization and quality assessment
- Visualization using:
  - Boxplots
  - Density plots
  - UMAP

### 3. Differential Expression Analysis
Differential gene expression analysis was performed using the **limma** package.

Filtering criteria:
- Adjusted p-value (**adj.p.val**) < 0.05
- |logFC| > 0.1

### 4. Data Visualization
Results were visualized using:
- **Volcano plots**
- **Heatmaps**
- **UMAP clustering**

### 5. Functional Enrichment Analysis
Significant genes were analyzed using **g:Profiler** to perform:

- **Gene Ontology (GO) enrichment analysis**
- **KEGG pathway enrichment analysis**

---

## Key Findings

### Platelet Activation and Hemostasis
GO and KEGG analyses highlighted significant enrichment in pathways related to **platelet activation and blood coagulation**, suggesting that platelets may contribute to inflammatory responses in CSU.

### Immune System Activation
Genes associated with **adaptive immune regulation**, particularly **T-cell activation**, were also identified, supporting the hypothesis that CSU involves **autoimmune mechanisms**.

### PI3K–AKT Signaling Pathway
The **PI3K–AKT signaling pathway** was significantly enriched, indicating its role in:

- Platelet hyper-responsiveness
- Sustained inflammatory responses

These findings suggest that **interactions between the coagulation system and the immune system** play a critical role in the pathophysiology of CSU.

---

## Tools and Packages

### Software
- R
- RStudio

### R Packages
- GEOquery
- limma
- AnnotationDbi
- ggplot2
- pheatmap
- dplyr
- umap
