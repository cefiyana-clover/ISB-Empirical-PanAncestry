# ISB-Empirical-PanAncestry
Empirically-validated computational psychiatry model for MDD using Pan-Ancestry genomics and clinical NCBI data.
# Bioenergetic Stability Index (ISB): Empirical Validation & Pan-Ancestry Population Genomics

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.19317017.svg)](https://doi.org/10.5281/zenodo.19317017)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Python 3.8+](https://img.shields.io/badge/python-3.8+-blue.svg)](https://www.python.org/downloads/)

## 1. Overview
This repository contains the official empirically-driven computational pipeline for the **Bioenergetic Stability Index (ISB)**. This research bridges theoretical thermodynamic modeling with real-world clinical and genomic data to reconstruct **Major Depressive Disorder (MDD)** as a systemic mitochondrial shutdown within the astrocytic network.

The core mechanism identifies a pathological bridge where peripheral allostatic loads—specifically **Gastroesophageal Reflux Disease (GERD)** and **Hypocapnia ($PaCO_2$)**—are communicated via vagal afferent hyperactivation, triggering a deterministic **Saddle-Node Bifurcation** in brain energy equilibrium.

## 2. Key Empirical Findings
Unlike synthetic models, this pipeline utilizes parameters extracted from Q1-scale clinical literature (NCBI PubMed Central) and Pan-Ancestry GWAS data. Key highlights include:
* **Clinical Distinction:** The model accurately separates healthy cohorts from clinical cohorts ($PaCO_2$ < 35 mmHg) using non-linear ODEs.
* **Genetic Vulnerability:** Under "Borderline Stress" ($PaCO_2$ = 39.2 mmHg), the simulation reveals an extreme genetic vulnerability in **East Asian (EAS)** populations (80.1% collapse rate) compared to European (2.5%) and African (1.4%) cohorts.
* **Eurocentric Bias Correction:** These results mathematically expose the danger of applying Western diagnostic thresholds to diverse global populations.

## 3. Repository Structure
* `/scripts/01_ncbi_data_miner.py`: Automated API crawler for clinical and GWAS metadata extraction.
* `/scripts/02_clinical_validation.py`: ODE-based simulation validating the ISB threshold (> 2.50) using real-world patient parameters.
* `/scripts/03_panancestry_megasim.py`: Mega-simulation (N=30,000) mapping Gene-Environment (GxE) interactions across EUR, EAS, and AFR ancestries.
* `/results/figures/`: High-resolution (300 DPI) Phase-Space scatter plots, KDE density plots, and Pan-Ancestry Violin plots.

## 4. Installation & Reproducibility
```bash
git clone [https://github.com/cefiyana-clover/ISB-Empirical-PanAncestry.git](https://github.com/cefiyana-clover/ISB-Empirical-PanAncestry.git)
cd ISB-Empirical-PanAncestry
pip install -r requirements.txt
5. Mathematical Architecture
​The system evaluates the trajectory of astrocytic ATP (A) and extracellular glutamate (G):
$$ \frac{dA}{dt} = P_{ATP}(PaCO_2) - C_{basal} - 4 \cdot V(A, G) $$
$$ \frac{dG}{dt} = R_{stress} - V(A, G) - k_{leak} \cdot G $$
​6. Citation
​If you utilize this framework, the code, or the empirical findings, please cite the associated Zenodo publication:
​Cefiyana. (2026). Empirical Data and Pan-Ancestry Validation Methodology of the Bioenergetic Stability Index (ISB). Zenodo. DOI: 10.5281/zenodo.19317017
​7. License
​This project is licensed under the MIT License.
