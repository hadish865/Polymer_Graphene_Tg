# Machine Learning Prediction of Glass Transition Temperature in Polymer/Graphene Nanocomposites
This repository contains the datasets and source codes associated with the article:
"Interaction-Based Molecular Descriptors for Machine Learning Prediction of Glass Transition Temperature in Polymer/Graphene Nanocomposites"

Authors: Bardia Afsordeh, Hadi Shirali  
Journal: Computational Materials Science

## Overview
The purpose of this repository is to provide full access to the data and computational tools used in the development and validation of machine learning models for predicting the glass transition temperature (Tg) of polymer/graphene nanocomposites. The codes implement interaction-based molecular descriptors and supervised machine learning workflows as reported in the associated publication.
This repository is intended to enable full reproducibility of the results presented in the paper.

## Repository Structure
- `data/`  
  Processed datasets used for model training and testing.

- `descriptors/`  
  Python codes for generating molecular and interaction-based descriptors.

- `results/`  
  Output files containing prediction results and performance metrics.

## Requirements

The following Python packages are required to reproduce the results:

- numpy  
- pandas  
- scikit-learn  
- xlsxwriter  
- rdkit  
- pillow  

The `pickle` module is part of the standard Python library and does not require separate installation.

## How to Run
From the root directory of the repository:

### Descriptor generation
python descriptors/Descriptors_test.py

### Machine learning model training and testing
```bash
python train_models.py

## Note
-The datasets provided in this repository are the processed versions used directly in the machine learning workflows described in the manuscript.

-All scripts were tested using standard Python environments. Minor adjustments may be required depending on the operating system and Python version.

## Reproducibility
All numerical results, figures, and performance metrics reported in the associated article can be reproduced using the data and scripts provided in this repository.

```bash
python descriptors/Descriptors_test.py
