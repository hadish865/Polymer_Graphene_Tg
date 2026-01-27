# Machine Learning Prediction of Glass Transition Temperature in Polymer/Graphene Nanocomposites

this repository contains the data and source codes associated with the paper:
"Interaction-Based Molecular Descriptors for machine learning Prediction of Glass Transition Temperature in
Polymer/Graphene Nanocomposites"

Authors: Bardia Afsordeh, Hadi Shirali

Journal: Computational Materials Science

## Contents
data/: processed data used for model training
descriptors/: Codes for generating molecular and interaction-based descriptors
models/: Training and testing scripts for machine learning models
results/: Prediction results

## Prerequisite
in order to run the train_models.py script (located at models/train_models.py) the following packages are required:
numpy (pip install numpy)
scikit-learn (pip install -U scikit-learn)
pandas (pip install pandas)
xlsxwriter (pip install xlsxwriter)
pickle (this package usually comes pre installed on recent pyhon versions)

in order to run the Descriptors_test.py script (located at descriptors\Descriptors_test.py) the following packages are required:
xlsxwriter (pip install xlsxwriter)
rdkit (pip install rdkit)
PIL (pip install pillow)

## how to run
to run the descriptors test on windows from the main directory:
python.exe descriptors\Descriptors_test.py

to run the machine learning models test on windows from the main directory:
python.exe models\train_models.py
