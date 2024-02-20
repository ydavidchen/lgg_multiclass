# Multiclass prediction of lower-grade glioma subtypes via gene promoter CpG island profiles: a tutorial

David Chen, Ph.D.

## Repository Structure

Each cohort (TCGA - American or GCN - German) has the following scripts tagged:

- dataprep.R: data curation 
- hclust.R: unsupervised hierarchical clustering
- pca.R: principal component analysis
- holdout.py: 80:20 hold-out experiment
- cv.py: 5-fold cross-validation experiment

Additional scripts include utils.R (R helpers) and ml_utils.py (Python ML helpers) shared across the project.

## Notes on Software Used

R 4.3.1 was used for preparing and curating the data from raw data sources as well as exploratory data analyses.

Python 3.12 was used for machine learning.  

## Other Notes

The curated datasets containing the computed CGI methylation values and clinical subtypes (CSV output of scripts tagged "_dataprep.R") were deposited in [Kaggle](https://www.kaggle.com/datasets/ydavidchen/lower-grade-glioma-cpg-islands-and-subtypes/).

The unsupervised and supervised ML analyses are excellent choices for teaching purposes and demos.


