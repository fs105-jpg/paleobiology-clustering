# Paleobiology Database Cluster Analysis

## Overview
This project analyzes fossil occurrence data from multiple geographic regions to identify patterns in fossil localities using unsupervised learning.

The workflow includes:
- merging and cleaning fossil datasets  
- engineering geological features (lithology, material, environment)  
- clustering using Gower distance and PAM  
- visualizing spatial and compositional patterns  

---

## Data
The data is derived from fossil occurrence records across multiple regions.

Due to size and source constraints:
- raw datasets are not included  
- a cleaned dataset is used for clustering  

Required file:
data/clusterReady_clean.csv

This dataset includes:
- latitude, longitude  
- lithology grouping (lith_combo)  
- fossil material grouping (material_grp)  
- environment grouping (env_group)  
- optional time variables (max_ma, min_ma)  

---

## Methods

### Data Cleaning
- standardized fossil composition into categories:
  - carbonate, phosphate, silica, organic, mixed, none  
- grouped lithology into:
  - carbonate, clastic, volcanic, silica/chemical, organic  
- categorized depositional environments:
  - marine, marginal marine, fluvial, lacustrine, terrestrial  

### Clustering
- distance metric: Gower distance (handles mixed data types)  
- algorithm: PAM (Partitioning Around Medoids)  
- tested k = 3–8  
- selected optimal k using silhouette scores  

### Visualization
- spatial clustering (latitude vs longitude)  
- cluster composition by:
  - geologic period  
  - environment  
  - lithology  
  - fossil material  

---

## Results
The analysis identifies clusters of fossil localities with distinct:
- geographic distributions  
- depositional environments  
- lithology and material composition  

---

## How to Run

1. Clone the repository:
git clone https://github.com/fs105-jpg/paleobiology-clustering.git

2. Add dataset:
data/clusterReady_clean.csv

3. Run:
paleobiology_clustering.R

---

## Key Takeaways
- Fossil localities cluster meaningfully based on geological features  
- Gower distance is effective for mixed categorical data  

---

## Author
Felicia Selbst  
Wellesley College — Data Science Major, Chemistry Minor
