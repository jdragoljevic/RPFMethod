# RPFMethod
Method for Sequence Similarity Analysis Based on the Position and Frequency of Statistically Significant Repeats

R-P/F method is the alignment-free method for biological sequence similarity analysis based on the position and frequency of statistically significant repeats.
R-P/F method is implemented in R programming language. It supports fetching data from the database for the set of sequences, calculating local based frequency entropy, sequence similarity calculation, clustering as well as validation.
This repository contains R scripts as well as test data to be used in R-P/F analysis.

## Requirements

The following software tools are required:
1. [StatRepeats](http://bioinfo.matf.bg.ac.rs/home/downloads.wafl?cat=Software&project=StatRepeats) for extracting four types of statistically significant repeat sequences

The following R packages are required:
1. [RODBC] (Ripley B., Lapsley M., RODBC: ODBC Database Access, 2017.) that implements Open Database Connectivity (ODBC) for retrieving and managing data from IBM Db2 database;
2. [parallel] (R Core Team, R: A Language and Environment for Statistical Computing, 2016.) for running a bigger number of computations in parallel; 
3. [Lsa] (Wild F., lsa: Latent Semantic Analysis, 2015.) for calculating cosine similarity matrix once local frequency entropy and sequence vectors were calculated as well as for printing dendrograms; 
4. [stats] (R Core Team, R: A Language and Environment for Statistical Computing, 2016.) for clustering using agglomeration method average (UPGMA);  
5. [d3heatmap] (Cheng J.,Galili T., d3heatmap: Interactive Heat Maps Using 'htmlwidgets' and 'D3.js', 2016.) package is used for visualizing results as heatmaps; 
6. [purr] (Henry L., Wickham H., purrr: Functional Programming Tools, 2017.)   for visualizing results of cluster analysis;
7. [factoextra] (Kassambara A., Mundt F., factoextra: Extract and Visualize the Results of Multivariate Data Analyses, 2017.) for visualizing results of cluster analysis;
8. [corrplot] (Wei T., Simko V. , R package "corrplot": Visualization of a Correlation Matrix, 2017.) for comparing resulted dendrograms;
9. [dendextend] (Galili T., Benjamini Y., et al., R package "dendextend": Extending 'dendrogram' Functionality in R, 2020.) for comparing resulted dendrograms;

## Options


```
	repeat_type
	--dn	: Search for direct non-complementary repeats
	--dc    : Search for direct complementary repeats
	--in    : Search for inverse non-complementary repeats
	--ic   	: Search for inverse complementary repeats
     
```
## Scripts

### Pull data from IBMDB2 database
To get data about repeats from the IBMDB2 database please use script: RPF_get_data_from_ibmdb2_db.r.

### Load mitochondrial genome DNA sequences and RNA viruses: Ebolavirus, Marburgvirus, and Betacoronavirus genus strains sequences used in the research
To load data used in the research please use script: RPF_load_data.r as well as data files attached:
| file name | description |
--- | --- 
| [RData_mitoh_3seq](https://github.com/jdragoljevic/RPFMethod/blob/main/datasets/RData_mitoh_3seq) | test data for 3 mitochondrial DNA sequence: V00662.1 - Human, D38113.1 - Chimpanzee and AJ001588.1 - Rabbit, dn repeat type |
| [RData_mitoh_dn](https://github.com/jdragoljevic/RPFMethod/blob/main/datasets/RData_mitoh_dn) | test data for complete mitochondrial DNA sequences from 45 different mammalian species with dn repeat type |
| [RData_mitoh_dc](https://github.com/jdragoljevic/RPFMethod/blob/main/datasets/RData_mitoh_dc) | test data for complete mitochondrial DNA sequences from 45 different mammalian species with dc repeat type |
| [RData_mitoh_in](https://github.com/jdragoljevic/RPFMethod/blob/main/datasets/RData_mitoh_in) | test data for complete mitochondrial DNA sequences from 45 different mammalian species with in repeat type |
| [RData_mitoh_ic](https://github.com/jdragoljevic/RPFMethod/blob/main/datasets/RData_mitoh_ic) | test data for complete mitochondrial DNA sequences from 45 different mammalian species with ic repeat type |
| [RData_viruses_dn](https://github.com/jdragoljevic/RPFMethod/blob/main/datasets/RData_viruses_dn) | test data for 383 RNA viruses dn repeat type |
| [RData_viruses_dc](https://github.com/jdragoljevic/RPFMethod/blob/main/datasets/RData_viruses_dc) | test data for 383 RNA viruses dc repeat type |
| [RRData_xylanase_in](https://github.com/jdragoljevic/RPFMethod/blob/main/datasets/RData_xylanase_in) | test data for 20 protein sequences |
| [RData_xylanase_dn](https://github.com/jdragoljevic/RPFMethod/blob/main/datasets/RData_xylanase_dn) | test data for 20 protein sequences |
| [RData_transfferin_in](https://github.com/jdragoljevic/RPFMethod/blob/main/datasets/RData_transfferin_in) | test data for 24 protein sequences |
| [RData_transfferin_dn](https://github.com/jdragoljevic/RPFMethod/blob/main/datasets/RData_transfferin_dn) | test data for 24 protein sequences |
| [RData_spike_coronavirus_in](https://github.com/jdragoljevic/RPFMethod/blob/main/datasets/RData_spike_coronavirus_in) | test data for 50 protein sequences |
| [RData_spike_coronavirus_dn](https://github.com/jdragoljevic/RPFMethod/blob/main/datasets/RData_spike_coronavirus_dn) | test data for 50 protein sequences |
| [RData_ND6_in](https://github.com/jdragoljevic/RPFMethod/blob/main/datasets/RData_ND6_in) | test data for 8 protein sequences |
| [RData_ND6_dn](https://github.com/jdragoljevic/RPFMethod/blob/main/datasets/RData_ND6_dn) | test data for 8 protein sequences |
| [RData_ND5_in](https://github.com/jdragoljevic/RPFMethod/blob/main/datasets/RData_ND5_in) | test data for 9 protein sequences |
| [RData_ND5_dn](https://github.com/jdragoljevic/RPFMethod/blob/main/datasets/RData_ND5_dn) | test data for 9 protein sequences |
| [RData_betaglobin_in](https://github.com/jdragoljevic/RPFMethod/blob/main/datasets/RData_betaglobin_in) | test data for 50 protein sequences |
| [RData_betaglobin_dn](https://github.com/jdragoljevic/RPFMethod/blob/main/datasets/RData_betaglobin_dn) | test data for 50 protein sequences |

### calculating local based frequency entropy, sequence similarity calculation, clustering as well as validation
please use script: RPF_similarity.r

## System Requirement

R environment

## Cite

DOI: 10.2174/1574893616999210805165628

## License

R-P/F method, it's implementation programs and data attached should be used only for purely academic purposes and commercial use is strictly prohibited.
Work should be cited (DOI: 10.2174/1574893616999210805165628).
Please contact [Jasmina Jovanovic](mailto:jasmina.dragoljevic@gmail.com) for any commercial license.


