**NOTE: This repository provided all source data and codes for generation of all results in the single-tooth ECC manuscript.**

## Title: Single-tooth resolved, whole-mouth prediction of early childhood caries via spatiotemporal variations of plaque microbiota

## Brief summary

Early Childhood Caries (ECC) remains a prevalent and serious health concern worldwide, impacting both individual health and public health systems. To further enhance the precision of ECC monitor and intervene at the single-tooth level, we longitudinally tracked single-tooth dental plaque microbiomes from children with and without ECC, providing the microbiome insights into the tooth-specific etiology of ECC, and developed a comprehensive machine-learning model predicting the ECC at the tooth. Specifically, this longitudinal study of 89 preschoolers over a 13-month period allowed us to characterize distinct microbiome trajectories specific to each tooth, which could then be used to predict ECC onset with unprecedented accuracy.

Overall, our findings revealed the spatiotemporal microbiota dynamics that differentiate healthy and caries-active children, including the anterior-posterior ecological gradients for maxillary teeth and spatial symmetry in healthy children, which were disrupted in the presence of caries. We further proposed the Spatial Microbial Indicators of Caries (Spatial-MiC), a diagnostic tool that incorporating both on-site microbiome and the whole-mouth microbial/clinical spatial features, with 98% accuracy, could identify ECC at single-tooth resolution and predict future ECC onset for clinically healthy teeth up to two months in advance with 93% accuracy. Our work represents the only study to date at the *single-tooth*, *full-mouth*, and *cohort scale*, providing a foundation for a *tooth-specific*, *high-resolution* approach to ECC prevention and a microbial map for understanding the onset and progression of caries at the single-tooth level. 

## Data availability

The metagenomic sequencing data used for the analysis presented in this study are available from Qiita (https://qiita.ucsd.edu/study/description/14341). 

## Key resources table
| REAGENT or RESOURCE | SOURCE | IDENTIFIER |
|-------|-------------------|--------------------------|
| **Biological samples** |
| Plaque samples | This study | N/A |
| **Software and algorithms** |
| FastQC v0.11.9 | |https://www.bioinformatics.babraham.ac.uk/projects/fastqc|
| KneadData v0.7.2 | | https://huttenhower.sph.harvard.edu/kneaddata/ |
| HUMAnN2 | Franzosa et al. | https://huttenhower.sph.harvard.edu/humann2/ |
| PICRUSt2 | Douglas et al. | https://huttenhower.sph.harvard.edu/picrust/ |
| Phylo-rPCA | Martino et al. | https://github.com/biocore/gemelli/blob/master/ipynb/tutorials/Phylogenetic-RPCA-moving-pictures.ipynb |
| rPCA | Martino et al. | https://github.com/biocore/gemelli/blob/master/ipynb/tutorials/RPCA-moving-pictures.ipynb |
| Python (v.3.9.6) | Python | https://www.python.org/downloads/release/python-396 |
| R (v.4.1.3) | R Foundation | https://www.r-project.org |

## ACKNOWLEDGEMENTS
This work was funded by Grants 32030003 and 31600099 from National Natural Science Foundation of China, and SKLOD2015OF07 from State Key Laboratory of Oral Diseases. SH thanks the funding support from the Health and Medical Research Fund (HMRF) of Hong Kong (10212276). FY thanks the funding support from National Natural Science Foundation of China (31300424, 81670979) and Qingdao Natural Science Foundation, China (24-4-4-zrjj-158-jch). FT thanks the funding support from General Program of Natural Science Foundation of Shandong Province (ZR2024MH235) and Taishan Scholar Award For Young Expert (tsqn201909126), as well as the funding supports from Qingdao Key Health Discipline Development Fund (2022-2024), Qingdao Clinical Research Center for Oral Diseases (22-3-7-lczx-7-nsh), and Shandong Provincial Key Medical and Health Discipline of Oral Medicine. 
