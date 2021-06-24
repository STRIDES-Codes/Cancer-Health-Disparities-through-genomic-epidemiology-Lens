# Cancer Health Disparities through a genomic epidemiology lens
### Team: David Adeleke (Leader), David Enoma (Mentor), Rony Arauz, Kathleen Morrill, Sheryse Taylor, Thai Tran, Erika Whitney

Clinical genomics-based tests have been widely employed to aid physicians in determining a patient's risk and to inform therapeutic decisions. However, they can potentially reinforce health disparities if the tests are not prognostically valid for patients of all ancestries. For example, a widely-used 21-gene test established by [Paik et al. 2004](https://doi.org/10.1634/theoncologist.12-6-631) uses tumor genomic profiling in breast cancer patients to assess recurrence with recurrence scores and make recommendations on the benefits of chemotherapy for invasive breast cancer cases that are ER+ and HER2- at in early stages and regardless of lymph node invasion. [Hoskins et al. 2021](https://doi.org/10.1001/jamaoncol.2020.7320) demonstrated that not only do black women have disproportionately higher risk recurrence scores but also experience worse outcomes when matched against non-Hispanic white women with similar diagnoses and recurrence scores, meaning that the most popular tumor profiling test has poor prognostic value in black women. Thus, we investigated whether there were any differences in the prognostic validity of gene panel tests for colon cancer patients of African descent. Currently, African Americans experience disproportionately high rates of colon cancer incidence, and mortality in comparison to non-Hispanic European Americans [Carethers. 2021](https://doi.org/10.1016/bs.acr.2021.02.007). The increased colon cancer burden in this specific patient population creates an urgent need for examination of the accuracy of commonly used prognostic tests. 

### Pipeline
(in progress)

1. Dataset: TCGA-COAD
2. Extract: Gene expression, CNV, ancestry
3. Stratify/infer race
4. Gene expression analysis
  * DEG in tumor vs healthy tissue in all patients
  * DEG in tumor vs healthy tissue by ethnicity / race
5. Gene set enrichment analysis for pathways of DEG
  * GSEA by race of specific pathways found in gene expression panels
7. Primary outcome
  * DEG versus 5-year survival in all patients
  * DEG versus 5-year survival by ethnicity / race
8. Stratify analysis by all genes, genes overlapping clinical tests
  * measure how much information each gene adds to outcomes in a by-ancestry manner
9. Compartive immune infiltration by race using CIBERSORT
  * immune cells associated with PFS/OS, by race



### Workflow Graphic

![Workflow Graphic](https://github.com/STRIDES-Codes/Cancer-Health-Disparities-through-genomic-epidemiology-Lens/blob/main/Screen%20Shot%202021-06-24%20at%2012.13.58%20PM.png)




# Table KM-survival genes after SA - White Population
![Workflow Graphic](https://github.com/STRIDES-Codes/Cancer-Health-Disparities-through-genomic-epidemiology-Lens/blob/main/Survival_Correlations_White.PNG)



## Table KM-survival genes after SA - Black Population
![Workflow Graphic](https://github.com/STRIDES-Codes/Cancer-Health-Disparities-through-genomic-epidemiology-Lens/blob/main/Survival_Correlations_black.PNG)




## Table KM-survival by races
![Workflow Graphic](https://github.com/STRIDES-Codes/Cancer-Health-Disparities-through-genomic-epidemiology-Lens/blob/main/survival_races.PNG)
## Findings

### List of packages used in analysis

### Link to code

