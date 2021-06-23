# Cancer Health Disparities through a genomic epidemiology lens
Clinical genomics-based tests used in cancer patients have the potential to reinforce health disparities, especially if the targets of tests do not offer prognostic validity for people of diverse ancestries. For example, a widely-used 21-gene test established by [Paik et al. 2004](https://doi.org/10.1634/theoncologist.12-6-631) uses tumor genomic profiling in breast cancer patients to assess recurrence with recurrence scores and make recommendations on the benefits of chemotherapy for invasive breast cancer cases that are ER+ and HER2- at most anatomical stages and regardless of lymph node invasion. [Hoskins et al. 2021](https://doi.org/10.1001/jamaoncol.2020.7320) demonstrated that not only do black women have disproportionately higher risk recurrence scores but also experience worse outcomes when matched against non-Hispanic white women with similar diagnoses and recurrence scores, meaning that the mist popular tumor profiling tests has poor prognostic value in black women. As such, we investigated whether any such divergence in the value of prognostic gene expression panels to colon cancer patients of African ancestry. Currently, African Americans experience disproportionately high rates of colon cancer incidence and mortality in comparison to NHP European Americans [Carethers. 2021](https://doi.org/10.1016/bs.acr.2021.02.007). The increased cancer burden in this specific patient population creates an urgent need for examination of the accuracy of commonly used prognostic tests. 

### Pipeline
(in progress)

1. Dataset: TCGA-COAD
2. Extract: Gene expression, CNV, ancestry
3. Stratify/infer race
4. Gene expression analysis
  * DEG in tumor vs healthy tissue in all patients
  * DEG in tumor vs healthy tissue by ethnicity / race
5. Gene set enrichment analysis for pathways of DEG
6. Primary outcome
  * DEG versus 5-year survival in all patients
  * DEG versus 5-year survival by ethnicity / race
7. Stratify analysis by all genes, genes overlapping clinical tests
  * measure how much information each gene adds to outcomes in a by-ancestry manner
8. ? Compartive immune infiltration by race using CIBERSORT
  * immune cells associated with PFS/OS, by race
  * immune clusters associated with PFS/OS, by race
