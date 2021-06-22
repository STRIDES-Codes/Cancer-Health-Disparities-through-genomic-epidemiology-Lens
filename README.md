# Cancer Health Disparities through a genomic epidemiology lens
Clinical genomics-based tests used in cancer patients have the potential to reinforce health disparities, especially if the targets of tests do not offer prognostic validity for people of diverse ancestries. For example, a widely-used 21-gene test established by [Paik et al. 2004](https://doi.org/10.1634/theoncologist.12-6-631) uses tumor genomic profiling in breast cancer patients to assess recurrence with recurrence scores and make recommendations on the benefits of chemotherapy for invasive breast cancer cases that are ER+ and HER2- at most anatomical stages and regardless of lymph node invasion. [Hoskins et al. 2021](https://doi.org/10.1001/jamaoncol.2020.7320) demonstrated that not only do black women have disproportionately higher risk recurrence scores but also experience worse outcomes when matched against non-Hispanic white women with similar diagnoses and recurrence scores, meaning that the mist popular tumor profiling tests has poor prognostic value in black women. Could there exist pervasive disparities in clinical tumor oncotyping tests that we might address through re-examining gene expression and copy number variation stratified by ancestry?

### Pipeline
(in progress)

1. select cancer (TCGA-COAD is colon)
2. pull TCGA data (gene expression and CNV) (TCGA) and metadata (ancestry) (TCGAA)
3. stratify by ancestry and produce summary statistics
4. ??? differential analyses ???
  * DEG in tumor vs healthy tissue in all patients
  * DEG in tumor vs healthy tissue by ethnicity / race
5. Gene set enrichment analysis for pathways of DEG
6. ??? outcome analyses ???
  * DEG versus 5-year survival in all patients
  * DEG versus 5-year survival by ethnicity / race
7. stratify analysis by all genes, genes overlapping clinical tests
  * measure how much information each gene adds to outcomes in a by-ancestry manner

Team
Sheryse
