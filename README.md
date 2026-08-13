# Rhythmate

## Steps to successfully run the **rhythmate** analysis

---

***Step 1***


Load the required data:   

``
- rhythmate_cancer_data <- readRDS("rhythmate_cancer_data.rds")
- rhythmate_normal_data <- readRDS("rhythmate_normal_data.rds")
- rhythmate_cancer_clinical <- readRDS("rhythmate_cancer_clinical.rds")
- rhythmate_normal_clinical <- readRDS("rhythmate_normal_clinical.rds")
- rhythmate_survival_data <- read.xlsx("rhythmate_survival_data.xlsx")
- rhythmate_cancer_qn <- readRDS("rhythmate_cancer_qn.rds")
- rhythmate_normal_qn <- readRDS("rhythmate_normal_qn.rds")
``

---

***Step2***

Defining your genes:
````
group1 <- toupper(c("Clock","Arntl","Npas2","Nfil3","Cry1","Cry2", "Per1","Per2","Per3", Nr1d1","Nr1d2" ,"Bhlhe40","Dbp" ,"Hlf" ,"Tef"))
group1 <- append(group1, "C1orf51", after = which(group1 == "CRY2"))
group2 <- toupper(c("Cdk1", "Ccna2", "Plk1", "Mcm2", "E2f1", "Ccnd1", "Cdk4", "Pcna", "Rpa1", "Rb1", "Cdc6", "Cdc20", "Aurka", "Cdc25c", "Cenpe", "Ccnb1"))
group3 <- c("MYC", "MYCN", "MKI67")
````


run ``run_rhythmate_normal_analysis`` in two objects named ***out_norm*** and ***out_norm_cluster***. </br> This function should be run twice (*depending on whether you have applied log2 to your expression data*): 
1. Once with apply_log2_transform = FALSE
 
    - *here we input quantile normalized data where log2 already applied*
    ````
    out_norm <- run_rhythmate_normal_analysis(rhythmate_normal_qn, rhythmate_normal_clinical, group1, group2, group3, iterations = 20, apply_log2_transform = FALSE)
    ````

2. Once with apply_log2_transform = TRUE 

    - *here we input NON-normalized data where log2 has not been applied*
    ````
    out_norm_cluster <- run_rhythmate_normal_analysis(rhythmate_normal_data, rhythmate_normal_clinical, group1, group2, group3, iterations = 20, apply_log2_transform = TRUE)
    ````

---

***Step3***

run ``run_rhythmate_analysis`` in two objects named ***out_can*** and ***out_can_cluster***.
</br> This function should be run twice (*depending on whether you have applied log2 to your expression data*): 
1. Once with apply_log2_transform = FALSE
 
    - *here we input quantile normalized data where log2 already applied*
    ````
    out_can <- run_rhythmate_analysis( rhythmate_cancer_qn, rhythmate_cancer_clinical,  out_norm$ref.mat.mean, out_norm$ref.ccd.cir.mean, out_norm$ref.ccd.cell.mean, group1, group2, group3,  apply_log2_transform = FALSE,  normal_reference_cross = out_norm$normal_dataframe,   normal_reference_cir = out_norm$normal_dataframe.cir,  normal_reference_cell = out_norm$normal_dataframe.cell)
    ````

2. Once with apply_log2_transform = TRUE 

    - *here we input NON-normalized data where log2 has not been applied*
    ````
    out_can_cluster <- run_rhythmate_analysis(rhythmate_cancer_data, rhythmate_cancer_clinical,   out_norm_cluster$ref.mat.mean, out_norm_cluster$ref.ccd.cir.mean, out_norm_cluster$ref.ccd.cell.mean,  group1, group2, group3,  apply_log2_transform = TRUE,   normal_reference_cross = out_norm_cluster$normal_dataframe,  normal_reference_cir = out_norm_cluster$normal_dataframe.cir,  normal_reference_cell = out_norm_cluster$normal_dataframe.cell)
    ````

---

***Step4***



run the following helper functions: 

````
rhythmate_survival_analysis 
rhythmate_export_cox_excel
````

---

***Step 5***

run

 ``run_rhythmate_clustering``

 ````
 run_rhythmate_clustering(out_norm, out_can,  rhythmate_survival_data = rhythmate_survival_data.xlsx",  out_norm_cluster = out_norm_cluster,  out_can_cluster = out_can_cluster)
  ````


There are 3 output files with Rhythmates' prefix:
> 1. A multipage pdf file containing:
>    - clustered heatmap
>    - expression box plot
>    - CCD plots
>    - survival plots
> 2. An excel file containing COX models 
> 3. An csv file inlcuding the cluster assignments