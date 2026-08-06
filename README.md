# Rhythmate

## Steps to successfully run the **rhythmate** analysis

---

***Step 1***


Load the required data:   

``
- rhythmate_cancer_data <- readRDS("rhythmate_cancer_data.rds")
- rhythmate_normal_data <- readRDS("rhythmate_normal_data.rds")
- rhythmate_cancer_clinical <- readRDS("rhythmate_cancer_clinical.rds")
- rhythmate_normal_clinical <- readRDS("clinical_normal_tcga.rds")
- rhythmate_survival_data <- read.xlsx("rhythmate_survival_data.xlsx")
``

---

***Step2***

run ``run_rhythmate_normal_analysis`` in an object named *out_norm* and define your desired gene names

``out_norm <- run_rhythmate_normal_analysis(rhythmate_normal_data, rhythmate_normal_clinical, group1, group2, group3, iterations = 20)``

---

***Step3***

run ``run_rhythmate_analysis`` in an object named *out_can*

``out_can <- run_rhythmate_analysis(rhythmate_cancer_data, rhythmate_cancer_clinical,
                                  out_norm$ref.mat.mean, out_norm$ref.ccd.cir.mean,
                                  out_norm$ref.ccd.cell.mean, group1, group2, group3)``

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


There are 3 output files with rhythmates' prefix:
> 1. A multipage pdf file containing:
>    - clustered heatmap
>    - expression box plot
>    - CCD plots
>    - survival plots
> 2. An excel file containing COX models 
> 3. An csv file inlcuding the cluster assignments