# GLMGE
## The prediction of gene expression by generalized linear regression model (GLM) and the prediction of co-regulated genes by graph embedding (GE)
<img src="https://github.com/Park-Sung-Joon/GLMGE/assets/52985953/afedbb99-cfcc-4564-b751-1b0d1b215bf0" width=80%>

### Requirements
perl and R libraries
+ library(MASS)
+ library(limma)
+ library(coefplot)
+ library(cowplot)
+ library(LSD)
+ library(ggplot2)

### What GLM does?
The input (TXT file) is separated into dependent variables (Y) and explanatory variables (X) and the regression model **Y = wX + e** is built.
1. Starting with a Full model (FM); X is all the explanatory variables in the input file.
2. AIC (Akaike's Information Criterion) reduces FM, which yields a reduced model **M1** and the set of removed features **XX** from the X, and N-fold cross-validation (CV) measures the predictive performance of **M1**.
3. **XX** is randomly shuffled. (searching the first-order effects)
4. Loop1: **M1** is extended with each of **XX**, which yields a trial model **M2**.
5. AIC and N-fold CV evaluate **M2**, which yields a new removed feature set **XX2**.
6. If **M2** improves **M1**, **M1** and **XX** are replaced with **M2** and **XX2**. Then, go to Step 3. :point_up_2:
7. Otherwise, go to Step 4. till all **XX** are tried. :point_up_2:
8. End of Loop1
9. **XX** is randomly shuffled. (searching the second-order effect)
10. Loop2: **M1** is extended with each of all the possible pairs in **XX**, which yields a trial model **M2**.
11. AIC and N-fold CV evaluate **M2**, which yields a new removed feature set **XX2**.
12. If **M2** improves **M1**, **M1** and **XX** are replaced with **M2** and **XX2**. Then, go to Step 3.:point_up_2:
13. Otherwise, go to Step 10. till all pairs of **XX** are tried. :point_up_2:
14. End of Loop2

### Command line
The following command line executes GLM.
```
%>perl GLMGE_v4.pl [INPUT] [OUTPUT] [ID] [Log10? 0|1] [RND_SEED] [SKIPgene_Filename|NULL] [SKIPcols: NULL|"A,B,C..."]
```
"GLMGE.R" must exist in the same directory where "GLMGE_v4.pl" is located.

### Input file format
The input file must be formatted as;
```
#[space]ITEM:[space]nameFeature1[space]nameFeature2[space]nameFeature4[space]nameFeature5...
GeneName[TAB]Expression[TAB]valueFeature1[space]valueFeature2[space]valueFeature3[space]valueFeature4...
GeneName[TAB]Expression[TAB]valueFeature1[space]valueFeature2[space]valueFeature3[space]valueFeature4...
GeneName[TAB]Expression[TAB]valueFeature1[space]valueFeature2[space]valueFeature3[space]valueFeature4...
```

### Example1: iterative single runs by removing outliers
```
%>perl GLMGE_v4.pl input_matrix first_result/result.txt ID 1 123456 NULL NULL
```
This creates the directory "first_result" and the output files.

<img src="https://github.com/Park-Sung-Joon/GLMGE/assets/52985953/61195c90-82a1-4c77-aa69-f47a8426bdba" width=50%>
+ (A) Pearson's correlation coefficient (R) between Observation (Obs.) and Prediction (Pred.) of the final regression model
+ (B) Regression coefficient (RC) of each feature in the full model and in the final model
+ (C) Red points representing outliers in the distribution of 90%, 95%, and 99% points.

You may want to remove a set of specific genes;
```
%>perl GLMGE_v4.pl input_matrix second_result/result.txt ID 1 123456 first_result/result.txt.outlier_99per.txt NULL
```
Here, the genes listed in "first_result/result.txt.outlier_99per.txt" are removed in this run. NOTE that the "first_result/result.txt.outlier_99per.txt" includes outlier genes found from the "first_result/" and the residuals were (+/- 2.58 in Zscore).

This is an example to run three times by removing outliers in each run.
```
%>perl GLMGE_v4.pl input_matrix first_result/result.txt ID 1 123456 NULL NULL
%>perl GLMGE_v4.pl input_matrix second_result/result.txt ID 1 123456 first_result/result.txt.outlier_99per.txt NULL
%>cat first_result/result.txt.outlier_99per.txt second_result/result.txt.outlier_99per.txt > outliers.txt
%>perl GLMGE_v4.pl input_matrix third_result/result.txt ID 1 123456 outliers.txt NULL
```
<img src="https://github.com/Park-Sung-Joon/GLMGE/assets/52985953/779fb0b2-863d-4341-b6f0-90e64f41cc86" width=50%>

Refer to the file **"example1/example1.run.sh"** for details.

### Example2: multiple runs with different random seeds
The "GLMGE_v4.pl" accepts a random seed. Using different seeds, you can run several times to get an ensemble of regression coefficients.
```
%>perl GLMGE_v4.pl input_matrix outdir/run_1/result.txt ID 1 123456 NULL NULL
%>perl GLMGE_v4.pl input_matrix outdir/run_2/result.txt ID 1 654321 NULL NULL
%>perl GLMGE_v4.pl input_matrix outdir/run_3/result.txt ID 1 9876543 NULL NULL
```
Then, the following command line gets merged results and stats; 
```
%>perl 2.analysis_v4.pl input_matrix outdir Merge_outdir
```
<img src="https://github.com/Park-Sung-Joon/GLMGE/assets/52985953/27596022-1e8b-4f5d-a7e4-df9a86de9b0c" width=50%>
+ (A) Pearson's correlation coefficient (R) between Observation and mean of Prediction of the final regression models
+ (B) Regression coefficient (RCs) of features after doing a two-tailed t-test and BH correction p-value (<0.05 adjusted p-value)
+ (C) Red points representing outliers in the distribution of 90%, 95%, and 99% points.

In addition, you can remove outliers and run again the regression;
```
%>perl GLMGE_v4.pl input_matrix first_outdir/run_1/result.txt ID 1 123456 NULL NULL
%>perl GLMGE_v4.pl input_matrix first_outdir/run_2/result.txt ID 1 654321 NULL NULL
%>perl GLMGE_v4.pl input_matrix first_outdir/run_3/result.txt ID 1 9876543 NULL NULL
%>perl 2.analysis_v4.pl input_matrix first_outdir Merge_first_outdir

%>perl GLMGE_v4.pl input_matrix second_outdir/run_1/result.txt ID 1 123456 Merge_first_outdir/Residuals/residual_stat.txt.outlier_99per.txt NULL
%>perl GLMGE_v4.pl input_matrix second_outdir/run_2/result.txt ID 1 654321 Merge_first_outdir/Residuals/residual_stat.txt.outlier_99per.txt NULL
%>perl GLMGE_v4.pl input_matrix second_outdir/run_3/result.txt ID 1 9876543 Merge_first_outdir/Residuals/residual_stat.txt.outlier_99per.txt NULL
%>perl 2.analysis_v4.pl input_matrix second_outdir Merge_second_outdir

%>cat Merge_first_outdir/Residuals/residual_stat.txt.outlier_99per.txt Merge_second_outdir/Residuals/residual_stat.txt.outlier_99per.txt > skipped_genes_in_3rd.txt
%>perl GLMGE_v4.pl input_matrix third_outdir/run_1/result.txt ID 1 123456 skipped_genes_in_3rd.txt NULL
%>perl GLMGE_v4.pl input_matrix third_outdir/run_2/result.txt ID 1 654321 skipped_genes_in_3rd.txt NULL
%>perl GLMGE_v4.pl input_matrix third_outdir/run_3/result.txt ID 1 9876543 skipped_genes_in_3rd.txt NULL
%>perl 2.analysis_v4.pl input_matrix third_outdir Merge_third_outdir
```

Refer to the file **"example2/example2.run.sh"** for details.

### What GE does?

LINE: Large-scale information network embedding (https://github.com/tangjianpku/LINE)
