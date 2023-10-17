# GLMGE
## The prediction of gene expression by generalized linear regression model (GLM) and the prediction of co-regulated genes by graph embedding (GE)
![image](https://github.com/Park-Sung-Joon/GLMGE/assets/52985953/afedbb99-cfcc-4564-b751-1b0d1b215bf0)

### What GLM does?
The input (TXT file) is separated into dependent variables (Y) and explanatory variables (X) and the regression model **Y = wX + e** is built.
1. Starting with a Full model (FM); X is all the explanatory variables in the input file.
2. AIC (Akaike's Information Criterion) reduces FM, which yields a reduced model **M1** and the set of removed features **XX** from the X, and N-fold cross-validation (CV) measures the predictive performance of **M1**.
3-1. **XX** is randomly shuffled.
4-2. **M1** is extended using a variable randomly selected from **XX**, which yields a trial model **M2**.
5. AIC and N-fold CV evaluate **M2**, which yields a new removed feature set **XX2**.
6. If **M2** improves **M1**, **M1** and **XX** are replaced with **M2** and **XX2**. Goto Step 3.:point_up_2:
7. Otherwise, **M1** is extended with each of all the possible pairs in **XX2**, which yields a trial model **M2**.
8. AIC and N-fold CV evaluate **M2**, which yields a new removed feature set **XX2**.
9. 
10. a variable is if the predictive power has been improved.
11. ested. We ran one set of this iterative feature selection procedure 10 times by specifying different random seeds and performed two-tailed t-tests with the ensemble of estimated regression coefficients. The Bonferroni correction was applied, and variables having < 0.05 were considered significant.


```
%>perl GLMGE_v4.pl [INPUT MATRIX] [OUTPUT Filename] [ID] [Log10? 0|1] [RND_SEED] [SKIPgene_Filename|NULL] [SKIPcols: NULL | "A,B,C..."

%>perl GLMGE_v4.pl matrix/test.matrix.txt first_result/result.txt TEST 1 123456 NULL NULL
```
