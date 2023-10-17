# GLMGE
## The prediction of gene expression by generalized linear regression model (GLM) and the prediction of co-regulated genes by graph embedding (GE)
![image](https://github.com/Park-Sung-Joon/GLMGE/assets/52985953/afedbb99-cfcc-4564-b751-1b0d1b215bf0)

### What GLM does?
The input (TXT file) is separated into dependent variables (Y) and explanatory variables (X) and the regression model **Y = wX + e** is built.
1. Starting with a Full model (FM); X is all the explanatory variables in the input file.
2. AIC (Akaike's Information Criterion) reduces FM, which yields a reduced model **M1** and the set of removed features **XX** from the X, and N-fold cross-validation (CV) measures the predictive performance of **M1**.
  
3. **XX** is randomly shuffled. (searching the first-order effects)
4. Loop1: **M1** is extended with each of **XX**, which yields a trial model **M2**.
5. AIC and N-fold CV evaluate **M2**, which yields a new removed feature set **XX2**.
6. If **M2** improves **M1**, **M1** and **XX** are replaced with **M2** and **XX2**. Then, goto Step 3. :point_up_2:
7. Otherwise, goto Step 4. till all **XX** are tried. :point_up_2: 
8. End of Loop1

9.**XX** is randomly shuffled. (searching the second-order effect)
10. Loop2: **M1** is extended with each of all the possible pairs in **XX**, which yields a trial model **M2**.
12. AIC and N-fold CV evaluate **M2**, which yields a new removed feature set **XX2**.
13. If **M2** improves **M1**, **M1** and **XX** are replaced with **M2** and **XX2**. Then, goto Step 3.:point_up_2:
14. Otherwise, goto Step 10. till all pairs of **XX** are tried. :point_up_2:
15. End of Loop2

```
%>perl GLMGE_v4.pl [INPUT MATRIX] [OUTPUT] [ID] [Log10? 0|1] [RND_SEED] [SKIPgene_Filename|NULL] [SKIPcols: NULL|"A,B,C..."]
%>perl GLMGE_v4.pl matrix/test.matrix.txt first_result/result.txt TEST 1 123456 NULL NULL
```
