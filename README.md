# GLMGE
## The prediction of gene expression by generalized linear regression model (GLM) and The prediction of co-regulated genes by graph embedding (GE)
![image](https://github.com/Park-Sung-Joon/GLMGE/assets/52985953/afedbb99-cfcc-4564-b751-1b0d1b215bf0)

### What GLM does?
The input (TXT file) is separated into dependent variables (Y) and explanatory variables (X). Y = wX + e

Full model ⇒ 
AIC (Akaike's Information Criterion) predicts the most adequate explanatory variables from a full regression model, which yields a reduced model M1 whose predictive power is evaluated by the 5-fold cross-validation (CV). Next, the M1 is extended using a variable randomly selected from the set of variables that AIC removed in the previous step, which yields a trial model M2. The AIC model selection and the 5-fold CV evaluate M2, and the current M1 is replaced with M2 if the predictive power has been improved. Otherwise, another M2 with a different random variable is tested. This procedure is repeated until all the variables removed have been tested. We ran one set of this iterative feature selection procedure 10 times by specifying different random seeds and performed two-tailed t-tests with the ensemble of estimated regression coefficients. The Bonferroni correction was applied, and variables having < 0.05 were considered significant.


```
%>perl GLMGE_v4.pl [INPUT MATRIX] [OUTPUT Filename] [ID] [Log10? 0|1] [RND_SEED] [SKIPgene_Filename|NULL] [SKIPcols: NULL | "A,B,C..."

%>perl GLMGE_v4.pl matrix/test.matrix.txt first_result/result.txt TEST 1 123456 NULL NULL
```
