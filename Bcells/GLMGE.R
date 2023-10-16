
### tested with R 3.6
### source("GLMGE.R")x

### STEP 1. read data from inputfile
## see example files to know the file format: gene expression feature_1 feature_2 ....
## colname    <- c("1_colname", "2_colname", "3_colname") # IDs for explanatory variables
## SKIP_GENES <- c("gene1", "gene2", "gene3", ...)
## SKIP_COLS  <- c("10_colname", "11_colname", ...) # IDs of explanatory variables for skipping
## logmode <- 1 or 0. If one, expression is converted to log10(Expression+1)
## shuff <-  1 or 0. the original Y - X pairs are randomly changed
## Model <- RG_ReadDB("score.matrix.txt", colname, 1, 0, NULL, NULL)
## Model$X, Model$Y, Model$D
RG_ReadDB <- function(inputfile, colname, logmode, shuff, SKIP_GENES, SKIP_COLS){

	#--------------------------------------#  
	# 1. READ the input file  
	# no header line
	# columns are seperated by TAB or SPACE
	# the first col. is gene names (used as row ID)
	D <- read.table(inputfile, row.names=1, header=FALSE)
	# in this moment, col= expression_col + feature_cols
	cat( "# DIM(input) loaded: row = ", dim(D)[1], " col(expression+features) = ", dim(D)[2], "\n" )

	# Update D by skipping genes
	if(!is.null(SKIP_GENES)){
	    for(i in 1:length(SKIP_GENES)){
	    	  idx <- rownames(D) == SKIP_GENES[i]
		  D   <- D[!idx,]
            }
	}
	cat( "# DIM(input) removing genes: row = ", dim(D)[1], " col(expression+features) = ", dim(D)[2], "\n" )
	# replace NaN with 0
	D[ is.na(D) ] <- 0
	#--------------------------------------#  

	#--------------------------------------#  
	# 2. Set up explanatory variables (X)
	X <- as.matrix( D[,2:ncol(D)] ) # without expression column
	if( ncol(X) != length(colname) ){
	    cat("# Length of the feature_col in input matirx was " , ncol(X), "\n" )
	    cat("# Length of the given colname was " , length(colname), "\n" )
	    cat("# Aborted. GLMGE.R (1)", "\n")
	    #cat( colnames(X), "\n")
	    #cat( colname, "\n")
	    return()	
	}
	rownames(X) <- rownames(D)
	colnames(X) <- c(colname)

	# Skip explanatory variables. Use colnames
	if(!is.null(SKIP_COLS)){
	   for(i in 1:length(SKIP_COLS)){
	    idx <- colnames(X) == SKIP_COLS[i]
	    X   <- as.matrix(X[,!idx])
	   }
        }  
	cat( "# After removing given genes and features...", "\n" )
	cat( "# Now, the input matrix: Y = ", dim(X)[1], " X = ", dim(X)[2], "\n" )
        cn <- colnames(X)

	# removed all zero ROWs
	T  <- apply( abs(X), 1, sum )
	X2 <- c()
	rn <- c()
	for(i in 1:length(T)){
	  if(T[i] > 0){
             X2 <- rbind(X2, X[names(T[i]),])
             rn <- c(rn, names(T[i]) ) # retained rownames
          }
        }	
	X <- X2 # update X
	colnames(X) <- cn
        rownames(X) <- rn
	remove(X2)

	# removed all zoro COLs
	X2 <- c(); cn <- c()
	T <- apply(abs(X), 2, sum)
	for(i in 1:length(T)){
	  if(T[i] > 0){
             X2 <- cbind(X2, X[,names(T[i])])
             cn <- c(cn, names(T[i]) ) # retained colnames
          }
        }	
	remove(X2)

	# retained rows
	X <- as.matrix( D[rn,2:ncol(D)] )
	colnames(X) <- c(colname)
	# retained cols
        X <- as.matrix(X[,cn])
	rownames(X) <- rn
	colnames(X) <- cn
	cat( "# After removing all zero cases...", "\n" )
	cat( "# Size of X (features): ", dim(X)[2], "\n" )
	#--------------------------------------#

	#--------------------------------------#
	# 3. Set up Y (name is ACT)
	IDEXPRE <- "ACT"
	Y <- as.matrix(D[rn,1])
	rownames(Y) <- rn
	colnames(Y) <- IDEXPRE
	cat( "# Size of Y (gene expression) = ", dim(Y)[1], "\n" )
	#--------------------------------------#

	#--------------------------------------#  
	# 4. Log10 scale ?
	if(logmode){
         #  Y[,1] <- log10(Y[,1])
	 Y[,1] <- log10(Y[,1]+1)
        }
	cat( "# Y in log-10 scale: ", logmode, "\n" )
	#--------------------------------------#  

	#--------------------------------------#  
	# 5. Shuffling X for Y
	if(shuff){
	   idx <- sample(1:nrow(X))
	   X2  <- X[idx,]
	   X   <- X2         # shuffled X
	   #========= 
	   rownames(X) <- rn # put original Y names for making random matrix
 	   colnames(X) <- cn # put original X names for making random matrix
	   #=========
	   remove(X2)
	}
	cat( "# Random shuffling of X: ", shuff, "\n" )
	#--------------------------------------#

	#--------------------------------------#  
	# 6. Shuffling for columns (X), explanatory variables
	# Y - X pair will not be changed.
	# This suffling is performed every time. (based on random seed)
	idx <- sample(1:ncol(X))
	X2 <- X[,idx]
	cat( "# X was shuffled for this run. Note, Y-X pairs were not changed", "\n" )
	#--------------------------------------#  

	#--------------------------------------#  
	# 7. Normalization (limma)
	library(limma)
	X <- normalizeQuantiles(X2, ties=TRUE)
	cat( "# X was normalized", "\n" )
	#--------------------------------------#  
	remove(X2)
	
        ### return the results. $D: imported matrix, Y: expression ("ACT"), X: features, rn: gene names, cn: feature names
	cat("# DIM(D) = ", dim(D), "\n") # given genes removed, expression + features
	cat("# Y = ", dim(Y)[1], "\n")
	cat("# X = ", dim(X)[2], "\n")
	cat("# Length of X and feature names were   " , ncol(X), " ", length(cn), "\n" )
	cat("# Length of Y and number of genes were " , nrow(Y), " ", length(rn), "\n" )

	if( ncol(X) != length(cn) ){
	    cat("# Length of X and feature names were " , ncol(X), " ", length(cn), "\n" )
	    cat("# Aborted. GLMGE.R (2)", "\n")
	    return()	
	}
	if( nrow(Y) != length(rn) ){
	    cat("# Length of Y and number of genes were " , nrow(Y), " ", length(rn), "\n" )
	    cat("# Aborted. GLMGE.R (3)", "\n")
	    return()	
	}

	return(	list(D=D, Y=Y, X=X, rn=rn, cn=cn) )
}



### STEP 2. Linear regression
## Y is the gene expression predicted
## X is the explanatory variables
RG_Regression <- function(Y, X){

	debug <- 0
	if(debug){
		cat("0. Y = ", dim(Y), "\n")
		cat("0. X = ", dim(X), "\n")
	}
	
	# 1. Regression Model (full model)
	RM <- lm( Y~ ., data=as.data.frame(X), singular.ok=TRUE)
	
	# NOTE: if there exist "X:B" in colnames(=feature names), this line causes error. Also "Nkx3-1" -> "`Nkx3-1`"
	# So, watch out !!!!
	# 2-1. remove singular "NA"
	# "cn" is the array with feature names and Residuals
        cn  <- attributes(anova(RM))$row.names
        # 2-2. remove "Residuals"
        cn  <- cn[1: length(cn)-1]
        # 2-3. remove backquote
        cn  <- gsub("`", "", cn)
	remove(RM)

	# 3. update X (this step is for removing abnormal features as done in step 2-1 ~ 2-3.
	X  <- as.matrix( X[ ,cn] )
	if(debug){
		cat("3. Y = ", dim(Y), "\n")
		cat("3. X = ", dim(X), "\n")
        }

# Normalization (not used)
#        p      <- ncol(X) #predictors used in this CV
#        bN     <- nrow(X) #genes used in this CV
#        one    <- matrix(1,1,bN)
#        mu     <- one %*% Y/bN  
#        meanx  <- drop(one %*% X)/bN
#        sY     <- scale(Y, drop(mu), FALSE)
#        sX     <- scale(X, meanx, FALSE)
#        normx  <- sqrt(drop(one %*% (sX^2))/(bN - 1))
#        sX     <- scale(sX, FALSE, normx)

	# 4-1. keep the original Y and X 
	sY <- Y; sX <- X;
	# 4-2. Run regression with full model of new X (=sX)
	RM   <- lm( sY~ ., data=as.data.frame(sX), singular.ok=TRUE)
        pred <- predict(RM) # vectors

	# for making consistency between prediction and observation
        Y2   <- as.matrix( sY[ names(pred), ]) # pick up predicted genes only
	X2   <- as.matrix( sX[ names(pred), ]) # pick up features used for this prediction
	rownames(Y2) <- names(pred)
	rownames(X2) <- names(pred)
	colnames(Y2) <- colnames(sY)
	colnames(X2) <- colnames(sX)

	# update the given X and Y (only for clarify)
	X <- X2; Y <- Y2
	if(debug){
		cat("4. Y = ", dim(Y), "\n")
		cat("4. X = ", dim(X), "\n")
	}
	remove(sY); remove(sX)

	OX <- X; OY <- Y; #not used

	# 5-1. run again (only for clarify)
	RM   <- lm( Y~. , data=as.data.frame(X), singular.ok=TRUE)
	pred <- predict(RM)

	# performance
	C    <- as.vector(cor(Y, pred))
	if(debug){
		cat("prediction correlation: ", C, "\n")
	}
	
	### return the results
	#  Y, X: the given matrix (original matrix)
	#  OY, OX: the matrix used in this regression. not used (Y is equal to OY, X is equal to OX)
	#  PRED: prediction of gene expression
	#  RM: reduced model
	#  C: performance (correlaton with Prediction and Observation)
	return( list(Y=Y, X=X, OY=Y, OX=X, C=C, PRED=pred, RM=RM) )
}


# look at the "RG_Regression" for detailes. Here use rlm instead of lm
RG_Regression_Rlm <- function(Y, X){

	RM <- rlm( Y~ ., data=as.data.frame(X), method="MM")
	
	#remove singular NA
        cn  <- attributes(anova(RM))$row.names
        #cn inlucdes 1 TFBSxx ..... last Residuls

        # remove "Residuals"
        cn  <- cn[1: length(cn)-1]
        # remove backquote
        cn  <- gsub("`", "", cn)
	remove(RM)

	# If there exist "X:B", this line shows error
	# Also "Nkx3-1" -> "`Nkx3-1`"	
	# So, how to remove quote!
	X  <- as.matrix( X[ ,cn] )
        
	sY <- Y; sX <- X;
	RM   <- rlm( sY~ ., data=as.data.frame(sX), method="MM", singular.ok=TRUE)
        pred <- predict(RM)
        Y2   <- as.matrix( sY[ names(pred), ])
	X2   <- as.matrix( sX[ names(pred), ])
	rownames(Y2) <- names(pred)
	rownames(X2) <- names(pred)
	colnames(Y2) <- colnames(sY)
	colnames(X2) <- colnames(sX)

	X<-X2; Y<-Y2;
	 
	OX <- X; OY <- Y;
	# run again
	RM   <- rlm( Y~. , data=as.data.frame(X), method="MM", singular.ok=TRUE)
	pred <- predict(RM)
	C    <- as.vector(cor(Y, pred))

	### return the results
	return( list(Y=Y, X=X, OY=Y, OX=X, C=C, PRED=pred, RM=RM) )
}



### STEP 3. Linear regression for finding reduced model
## Model must be the output from "RG_Regression()"
RG_ReduceModel <- function(Model){
	debug <- 0
	if(debug){ cat("In RG_ReduceModel", "\n") }
	       
	X <- Model$X # input features
	Y <- Model$Y # input expression predicted

	if(debug){ cat("0. Y = ", dim(Y), "\n") }
	if(debug){ cat("0. X = ", dim(X), "\n") }
	
	#----------------------#
	# 1. Regression Model with the input current model (= full model)
	RM  <- lm( Y~. , data=as.data.frame(X), singular.ok=TRUE)
	all_features <- colnames(RM$model)[2:ncol(RM$model)]
	#----------------------#

	#----------------------#
	# 2.Reduction by AIC (This takes time)
	# the number of degrees of freedom (k=2). if k=log(n), this is BIC
	library(MASS)
	#AIC <- stepAIC(RM, trace=0, direction="both")

	# speedy version
	AIC <- stepAIC(RM, trace=0, direction="backward", steps=10)

	# AICed Regression Model
	AM <- AIC$model

	# Final AIC score
	SCORE <- AIC$anova$AIC[ length(AIC$anova$AIC) ]
	if(debug){ cat( SCORE, "\n") }

	#----------------------#
	# 3. Find retained features after AIC
	selected_features <- colnames(AM)[ 2:ncol(AM) ]
	if(debug){ cat( "selected ", selected_features , "\n") }

	# 4. Find removed features after AIC
	removed_features    <- c()
	for( i in 1:length(all_features) ){
	  f <- 1
	  for( j in 1:length(selected_features) ){
	   if( all_features[i] == selected_features[j] ){
	     f <- 0; break
	   }
	  }

	  if(f){
	    removed_features <- c(removed_features, all_features[i])
	  }
	}
	
	if(debug){ cat("removed ", removed_features, "\n") }
	#----------------------#

	### return the results
	# AIC.     Information of AIC reduction
	# SCORE.   The final AIC score
	# TERM.    Retained features by AIC
	# REMOVED. Removed  features after AIC
	return( list(AIC=AIC, SCORE=SCORE, TERM=selected_features, REMOVED=removed_features) )
}


### STEP 4. Cross Validation
cv.park <- function (x, y, fold = 10, run=10){
    x  <- as.matrix(x)
    n  <- nrow(x)
    y  <- as.matrix(y)

    MAX_TC <- 0.0  # maximal testing correlation
    MAX_XY <- 0.0  

    for(i in 1:run){ # iterative run of x-fold CVs

    # rank-deficient fit: collinear
    TOTAL_TC <- 0.0

    foldi <- split( sample(1:n), rep(1:fold, length=n) )

    for (j in 1:fold) {    # x Blocks
        omit <- foldi[[j]] # Testing block

	# Learning
 	bX   <- as.matrix( x[-omit,] ) # features in leaning

	# Find all zero vectors. If all zero column exists, remove this column
	zC   <- vector()
	for(k in 1:ncol(bX)){
	  if( sum(bX[,k]) == 0 ){ zC <- cbind(zC,k) }
	}

        # Normalization
        bY   <- as.matrix(y[-omit,]) # gene expression
        bX   <- as.matrix(x[-omit,]) # features in learning
        newX <- as.matrix(x[omit, ]) # features in testing
        if(length(zC)>0){            # remove all zero columns
        	bX   <- as.matrix(x[-omit,-zC]) # for leaning
	        newX <- as.matrix(x[omit,-zC])  # for testing 
        }

        p     <- ncol(bX) # predictors used in this CV
        bN    <- nrow(bX) # genes used in this CV

        one   <- matrix(1,1,bN)
        mu    <- one %*% bY/bN	
        meanx <- drop(one %*% bX)/bN

	sY     <- scale(bY, drop(mu), FALSE)
        sX     <- scale(bX, meanx, FALSE)
        normx  <- sqrt(drop(one %*% (sX^2))/(bN - 1))
        sX     <- scale(sX, FALSE, normx)
        newX   <- scale(newX, meanx, normx)
        #=============================================#

	#===========================================#
	# Learning using x-1 blocks
	RG <- RG_Regression(sY, sX)
	R  <- RG$RM
	LC <- cor(sY, predict(R))          # learning correlation
	#R2 <- summary(R)$adj.r.squared    # learining R^2
	#===========================================#
       
	#===========================================#
        # Testing
	pred <- predict.lm(R, as.data.frame(newX) )  # testing. Prediction using one block
        pred <- pred + matrix(1, nrow(newX),1) %*% mu
	#===========================================#

	#===========================================#
	# R^2 = sigma(pred - mean)^2 / sigma(obs  - mean)^2
	m  <- mean(y[omit,])
	R2 <- sum( (pred - m)^2 ) / sum( (y[omit,] - m)^2 ) #raw R^2
	#===========================================#

	#===========================================#
	# Adjusted R^2
	p <- nrow(newX) #samples
	q <- ncol(newX) #predictors
	if(q != (p -1)){
		AR2 <- 1.0 - ( (p-1)/(p-q-1) * (1.0-R2) )
	}else{
		AR2 <- 0.0
	}
	#===========================================#

        #-----------------------------# 
	TC       <- cor(y[omit,], pred) # testing correlation
        TOTAL_TC <- TOTAL_TC + TC       # total of testing correlation
        #-----------------------------# 

        #-----------------------------# 
        TWO <- cbind(y[omit,], pred)
        if(j == 1){
           XY <- TWO
        }else{
           XY <- rbind(XY,TWO)
        }
        #-----------------------------# 


	if(i == 1 && j == 1){
	  x_total        <- ncol(newX)
	  y_total        <- nrow(newX)
	  learning_total <- LC
	  testing_total  <- TC
	  rsquare        <- R2
	  arsquare       <- AR2
	}else{
	  x_total        <- cbind(x_total, ncol(newX))
	  y_total        <- cbind(y_total, nrow(newX))
	  learning_total <- cbind(learning_total, LC)
	  testing_total  <- cbind(testing_total, TC)
	  rsquare        <- cbind(rsquare, R2)
	  arsquare       <- cbind(arsquare, AR2)
	}
   }
   #end of this CV	

   #--------------------------------#
   # for the case of xruns of CVs
   # take most significant run
   TOTAL_TC <- TOTAL_TC / fold
   if(MAX_TC < TOTAL_TC){
        MAX_TC <- TOTAL_TC  #maximal testing correation coefficient
        MAX_XY <- XY        #XY data (obs. vs. pred)
   }
   #--------------------------------#

  }

   xavg   <- mean(x_total)
   yavg   <- mean(y_total)
   lavg   <- mean(learning_total)
   lsd    <- sd(as.vector(learning_total))	
   tavg   <- mean(testing_total)
   tsd    <- sd(as.vector(testing_total))	
   r2avg  <- mean(rsquare)
   r2sd   <- sd(as.vector(rsquare) )

   ar2avg  <- mean(arsquare)
   ar2sd   <- sd(as.vector(arsquare) )
	
   ## return the results
   return( list(LC=lavg,LSD=lsd,TC=tavg,TSD=tsd, R2=r2avg,R2SD=r2sd, AR2=ar2avg,AR2SD=ar2sd, RUN=run, nY=yavg, nX=xavg, XY=MAX_XY, XYC=MAX_TC) )
}




#============================================================#
### STEP 5. Compare two models and update
UpdateModel <- function(BESTR, BESTTCV, TR, ID, IDT, FOLD, keta){

    temp1 <- round(TR$C, keta)
    temp2 <- round(BESTR$C, keta)

    flag <- 0
    if( temp1 >= temp2){ flag = 1 } #improved
	    
    # Check the improvement by x-Fold cross-validation 
    if(flag){
  	variables <- colnames(TR$RM$model)[2:length(colnames(TR$RM$model))]
	names(variables) <- variables

	if( !is.na(variables[IDT]) ){ #this model includes this valiable
	    # Do CV at here
	    TCV   <- cv.park(TR$OX, TR$OY, FOLD, 1)

   	    temp3 <- round( TCV$TC, keta )
	    temp4 <- round( BESTTCV$TC, keta )
	    # This is very important. 
	    # If TC was improved, anyway accept it
	    # Even if TC was tie, C was imporoved, then accept it
	    if( (temp1 == temp2) && (temp3 <= temp4)){
		return( list(BESTR=NULL, BESTTCV=NULL, BESTID=NULL, BESTIDT=NULL,updated=0) )
		#next
	    } #nothing to do
		    
	    # Update the BEST model (current model)
	    return( list(BESTR=TR, BESTTCV=TCV, BESTID=ID, BESTIDT=IDT, updated=1) )
       }
    }

    # end of if( !is.na(variables[IDT]) )
    return( list(BESTR=NULL, BESTTCV=NULL, BESTID=NULL, BESTIDT=NULL,updated=0) )
}
#============================================================#
