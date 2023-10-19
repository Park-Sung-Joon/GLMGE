
library(MASS)
library(limma)
library(coefplot)
library(cowplot)

# load GLMGE functions
source("../src/GLMGE.R")

# set random seed
set.seed(123456, kind="Mersenne-Twister")

# Result file
out <- "fourth_result/result.txt"

# Genes and Features to be removed
GIVEN_SKIP_GENES <- NULL
SKIP_COLS <- NULL
if( 0 == 1){
    SKIP_COLS <- c(NULL) # "AAA,CCC,BBB"
}

# Gene list given:  one gene per one line
if( "skipped_genes_in_4th.txt" != "NULL" ){
    # read the file
    tmp <- read.table("skipped_genes_in_4th.txt", header=FALSE, as.is=TRUE)
    # use the first column only (= gene name)
    GIVEN_SKIP_GENES <- unique( tmp[,1] )
    remove(tmp)
}

#===================================================#
# Step 0: PREPARE INPUT DATA, not normalized
#         BR (before AIC reduction)
#===================================================#
coln  <- c("TFBS1","TFBS2","TFBS3","TFBS4","TFBS5","TFBS6","TFBS7","TFBS8","TFBS9","TFBS10","TFBS11","TFBS12","TFBS13","TFBS14","TFBS15","TFBS16","TFBS17","TFBS18","TFBS19","TFBS20","TFBS21","TFBS22","TFBS23","TFBS24","TFBS25","TFBS26","TFBS27","TFBS28","TFBS29","TFBS30","TFBS31","TFBS32","TFBS33","TFBS34","TFBS35","TFBS36","TFBS37","TFBS38","TFBS39","TFBS40","TFBS41","TFBS42","TFBS43","TFBS44","TFBS45","TFBS46","TFBS47","TFBS48","TFBS49","TFBS50","TFBS51","TFBS52","TFBS53","TFBS54","TFBS55","TFBS56","TFBS57","TFBS58","TFBS59","TFBS60","TFBS61","TFBS62","TFBS63","TFBS64","TFBS65","TFBS66","TFBS67","TFBS68","TFBS69","TFBS70","TFBS71","TFBS72","TFBS73","TFBS74","TFBS75","TFBS76","TFBS77","TFBS78","TFBS79","TFBS80","TFBS81","TFBS82","TFBS83","TFBS84","TFBS85","TFBS86","TFBS87","TFBS88","TFBS89","TFBS90","TFBS91","TFBS92","TFBS93","TFBS94","TFBS95","TFBS96","TFBS97","TFBS98","TFBS99","TFBS100") #feature names

FOLD  <- 5           #X-fold CVs
ifelse(1 == 1, logMode <- T, logMode <- F) #expression in log-10 scale


sink(out, append=FALSE)
cat("# INPUT FILE: " , "../src/matrix/example.matrix.txt", "\n")
cat("# OUT FILE: ", out, "\n")
cat("# RND SEED: ", "123456", "\n")
cat("# GIVEN_SKIP_GENES: ", length(GIVEN_SKIP_GENES), "\n")
cat("# GIVEN_SKIP_COLS:  ", length(SKIP_COLS), "\n") 

cat("# Number of Features: ", length(coln), "\n")
cat("# x-FOLD CVs: ", FOLD, "\n")
cat("# log-10 scale: ", "1", "\n")
sink()

cnt_try <- 0
repeat{
  # "DB" keeps the original input data (X is randomly shuffled for this run, but Y-X pairs are retained)
  # NOTE: If colname includes "-" and other, anova automatically adds backquote `XXX` in some case.
  #       So, simpler name of features will be better.  
  DB <- RG_ReadDB("../src/matrix/example.matrix.txt", coln, logMode, 0, GIVEN_SKIP_GENES, SKIP_COLS)

  length(coln)
  nrow(DB$D); ncol(DB$D)
  nrow(DB$Y); ncol(DB$Y)
  nrow(DB$X); ncol(DB$X)

  # Before AIC selection. This is the first look of this input file.
  BR      <- RG_Regression(DB$Y, DB$X) 
  cnt_try <- cnt_try + 1

  if(cnt_try > 10){ break } #Stop
  if(BR$C > 0.1){ break }  #Stop

  # if the BR model looks bad,,, run "RG_ReadDB" again.
  # So, the function shuffles X again, which results differentially.
  # run <=10, if the low C is still .... this GLMGE will not be meaningful.
}
#===================================================#


#===================================================#
sink(out, append=TRUE)
cat("# Step 1: Reduce BR model", "\n");
sink()
#         AR (after AIC reduction)
#===================================================#
cnt_try <- 0
repeat{

   # AIC selection, BR includes all TFs (excepting singular TFs)
   M <- RG_ReduceModel(BR)         #model 

   # Build X from the reduced model M
   X <- as.matrix(DB$X[,M$TERM]) #retained features in M and their orig data in DB
   colnames(X) <- M$TERM
   rownames(X) <- rownames(DB$Y)  #orig data

   # Y from the reduced model will be DB
   # Now, Regression again with X and Y of the reduced model.
   # This is just for clerify
   AR <- RG_Regression(DB$Y, X)
   # Then, run x-fold Cross-validationsl (OX and OY will be equal to X and Y).
   #CV <- cv.park(AR$X, AR$Y, FOLD, 1)
   CV <- cv.park(AR$OX, AR$OY, FOLD, 1)

   cnt_try <- cnt_try + 1
   if(cnt_try > 10){ break }
   if(AR$C > 0.1){ break }
   
   # Try again due to something wrong happened
   # read the input file by shuffling X
   DB <- RG_ReadDB("../src/matrix/example.matrix.txt", coln, logMode, 0, GIVEN_SKIP_GENES, SKIP_COLS)
   # and generate BR model       
   BR <- RG_Regression(DB$Y, DB$X) 
}

### PRED is vector, Y is matrix, X is matrix
### names(AR$PRED)
### rownames(AR$Y)

#===================================================#
sink(out, append=TRUE)
cat("# Step 2: Loss and Find. Back the removed feature to AR model", "\n")
sink()
#===================================================#

# Features removed by AIC reduction
removed  <- M$REMOVED
# number of the removed features. a stop condition
stop <- length(removed)

X  <- AR$OX
Y  <- AR$OY
rn <- rownames(AR$OX)
cn <- colnames(AR$OX)

# the 1st "Loss and Find" procedure
for( i in 1:stop){
   NX  <- as.matrix(DB$X[,removed[i]]) #this removed feature vector for all the genes
   ID  <- removed[i]
   rownames(NX) <- rownames(DB$X) # rowName = Gene IDs
   NX  <- NX[rn,]                  # pick up genes presented at AR model

   TX <- cbind(X, NX)              # extended model with this removed feature 
   colnames(TX) <- c(cn, ID)       # add this colname 
   rownames(TX) <- rn

   # Regression -> AIC
   TR <- RG_Regression(Y, TX)  # TX is the trial model with this removed feature
   TM <- RG_ReduceModel(TR)    # Reduce TR model again by AIC

   TX2 <- TX[, TM$TERM]
   colnames(TX2) <- TM$TERM
   rownames(TX2) <- rn
   TR  <- RG_Regression(Y, TX2)            # Run again for clarify
   TCV <- cv.park(TR$OX, TR$OY, FOLD, 1) # X-fold CVs.

   # Check the CV correlation of Trial Model (TR) and Update the current model (AR)
   if(TR$C >= AR$C){
      if(TCV$TC > CV$TC){
         AR <- TR
         CV <- TCV
      }
  }

} ## for i 1 to STOP

# Output. BR model. Before AIC reduction
sink(out, append=TRUE)
summary(BR$RM)
sink()
write.table(paste("#Before AIC: ", "R= ", round(BR$C,2), " Y= ", nrow(BR$Y), " X= ", ncol(BR$X), sep=""), file=out,sep="",quote=FALSE,append=TRUE, row.names=FALSE, col.names=FALSE)

# Output. AR model. After AIC reduction. Appended mode
sink(out, append=TRUE)
summary(AR$RM)
sink()

write.table(paste("#After AIC: ", "R= ", round(AR$C,4), " Y= ", nrow(AR$Y), " X= ", ncol(AR$X), sep=""),
            file=out,sep="",quote=FALSE,append=TRUE, row.names=FALSE, col.names=FALSE)

write.table(paste("#After AIC CV: ", "LC= ", round(CV$LC,4), " TC= ", round(CV$TC,4), sep=""),
            file=out,sep="",quote=FALSE,append=TRUE, row.names=FALSE, col.names=FALSE)
#==  End of 1st loss and find procedure =======================#



# the 2nd "Loss and Find" procedure

# Find the removed features from AR model
variables <- colnames(AR$RM$model)[ 2:length(colnames(AR$RM$model)) ]
removed   <- c()
for(n in 1:length(DB$cn)){

    f <- 1
    for(m in 1:length(variables) ){
        #already exist in the model?
        if(DB$cn[n] == variables[m]){ f<-0; break }
    }

    if(f){ # not found from the current AR model

	#--------------------------------------------------
	# Comment out this, if the feature names have a pefix, e.g., "TFBS1", "TFBS2", "TFBS10".
	#--------------------------------------------------
	#if( length ( grep("TFBS[0-9]+", DB$cn[n]) ) > 0 ){
	#    removed <- c(removed, DB$cn[n])
	#}

	# default. pick up the removed features
	removed <- c( removed, DB$cn[n] )
    }
}
#------------------#

# Start configuration

#========== These are options. Shuffle promoter TFBSs and suffle lrc TFBSs
# 1. select intraTFBS (TFBS found from the long-range contacted regions)
#idx   <- grep("^intraTFBS", removed)
#INTRA <- c()
#if(length(idx)){
  # Random shuffling
#  INTRA <- sample( removed[idx], replace=FALSE)
#}

# 2. select not intraTFBS (= TFBS found from promoters)
#idx  <- grep("^intraTFBS", removed, invert=TRUE)
#TFBS <- c()
#if(length(idx)){
  # Random shuffling
#  TFBS <- sample( removed[idx], replace=FALSE)
#}

# 3. Merge two randomly shuffled TFBSs arrays
#removed <- c( TFBS, INTRA) #TFBS first
#================ end of the option

# Default. random shuffling of removed features
removed <- sample( removed, replace=FALSE )

#============================================================#
#============================================================#
## Iterative finding single or pairwise combinatorial features from the removed featurs
#============================================================#
keta    <- 3 # number of digits after decimal point 
TRIED   <- c(); ID_POOL <- c()
k       <-  1;  LOOP    <- 1
#============================================================#

while(LOOP){

    # default stop condition
    LOOP <- 0
    stop <- length(removed)

    updated <- 0 # reset
    for( i in 1:stop){
	BESTR   <- AR; BESTTCV <- CV
        BESTID  <- BESTIDT <- c()
	
	#TFs already tried are never going to try again
	TRIED   <- unique( c(TRIED, removed[i]) )

	# Set a model from the current AR model
	X  <- AR$OX; Y  <- AR$OY
	rn <- rownames(AR$OX)
	cn <- colnames(AR$OX)

	#-------- Generate new input matrix ---------------#
	# feature vector coming from input file
	NX <- as.matrix( DB$X[,removed[i]] )
		
	# orignal pair name and temporal name (not used)
	ID  <- removed[i]; IDT <- ID

	# select only target genes
	rownames(NX) <- rownames(DB$X)
	NX <- NX[rn,]

	# merge columns with a new column (=new feature)
	TX <- cbind(X, NX)
	colnames(TX) <- c(cn, IDT)
	rownames(TX) <- rn
	#--------- done --------#

	# Regression and AIC reduction
	TR <- RG_Regression(Y, TX)
	TM <- RG_ReduceModel(TR)

	TX2 <- as.matrix( TX[,TM$TERM] )
	colnames(TX2) <- TM$TERM
	rownames(TX2) <- rn
	TR  <- RG_Regression(Y, TX2)

	CHK <- UpdateModel(BESTR, BESTTCV, TR, ID, IDT, FOLD, keta)
	if(CHK$updated == 1){ # rarely updated 
	    updated <- 1  # do here
	    AR <- BESTR   <- CHK$BESTR
	    CV <- BESTTCV <- CHK$BESTTCV
	    BESTID  <- CHK$BESTID
	    BESTIDT <- CHK$BESTIDT

	    # Write information to output file
	    write.table( paste0("#  Update local: ", "R= ", round(TR$C,keta), " V= ", ID, " (", IDT, ")", 
			      " CV_LC= ", round(TCV$LC,keta), " CV_TC= ", round(TCV$TC,keta) ),
			 file=out,sep="",quote=FALSE,append=TRUE,
			 row.names=FALSE, col.names=FALSE
	    )
       }
       remove(CHK)
    } ## for i 1 to STOP

    if( updated == 0){ # not updated so far. try to find a pairwise one

     if(stop < 2){ break }

     # Here we go for finding a pairwise combination of features
     for( i in 1: (stop-1) ){
	# Information of the best model in this loop is saved to BESTxx
	BESTR   <- AR
        BESTTCV <- CV
        BESTID  <- BESTIDT <- c()
	
	#TFs already tried are never going to try again
	TRIED   <- unique( c(TRIED, removed[i]) )

	# Set a model from the current AR model
	X  <- AR$OX
	Y  <- AR$OY
	rn <- rownames(AR$OX)
	cn <- colnames(AR$OX)

	updated <- 0 # reset
	for( j in (i+1):stop){

	    #-------- Generate new input matrix ---------------#
	    # feature vector coming from input file
	    NX <- as.matrix(DB$X[,removed[i]] * DB$X[,removed[j]])
		
	    # orignal pair name
            ID  <- paste0( removed[i],":",removed[j] )
	    # temporal pair name
	    IDT <- paste0("PARKSJ", k)

	    # select only target genes
	    rownames(NX) <- rownames(DB$X)
	    NX <- NX[rn,]

	    # merge columns with a new column (=new feature)
	    TX <- cbind(X, NX)
	    colnames(TX) <- c(cn, IDT)
	    rownames(TX) <- rn
	    #--------- done --------#

	    # Regression and AIC reduction
	    TR <- RG_Regression(Y, TX)
	    TM <- RG_ReduceModel(TR)

	    TX2 <- as.matrix( TX[,TM$TERM] )
	    colnames(TX2) <- TM$TERM
	    rownames(TX2) <- rn
	    TR  <- RG_Regression(Y, TX2)

	    CHK <- UpdateModel(BESTR, BESTTCV, TR, ID, IDT, FOLD, keta)
	    #updated <- CHK$updated ### no, no, do not do like this
	    #if(updated == 1){
	    if(CHK$updated == 1){
		updated <- 1  # do here
		k = k + 1
	        BESTR   <- CHK$BESTR
		BESTTCV <- CHK$BESTTCV
		BESTID  <- CHK$BESTID
		BESTIDT <- CHK$BESTIDT

	        # Write information to output file
	        write.table( paste0("#  Update local: ", "R= ", round(TR$C,keta), " V= ", ID, " (", IDT, ")", 
				      " CV_LC= ", round(TCV$LC,keta), " CV_TC= ", round(TCV$TC,keta)),
				file=out,sep="",quote=FALSE,append=TRUE,
				row.names=FALSE, col.names=FALSE
		 )
	   }
	   remove(CHK)
     } # Next J. Try all Js even updated==1
     # End of the searching J for the feature I

	#=========================#
	# break I loop, and reset parameters and loop again	
	if( updated ){
	    break #!!important
	}
	#=========================#

    } # Next Feature i

    } #   if( updated == 0)

    # If, the model was updated, reset TF sets
    if(updated){

      if( length( grep(":", BESTID)) > 0 ){
         ID_POOL <- rbind(ID_POOL, BESTID)
         rownames(ID_POOL)[ nrow(ID_POOL)] <- BESTIDT
      }
      ID_POOL <- unique( ID_POOL )
 
      AR  <- BESTR
      CV  <- BESTTCV
      write.table(paste0("#Update: ", "R= ", round(AR$C,keta), " Y= ", nrow(AR$Y), " X= ", ncol(AR$X),
                       " Added= ", BESTID, " (", BESTIDT, ")",
                       " CV_LC= ", round(CV$LC,keta), " CV_TC= ", round(CV$TC,keta)),
                 file=out,sep="",quote=FALSE,append=TRUE, row.names=FALSE, col.names=FALSE)

     #----------------------------------#
     # Get variables used in this finding
     variables <- colnames(AR$RM$model)[2:length(colnames(AR$RM$model))]
     names(variables) <- variables

     TID_POOL <- c()            # clear ID pool
     this_selected <- variables # default
    
     idx <- grep("PARKSJ", variables)
     if(length(idx) > 0){ # if a new feature found

       this_selected <- variables[ -idx ] # selected feature names without "PARKKSJ"
    
       for(n in 1:length(idx)){
          M <- variables[ idx[n] ] # m must be in ID_POOL

          # Keep this information
          TID_POOL      <- rbind(TID_POOL, as.matrix(ID_POOL[M, ]))

          # Get each TF in this feature pair. update selected feature names
          this_selected <- c(this_selected, unlist( strsplit(ID_POOL[M,1], ":") )[1])
          this_selected <- c(this_selected, unlist( strsplit(ID_POOL[M,1], ":") )[2])
        } #next N

        ID_POOL <- TID_POOL #update ID POOL
     }
     names(this_selected) <- c()
     #------------------------------#

     # Update the featurs already tried
     this_selected <- unique( c(this_selected, TRIED))
   
     # Reset "removed" array. This is used in the next round.
     removed <- c()
     for(n in 1:length(DB$cn)){

	 f <- 1 # not found in "this_selected" array
	 for(M in 1:length(this_selected)){
	     if(DB$cn[n] == this_selected[M]){ f<-0; break }
	 }

	 if(f){
	     # This is an option
	     # update remained features. xxxxTFBSyyyy only
	     #if( length ( grep("TFBS[0-9]+", DB$cn[n]) ) > 0 ){
	     #	 removed <- c(removed, DB$cn[n])
	     #}
	     
	     # Defalt mode
	     removed <- c(removed, DB$cn[n])
	 }
     }
     #------------------#


     # If features still left,
     # Random shuffling again!!!! So, next loop use different order
     if(length(removed) > 1){

	 # This is options;
	 #idx <- grep("^intraTFBS", removed)
	 #INTRA <- c()
	 #if(length(idx)){
	 #    INTRA <- sample( removed[idx], replace=FALSE)
         #}

	 #idx <- grep("^intraTFBS", removed, invert=TRUE)
	 #TFBS <- c()
	 #if(length(idx)){
	 #    TFBS <- sample( removed[idx], replace=FALSE)
         #}
	 #removed <- c( TFBS, INTRA) #TFBS first
	 
	 # Default mode.
	 removed <- sample(removed, replace=FALSE)
         # Loop again
	 LOOP <- 1

	 # skip other Ith. back to while(LOOP)
	 #break #!!important
     }

    } #----------- End of Update -----------#
  
 

 # End ?
 if(LOOP == 0){ break }
}
#warnings()
write.table( paste0("# DONE: ", "R= ", round(AR$C,keta), " Y= ", nrow(AR$Y), " X= ", ncol(AR$X),
		  " CV_LC= ", round(CV$LC,keta), " CV_TC= ", round(CV$TC,keta)),
	     file=out,sep="",quote=FALSE,append=TRUE, row.names=FALSE, col.names=FALSE)
#==  End of 2nd loss and find procedure =======================#


#===================================================#
# Step 3: Final step 
#===================================================#

# Write the result. Appended mode
sink(out, append=TRUE)
summary(AR$RM)

if( length(ID_POOL) > 0 ){
   colnames(ID_POOL) <- "TF_Pair"
   (ID_POOL)
}
sink()

write.table(paste("#After Try: ", "R= ", round(AR$C,keta), " Y= ", nrow(AR$Y), " X= ", ncol(AR$X), sep=""),
            file=out,sep="",quote=FALSE,append=TRUE, row.names=FALSE, col.names=FALSE)

# R object (AR model)
save(AR, file=paste("fourth_result/result.txt", "_AR.RData", sep=""))

# x-fold Cross-validation
CV <- cv.park(AR$OX, AR$OY, FOLD, 10)
# R object (CV data)
save(CV, file=paste("fourth_result/result.txt", "_CV.RData", sep=""))

write.table(paste("#After Try CV: ", "LC= ", round(CV$LC,keta), " TC= ", round(CV$TC,keta), sep=""),
            file=out,sep="",quote=FALSE,append=TRUE, row.names=FALSE, col.names=FALSE)
#===================================================#

# Residual: fourth_result/result.txt.residual.txt is created
residual <- cbind(AR$Y, AR$PRED, (AR$Y - AR$PRED) )
colnames(residual) <- c("Observed", "Predicted", "Residual")

quantile(abs(residual[,1]))
quantile(abs(residual[,2]))
quantile(abs(residual[,3]))

temp <- cbind(rownames(AR$Y), round(residual,6))
colnames(temp) <- c("Gene","Observed", "Predicted", "Residual")
write.table( temp, file="fourth_result/result.txt.residual.txt", sep="\t", quote=FALSE, append=FALSE, row.names=FALSE,col.names=TRUE)



# pick up: +/- 1.65 90%, +/- 1.96 95%, +/- 2.58 99%
zthres <- c(1.65, 1.96, 2.58)
zid <- c("90per", "95per", "99per")
resiMean <- mean( residual[, 3] )
resiSD   <- sd( residual[, 3] )
zscore   <- (residual[,3]  - resiMean ) / resiSD
names(zscore) <- rownames(AR$Y)

#--- Plots ----------------------#
postscript("fourth_result/result.txt.ps")
op <- par(cex=1.5,pty="s")

maxx<-max(AR$PRED)
minx<-min(AR$PRED)
maxx<-max(maxx, AR$Y)
minx<-min(minx, AR$Y)

plot( AR$PRED, AR$Y, xlab="Pred.", ylab="Obs.",
     col="black", xlim=c(minx,maxx), ylim=c(minx,maxx), main="TEST")
#abline(h=0, v=0, col="black")
abline(0,1, col="red", lwd=2)
text(minx, maxx, paste0("R=",round(AR$C,3)), adj=c(0,1) )
#-------------------------#

plot_grid( nrow=1, ncol=2, 
	   coefplot(BR$RM, title="Before Reduced",intercept=FALSE, decreasing=FALSE, sort="magnitude"),
	   coefplot(AR$RM, title="After Reduced" ,intercept=FALSE, decreasing=FALSE, sort="magnitude")
    )
par(op)


op <- par(cex=1.2, pty="s", mfrow=c(1,3) )
for( i in 1:length(zthres)){
    cut <- zthres[i]
    idx <- abs(zscore) >= cut
    filename <- paste0("fourth_result/result.txt.outlier", "_" , zid[i] , ".txt")

    plot( AR$PRED, AR$Y, xlab="Pred.", ylab="Obs.",
	  col="black", xlim=c(minx,maxx), ylim=c(minx,maxx),
	  main=paste0("outliers, zscore >= ",cut, " (", sum(idx), " genes)") )
    
    if( sum(idx) ){
	    genes <- names(zscore)[idx]
	    #write.table( genes,file=filename, quote=FALSE, append=FALSE, row.names=FALSE,col.names=FALSE )
	    sink(filename, append=FALSE)
	    cat("# Zscore cut: >=", cut, "\n")
	    cat("# Gene Prediction Observation" , "\n")
	    sink()

	    t <- cbind(genes, round(AR$PRED[genes],6) )
	    t <- cbind(t, round(AR$Y[genes,],6) )
            write.table( t, file=filename, quote=FALSE, append=TRUE, row.names=FALSE,col.names=FALSE, sep=" ")

	    points( AR$PRED[genes], AR$Y[genes,], col="red", cex = 1.2, pch=16)
    }
    abline(0,1, col="red", lwd=2)
    text(minx, maxx, paste0("R=",round(AR$C,3)), adj=c(0,1) )
}
par(op)


op <- par(cex=1.5,pty="s")
qqnorm(residual[,1], main="Normal Q-Q plot (Observed)", col="black", pch=1, frame=FALSE)
qqline( residual[,1], col="red", lwd=2)

qqnorm(residual[,2], main="Normal Q-Q plot (Predicted)", col="black", pch=1, frame=FALSE)
qqline( residual[,2], col="red", lwd=2)

qqnorm(residual[,3], main="Normal Q-Q plot (Residual)", col="black", pch=1, frame=FALSE)
qqline( residual[,3], col="red", lwd=2)
#-------------------------- End of plots ------------


#===== End of GLMGE

dev.off()
q()
