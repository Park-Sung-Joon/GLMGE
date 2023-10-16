#!/usr/bin/env perl

###  GLMGE version 4.0 (perl script)

# STEP 1. Find $R_function from the same directory of this perl script
$R_function = "GLMGE.R";

undef(@data);
@data = split(/\//, $0);
pop(@data);
if(scalar(@data) > 0){
    $FUNC = join("/", @data) . "/" . $R_function;
}else{
    $FUNC = "./". $R_function;
}

if(-e $FUNC){
    $flag = 1;
}else{
    $flag = 0;
}
if(!$flag){
    print "\t*** " . $FUNC . " is required\n";
    exit;
}else{
    print "\t*** " . $FUNC . " was found\n";
}


#NOTE xRUN= run times for a X-Cross Validation. not used at here 
if(@ARGV != 7){
    print "%perl " . $0. " [INPUT MATRIX] [OUTPUT Filename] [ID] [Log10? 0|1] [RND_SEED] [SKIPgene_Filename|NULL] [SKIPcols: NULL | \"A,B,C...\"\n";
    print "\t" . "INPUT MATRIX formatted \"gene_name expression feature_1 feature_2 ...\" \n";
    print "\t" . "OUTPUT Filename created by this script" . "\n";
    print "\t" . "ID used for plot title, filename, etc " . "\n";
    print "\t" . "Log10: 0 or 1 converting the expression in log-10 scale" . "\n"; 
    print "\t" . "RND_SEED: a seed for random number generator" . "\n"; 
    print "\t" . "SKIPgene_Filename: list of gene names to be skip. or null" . "\n"; 
    print "\t" . "SKIPcols: list of feature names to be skip (comma seperation)." . "\n"; 

    if(!$flag){
	print "\t*** " . $FUNC . " is required\n";
    }else{
	print "\t*** " . $FUNC . " was found\n";
    }
    exit;
}

$TARGET_FILE = $ARGV[0];
if(&Chk_File($TARGET_FILE) == 0){ exit;}

$OUT_FILE = $ARGV[1];
$ID       = $ARGV[2]; $ID =~ s/\s+//g;
$LOG10    = $ARGV[3]; # 0 or 1
$SEED     = $ARGV[4]; # 12345
$SKIP_GENEFILE = $ARGV[5]; # filename: one gene name per one line
$SKIP_COLS     = $ARGV[6]; # filename: one feature name per one line
$XFOLD = 5; #5-fold Cross-validation

#--------------
# Set column names for skipping
$SKIP_COL_USE = 0;
$CNT_SKIPCOL = 0;
if($SKIP_COLS ne "NULL"){
    $SKIP_COL_USE = 1;

    undef(@data);
    @data = split(/\,/, $SKIP_COLS);

    $SKIP_COLS = ""; # make vector for R
    foreach $each (@data){
	$SKIP_COLS .= "\"" . $each . "\",";
	$CNT_SKIPCOL++;
    }
    $SKIP_COLS =~ s/\,$//;
}
#--------------

# Read feature names from ITEM line. The input file must have to line. See the example inputfile
# \#ITEM: feature_1 feature_2 feature3 ....
undef(@ITEM);
@ITEM = &GET_ITEM($TARGET_FILE);
print scalar(@ITEM) . " items loaded from " . $TARGET_FILE . "\n";

## comment out
##print join(" ", @ITEM) . "\n";

#================================#
# Genes to be skip
if($SKIP_GENEFILE ne "NULL"){
    if(!-e $SKIP_GENEFILE){
	print "Not found " . $SKIP_GENEFILE . "\n";
	print "Aborted.\n";
	exit;
	#$SKIP_GENEFILE = "NULL";
    }else{

	# set up gene names to be skip
	undef(@data);
	@data = `grep -v "\#" $SKIP_GENEFILE | cut -f 1`;
	chop(@data);
	if(scalar(@data) < 1){
	    $SKIP_GENEFILE = "NULL";
	}
	print "SKIP GENEs: " . scalar(@data) . "\n";
	undef(@data);
    }
}
#================================#

print "SKIP GeneFile: "  . $SKIP_GENEFILE . "\n";
print "SKIP_COL_USE: "   . $SKIP_COL_USE . "\n";
print "SKIP FEATRUREs: " . $CNT_SKIPCOL . "\n";
print "xFOLD: " . $XFOLD . "\n";


# Now.... Let's GO
# Files creation
# 1. OUT_FILE = trace of modeling
# 2. OUT_FILE.residual.txt = residuals (prediction and observation)
# 3. OUT_FILE.r.cmd = R script command lines
# 4. OUT_FILE.r.out = output from r.cmd 
# 5. OUT_FILE.ps = postscript file showing plots

unlink( $OUT_FILE );
&RUN_GLMGE($R_function, $TARGET_FILE, $OUT_FILE, \@ITEM, $LOG10, $XFOLD, $ID, $SEED, $SKIP_GENEFILE, $SKIP_COLS, $SKIP_COL_USE);

# these commands will be helpful
#%>tail -f OUT_FILE
#%>tail -f OUT_FILE.r.out
#%>watch -n 10 OUT_FILE

exit;

sub RUN_GLMGE{
    my($func, $input, $outfile, $ref_item, $log, $xfold, $id, $seed, $skipGeneFile, $skipCols, $skipCols_use) = @_;
    my $r_cmd   = $outfile . ".R.cmd";
    my $r_out   = $outfile . ".R.out";
    my $ps_file = $outfile . ".ps";
    my $resid_file = $outfile . ".residual.txt";

    my $min_features = 10;
    unlink($ps_file);
    unlink($outfile);
    unlink($resid_file);

    #ITEM name (=colname == feature names) for R
    if(@$ref_item < $min_features){
	print "Number of features less than " . $min_features . "\n";
	print "Aborted\n";
	exit;
    }

    $colname = "";
    for($i=0;$i<@$ref_item; $i++){
	$colname .= "\"" . $$ref_item[$i] . "\",";
    }
    $colname =~ s/\,$//;
    $skipCols =~ s/\s+//g;
    $seed =~ s/\s+//g;
    $skipGeneFile =~ s/\s+//g;
    $id =~ s/\s+//g;

#===================== From Here, generate R script of $r_cmd 
    open(R, ">$r_cmd");
    print R <<EOF;

library(MASS)
library(limma)
library(coefplot)
library(cowplot)

# load GLMGE functions
source("$func")

# set random seed
set.seed($seed, kind="Mersenne-Twister")

# Result file
out <- "$outfile"

# Genes and Features to be removed
GIVEN_SKIP_GENES <- NULL
SKIP_COLS <- NULL
if( $skipCols_use == 1){
    SKIP_COLS <- c($skipCols) # "AAA,CCC,BBB"
}

# Gene list given:  one gene per one line
if( "$skipGeneFile" != "NULL" ){
    tmp <- read.table("$skipGeneFile", header=FALSE, as.is=TRUE)
    #SKIP_GENES <- unique( c(SKIP_GENES, tmp[,1]) )
    # use only the 1st column	
    GIVEN_SKIP_GENES <- unique( tmp[,1] )
    remove(tmp)
}

#===================================================#
# Step 0: PREPARE INPUT DATA, not normalized
#         BR (before AIC reduction)
#===================================================#
coln  <- c($colname) #feature names

FOLD  <- 5           #X-fold CVs
ifelse($log == 1, logMode <- T, logMode <- F) #expression in log-10 scale


sink(out, append=FALSE)
cat("# INPUT FILE: " , "$input", "\\n")
cat("# OUT FILE: ", out, "\\n")
cat("# RND SEED: ", "$seed", "\\n")
cat("# GIVEN_SKIP_GENES: ", length(GIVEN_SKIP_GENES), "\\n")
cat("# GIVEN_SKIP_COLS:  ", length(SKIP_COLS), "\\n") 

cat("# Number of Features: ", length(coln), "\\n")
cat("# x-FOLD CVs: ", FOLD, "\\n")
cat("# log-10 scale: ", "$log", "\\n")
sink()

cnt_try <- 0
repeat{
  # "DB" keeps the original input data (X is randomly shuffled for this run, but Y-X pairs are retained)
  # NOTE: If colname includes "-" and other, anova automatically adds backquote `XXX` in some case.
  #       So, simpler name of features will be better.  
  DB <- RG_ReadDB("$input", coln, logMode, 0, GIVEN_SKIP_GENES, SKIP_COLS)

  length(coln)
  nrow(DB\$D); ncol(DB\$D)
  nrow(DB\$Y); ncol(DB\$Y)
  nrow(DB\$X); ncol(DB\$X)

  # Before AIC selection. This is the first look of this input file.
  BR      <- RG_Regression(DB\$Y, DB\$X) 
  cnt_try <- cnt_try + 1

  if(cnt_try > 10){ break } #Stop
  if(BR\$C > 0.1){ break }  #Stop

  # if the BR model looks bad,,, run "RG_ReadDB" again.
  # So, the function shuffles X again, which results differentially.
  # run <=10, if the low C is still .... this GLMGE will not be meaningful.
}
#===================================================#


#===================================================#
sink(out, append=TRUE)
cat("# Step 1: Reduce BR model", "\\n");
sink()
#         AR (after AIC reduction)
#===================================================#
cnt_try <- 0
repeat{

   # AIC selection, BR includes all TFs (excepting singular TFs)
   M <- RG_ReduceModel(BR)         #model 

   # Build X from the reduced model M
   X <- as.matrix(DB\$X[,M\$TERM]) #retained features in M and their orig data in DB
   colnames(X) <- M\$TERM
   rownames(X) <- rownames(DB\$Y)  #orig data

   # Y from the reduced model will be DB$Y
   # Now, Regression again with X and Y of the reduced model.
   # This is just for clerify
   AR <- RG_Regression(DB\$Y, X)
   # Then, run x-fold Cross-validationsl (OX and OY will be equal to X and Y).
   #CV <- cv.park(AR\$X, AR\$Y, FOLD, 1)
   CV <- cv.park(AR\$OX, AR\$OY, FOLD, 1)

   cnt_try <- cnt_try + 1
   if(cnt_try > 10){ break }
   if(AR\$C > 0.1){ break }
   
   # Try again due to something wrong happened
   # read the input file by shuffling X
   DB <- RG_ReadDB("$input", coln, logMode, 0, GIVEN_SKIP_GENES, SKIP_COLS)
   # and generate BR model       
   BR <- RG_Regression(DB\$Y, DB\$X) 
}


#===================================================#
sink(out, append=TRUE)
cat("# Step 2: Loss and Find. Back the removed feature to AR model", "\\n")
sink()
#===================================================#

# Features removed by AIC reduction
removed  <- M\$REMOVED
# number of the removed features. a stop condition
stop <- length(removed)

X  <- AR\$OX
Y  <- AR\$OY
rn <- rownames(AR\$OX)
cn <- colnames(AR\$OX)

# the 1st "Loss and Find" procedure
for( i in 1:stop){
   NX  <- as.matrix(DB\$X[,removed[i]]) #this removed feature vector for all the genes
   ID  <- removed[i]
   rownames(NX) <- rownames(DB\$X) # rowName = Gene IDs
   NX  <- NX[rn,]                  # pick up genes presented at AR model

   TX <- cbind(X, NX)              # extended model with this removed feature 
   colnames(TX) <- c(cn, ID)       # add this colname 
   rownames(TX) <- rn

   # Regression -> AIC
   TR <- RG_Regression(Y, TX)  # TX is the trial model with this removed feature
   TM <- RG_ReduceModel(TR)    # Reduce TR model again by AIC

   TX2 <- TX[, TM\$TERM]
   colnames(TX2) <- TM\$TERM
   rownames(TX2) <- rn
   TR  <- RG_Regression(Y, TX2)            # Run again for clarify
   TCV <- cv.park(TR\$OX, TR\$OY, FOLD, 1) # X-fold CVs.

   # Check the CV correlation of Trial Model (TR) and Update the current model (AR)
   if(TR\$C >= AR\$C){
      if(TCV\$TC > CV\$TC){
         AR <- TR
         CV <- TCV
      }
  }

} ## for i 1 to STOP

# Output. BR model. Before AIC reduction
sink(out, append=TRUE)
summary(BR\$RM)
sink()
write.table(paste("\#Before AIC: ", "R= ", round(BR\$C,2), " Y= ", nrow(BR\$Y), " X= ", ncol(BR\$X), sep=""), file=out,sep="",quote=FALSE,append=TRUE, row.names=FALSE, col.names=FALSE)

# Output. AR model. After AIC reduction. Appended mode
sink(out, append=TRUE)
summary(AR\$RM)
sink()

write.table(paste("\#After AIC: ", "R= ", round(AR\$C,4), " Y= ", nrow(AR\$Y), " X= ", ncol(AR\$X), sep=""),
            file=out,sep="",quote=FALSE,append=TRUE, row.names=FALSE, col.names=FALSE)

write.table(paste("\#After AIC CV: ", "LC= ", round(CV\$LC,4), " TC= ", round(CV\$TC,4), sep=""),
            file=out,sep="",quote=FALSE,append=TRUE, row.names=FALSE, col.names=FALSE)
#==  End of 1st loss and find procedure =======================#



# the 2nd "Loss and Find" procedure

# Find the removed features from AR model
variables <- colnames(AR\$RM\$model)[ 2:length(colnames(AR\$RM\$model)) ]
removed   <- c()
for(n in 1:length(DB\$cn)){

    f <- 1
    for(m in 1:length(variables) ){
        #already exist in the model?
        if(DB\$cn[n] == variables[m]){ f<-0; break }
    }

    if(f){ # not found from the current AR model

	#--------------------------------------------------
	# Comment out this, if the feature names have a pefix, e.g., "TFBS1", "TFBS2", "TFBS10".
	#--------------------------------------------------
	if( length ( grep("TFBS[0-9]+", DB\$cn[n]) ) > 0 ){
	    removed <- c(removed, DB\$cn[n])
	}

	# default. pick up the removed features
	#removed <- c( removed, DB\$cn[n] )
    }
}
#------------------#

# Start configuration

#========== These are options. Shuffle promoter TFBSs and suffle lrc TFBSs
# 1. select intraTFBS (TFBS found from the long-range contacted regions)
idx   <- grep("^intraTFBS", removed)
INTRA <- c()
if(length(idx)){
  # Random shuffling
  INTRA <- sample( removed[idx], replace=FALSE)
}

# 2. select not intraTFBS (= TFBS found from promoters)
idx  <- grep("^intraTFBS", removed, invert=TRUE)
TFBS <- c()
if(length(idx)){
  # Random shuffling
  TFBS <- sample( removed[idx], replace=FALSE)
}

# 3. Merge two randomly shuffled TFBSs arrays
removed <- c( TFBS, INTRA) #TFBS first
#================ end of the option

# Default. random shuffling of removed features
#removed <- sample( removed, replace=FALSE )

keta    <- 3 # number of digits after decimal point 
TRIED   <- c(); ID_POOL <- c()
k       <-  1;  LOOP    <- 1

## Find pairwise combinatorial features from the removed featurs
while(LOOP){

    # default stop condition
    LOOP <- 0
    stop <- length(removed)
    if(stop < 2){ break }

    #cat("X: ", colnames(DB\$X), "\\n" )
    #cat("removed: ", removed, "\\n")
    #cat("stop: ", stop, "\\n")
	
    # Here we go for finding a pairwise combination of features
    for( i in 1: (stop-1) ){
	# Information of the best model in this loop is saved to BESTxx
	BESTR   <- AR
        BESTTCV <- CV
        BESTID  <- BESTIDT <- c()
	
	#TFs already tried are never going to try again
	TRIED   <- unique( c(TRIED, removed[i]))

	# Set a model from the current AR model
	X  <- AR\$OX
	Y  <- AR\$OY
	rn <- rownames(AR\$OX)
	cn <- colnames(AR\$OX)

	updated <- 0
	for( j in (i+1):stop){
	    #cat("removed: ", removed, "\\n")
	    #cat("i x j ", i, " ", removed[i], " j" , j, " ", removed[j], "\\n")

	    #search pair only TFBSne:TFBS or Pax6:Pax1
            #this1 <- removed[i]; this2 <- removed[j]
	    #m1 <- ifelse( is.na(charmatch("TFBS", this1)), TRUE, FALSE)
	    #m2 <- ifelse( is.na(charmatch("TFBS", this2)), TRUE, FALSE)
	    #if( !xor(m1, m2) == FALSE ){ next; }
	    
	    #-------- Generate new input matrix ---------------#
	    # feature vector coming from input file
	    if( is.null( DB\$X[,removed[j]] ) ){ next }
	    if( is.na( DB\$X[,removed[j]] ) ){ next }
	    
	    NX <- as.matrix(DB\$X[,removed[i]] * DB\$X[,removed[j]])
		
	    # orignal pair name
            ID  <- paste( removed[i],":",removed[j],sep="")
	    # temporal pair name
	    IDT <- paste("PARKSJ", k , sep="")

	   # select only target genes
	   rownames(NX) <- rownames(DB\$X)
	   NX <- NX[rn,]

	   # merge columns with a new column (=new feature)
	   TX <- cbind(X, NX)
	   colnames(TX) <- c(cn, IDT)
	   rownames(TX) <- rn
	   #--------- done --------#

	   # Regression and AIC reduction
	   TR <- RG_Regression(Y, TX)
	   TM <- RG_ReduceModel(TR)

	   TX2 <- as.matrix( TX[,TM\$TERM] )
	   colnames(TX2) <- TM\$TERM
	   rownames(TX2) <- rn
	   TR  <- RG_Regression(Y, TX2)

	   ### Update because this is the best partner for TFi
	   temp1 <- round(TR\$C, keta)
	   temp2 <- round(BESTR\$C, keta)

	   flag <- 0
	   if( temp1 >= temp2){ flag = 1 } #improved
	    
	   # Check the improvement by x-Fold cross-validation 
	   if(flag){
	  	variables <- colnames(TR\$RM\$model)[2:length(colnames(TR\$RM\$model))]
		names(variables) <- variables

		if( !is.na(variables[IDT]) ){ #this model includes this valiable
		    # Do CV at here
		    TCV   <- cv.park(TR\$OX, TR\$OY, FOLD, 1)
	   	    temp3 <- round( TCV\$TC, keta )
		    temp4 <- round( BESTTCV\$TC, keta )
		    # This is very important. 
		    # If TC was improved, anyway accept it
		    # Even if TC was tie, C was imporoved, then accept it
		    if( (temp1 == temp2) && (temp3 <= temp4)){ next; } #nothing to do
		    
		    # Update the BEST model
		    updated <- 1
		    BESTR   <- TR
		    BESTTCV <- TCV
			
		    BESTID  <- ID
		    BESTIDT <- IDT
		    k = k + 1
         
		    # Write information to output file
		    write.table(paste("\#  Update local: ", "R= ", round(TR\$C,keta), " V= ", ID, " (", IDT, ")", 
				      " CV_LC= ", round(TCV\$LC,keta), " CV_TC= ", round(TCV\$TC,keta),  sep=""),
				file=out,sep="",quote=FALSE,append=TRUE,
				row.names=FALSE, col.names=FALSE
		    )

	        } # end of if( !is.na(variables[IDT]) )
	    } # end of flag
	    ###=================================#

    } # Next J

    # End of the searching J for the feature I

    # If, the model was updated, reset TF sets
    if(updated){

      if( length( grep(":", BESTID)) > 0 ){
         ID_POOL <- rbind(ID_POOL, BESTID)
         rownames(ID_POOL)[ nrow(ID_POOL)] <- BESTIDT
      }
      ID_POOL <- unique( ID_POOL)
 
      AR  <- BESTR
      CV  <- BESTTCV
      write.table(paste("\#Update: ", "R= ", round(AR\$C,keta), " Y= ", nrow(AR\$Y), " X= ", ncol(AR\$X),
                       " Added= ", BESTID, " (", BESTIDT, ")",
                       " CV_LC= ", round(CV\$LC,keta), " CV_TC= ", round(CV\$TC,keta), sep=""),
                 file=out,sep="",quote=FALSE,append=TRUE, row.names=FALSE, col.names=FALSE)

     #----------------------------------#
     # Get variables used in this finding
     variables <- colnames(AR\$RM\$model)[2:length(colnames(AR\$RM\$model))]
     names(variables) <- variables

     idx       <- grep("PARKSJ", variables)

     TID_POOL <- c()            # clear ID pool
     this_selected <- variables # default
    
     if(length(idx) > 0){ # if a new feature found
       this_selected <- variables[ -idx ] #remove "PARKKSJ" and retrieve each TF
    
       for(n in 1:length(idx)){
          M <- variables[ idx[n] ] # m must be in ID_POOL

          # Keep this information
          TID_POOL      <- rbind(TID_POOL, as.matrix(ID_POOL[M, ]))

          # Get each TF in this complex. Add 
          # this_selected <- c(this_selected, unlist( strsplit(ID_POOL[m,1], ":") ) )
          this_selected <- c(this_selected, unlist( strsplit(ID_POOL[M,1], ":") )[1])
          this_selected <- c(this_selected, unlist( strsplit(ID_POOL[M,1], ":") )[2])
        } #next N

        ID_POOL <- TID_POOL
     }
     names(this_selected) <- c()
     #------------------------------#

     # Update the featurs already tried
     this_selected <- unique( c(this_selected, TRIED))
   
     # Reset "removed" array. This is used in the next round.
     removed <- c()
     for(n in 1:length(DB\$cn)){
	 f <- 1 # not found in "this_selected" array
	 for(M in 1:length(this_selected)){
	     if(DB\$cn[n] == this_selected[M]){ f<-0; break }
	 }

	 if(f){
	     # This is an option
	     # update remained features. xxxxTFBSyyyy only
	     if( length ( grep("TFBS[0-9]+", DB\$cn[n]) ) > 0 ){
	     	 removed <- c(removed, DB\$cn[n])
	     }
	     
	     # Defalt mode
	     #removed <- c(removed, DB\$cn[n])
	 }
     }
     #------------------#


     # If features still left,
     # Random shuffling again!!!! So, next loop use different order
     if(length(removed) > 1){

	 # This is options;
	 idx <- grep("^intraTFBS", removed)
	 INTRA <- c()
	 if(length(idx)){
	     INTRA <- sample( removed[idx], replace=FALSE)
         }

	 idx <- grep("^intraTFBS", removed, invert=TRUE)
	 TFBS <- c()
	 if(length(idx)){
	     TFBS <- sample( removed[idx], replace=FALSE)
         }
	 removed <- c( TFBS, INTRA) #TFBS first
	 
	 # Default mode.
	 #removed <- sample(removed, replace=FALSE)
         # Loop again
	 LOOP <- 1

	 # skip other Ith. back to while(LOOP)
	 break #!!important
     }

    } #----------- End of Update -----------#
  
 } # Next Feature i
 

 # End ?
 if(LOOP == 0){ break }

}
#warnings()
#==  End of 2nd loss and find procedure =======================#


#===================================================#
# Step 3: Final step 
#===================================================#

# Write the result. Appended mode
sink(out, append=TRUE)
summary(AR\$RM)

if(length(ID_POOL) > 0){
   colnames(ID_POOL) <- "TF_Pair"
   (ID_POOL)
}
sink()

write.table(paste("\#After Try: ", "R= ", round(AR\$C,keta), " Y= ", nrow(AR\$Y), " X= ", ncol(AR\$X), sep=""),
            file=out,sep="",quote=FALSE,append=TRUE, row.names=FALSE, col.names=FALSE)

# R object (AR model)
save(AR, file=paste("$outfile", "_AR.RData", sep=""))

# x-fold Cross-validation
CV <- cv.park(AR\$OX, AR\$OY, FOLD, 10)
# R object (CV data)
save(CV, file=paste("$outfile", "_CV.RData", sep=""))

write.table(paste("\#After Try CV: ", "LC= ", round(CV\$LC,keta), " TC= ", round(CV\$TC,keta), sep=""),
            file=out,sep="",quote=FALSE,append=TRUE, row.names=FALSE, col.names=FALSE)
#===================================================#

# Residual: $resid_file is created
residual <- cbind(AR\$Y, AR\$PRED, (AR\$Y - AR\$PRED) )
colnames(residual) <- c("Observed", "Predicted", "Residual")

quantile(abs(residual[,1]))
quantile(abs(residual[,2]))
quantile(abs(residual[,3]))

residual <- cbind(rownames(AR\$Y), round(residual,6))
colnames(residual) <- c("Gene","Observed", "Predicted", "Residual")
write.table(residual,file="$resid_file", sep="\\t", quote=FALSE, append=FALSE, row.names=FALSE,col.names=TRUE)


#-------------------------#
postscript("$ps_file")
par(cex=1.5,pty="s")

maxx<-max(AR\$PRED)
minx<-min(AR\$PRED)
maxx<-max(maxx, AR\$Y)
minx<-min(minx, AR\$Y)

plot(AR\$PRED, AR\$Y, xlab="Pred.", ylab="Obs.",
     col="black", xlim=c(minx,maxx), ylim=c(minx,maxx), main="$id")
#abline(h=0, v=0, col="black")
abline(0,1, col="red", lwd=2)
text(minx, maxx, paste0("R=",round(AR\$C,3)), adj=c(0,1) )
#-------------------------#

plot_grid( nrow=1, ncol=2, 
	   coefplot(BR\$RM, title="Before Reduced",intercept=FALSE, decreasing=FALSE, sort="magnitude"),
	   coefplot(AR\$RM, title="After Reduced" ,intercept=FALSE, decreasing=FALSE, sort="magnitude")
    )

#===== End of GLMGE

dev.off()
q()
EOF
        
close(R);

# Run this R script
system("R CMD BATCH --no-save $r_cmd $r_out");
#$line = `tail -3 $r_out`;
#print $line;
#system("cat $r_out");

unlink($r_cmd);
unlink($r_out);

}


sub Chk_File{
    my $file = @_[0];
    
    if(!-e $file){
        print "## File: \" " . $file . "\" does not exist\n";
        return 0;
    }

    return 1;
}

sub GET_ITEM{
    my $file = @_[0];
    my ($line, @ITEM);

    # space seperated 
    # ITEM: A B C D

    $line = `grep \"ITEM\:\" $file`;
    chop($line);
    if(length($line) < 1){
	print "no \"# ITEM: \" record found from " . $file . "\n";
	exit;
    }

    $line =~ /^.*ITEM\:(.*)/;
    $line = $1;
    $line =~ s/^\s+//;
    $line =~ s/\s+$//;
    $line =~ s/\s+/\^/g; # merge multiple spaces to one

    # get feature name
    undef(@ITEM);
    @ITEM = split(/\^/, $line);

    if(scalar(@ITEM) < 1){
	print "failed to get feature names from \"ITEM: \" record\n";
	exit;
    }

    return @ITEM;
}
