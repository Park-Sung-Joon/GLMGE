#!/usr/bin/env perl

if( @ARGV != 3 ){
    print "%perl " . $0 . " [SCORE_MATRIX_FILE] [INPUT_ROOT_DIR] [OUT_DIR]\n";
    exit;
}

$SCORE_MATRIX = shift(@ARGV);
$INDIR        = shift(@ARGV);
$OUTDIR       = shift(@ARGV);

if( &Chk_File($SCORE_MATRIX) == 0 ){ exit; }
if( !defined($INDIR = &Chk_Dir($INDIR)) ){ exit; }

#-------- Check. The name of final result file should be "result.txt" 
$runResult_File = "result.txt";
print "Starting \"$runResult_File\" from " . $INDIR . "\n";

undef(@res_files);
@res_files = glob( $INDIR . "*/" . $runResult_File);
if(scalar(@res_files) < 1){
    print "No result files were found from " . $INDIR . "\n";
    print "Aborted\n";
    exit;
}
$howmany_runs = scalar(@res_files);
print $howmany_runs . " result files were found\n";
#------------------- End of Check


# Create Directories
$OUTDIR        = &myMKDIR( $OUTDIR );
$TFBS_OUTDIR   = &myMKDIR( $OUTDIR . "Features" );
$RESID_OUTDIR  = &myMKDIR( $OUTDIR . "Residuals");


#-----------------------------------
print "1. Get results from multiple runs\n";
#-----------------------------------

undef(%GLM); 
#$GLM{ run_id }{ "RSE" | "DF" }
#$GLM{ run_id }{ "X" | "Y" }
#$GLM{ run_id }{ "PARKSJ1" } = TFBS1:TFBS2
#$GLM{ run_id }{ "Feature" }{ "TFBS1" }   = Regression Coeff
#$GLM{ run_id }{ "CV_LC" | "CV_TC" }      = CV
&Read_Results( \@res_files, \%GLM);

$howmany_runs = 0;
$total_y = 0;
undef(@UNIQ_Features);
foreach $key (sort{$a cmp $b} keys %GLM ){
    $howmany_runs++;
    # pool all features found from all runs
    foreach $feature (sort{$a cmp $b} keys %{ $GLM{$key}{"Feature"} }){
	push( @UNIQ_Features, $feature );
    }
    $total_y = $total_y + $GLM{ $key }{ "Y"};
}
$total_y = sprintf("%.1f", $total_y / $howmany_runs);

# find unique features
@UNIQ_Features = &Get_Unique_Elements_Sort( @UNIQ_Features );
print "\t" . "**** " . $howmany_runs . " runs done\n";
print "\t" . "**** " . $total_y . " genes found\n";
print "\t" . "**** " . scalar(@UNIQ_Features) . " features found\n";
#======================

#-----------------------------------
print "2. Write the ensemble of feature RCs\n";
#-----------------------------------

# Feature name has been converted to {intra_}Transfac
$RESULT_features_RC = $TFBS_OUTDIR . "feature_RCs.txt";
unlink($RESULT_features_RC);
&Write_Table2( \@UNIQ_Features, \%GLM, $RESULT_features_RC );
if( -e  $RESULT_features_RC ){
    print "\t" . $RESULT_features_RC . " created\n";
}
#-----------------------------------


#-----------------------------------
print "3. Write Feature Stats\n";
#-----------------------------------

$T_THRES = 0.05; # cut for q-value
#feature_name, freq, meanofRC, sd,,, pv, qv......
# the feature_name is {intra_}Transfac

$RESULT_features_STAT = $TFBS_OUTDIR . "feature_stat.txt";
$ps_file = $RESULT_features_STAT . ".ps";
unlink($RESULT_features_STAT);
unlink($ps_file);

if( $howmany_runs > 1){ # more than 2 runs
    &Get_Feature_Stat( $RESULT_features_RC, $RESULT_features_STAT, $ps_file, $T_THRES);

    if( -e  $RESULT_features_STAT){
	print "\t" . $RESULT_features_STAT . " created\n";
    }

    if( !-e  $ps_file){
	print "\t" . "Failed " . $RESULT_features_STAT . "\n";
	print "Aborted\n";
	exit;
    }
    print "\t" . $ps_file. " created\n";
}else{
    print "\t" . "**** " . "\tSKIP creating " . $RESULT_features_STAT . " (runs < 2)\n";
}


#-----------------------------------
print "4. Write Residual Stats\n";
#-----------------------------------
$runResidual_File = $runResult_File . ".residual.txt";
undef(@resid_files);
@resid_files = glob( $INDIR . "*/" . $runResidual_File);
if(scalar(@resid_files) < 1){
    print "No residual files were found from " . $INDIR . "\n";
    print "Aborted\n";
    exit;
}

$RESULT_residuals_VEC  = $RESID_OUTDIR . "residual_vector.txt";
$RESULT_residuals_STAT = $RESID_OUTDIR . "residual_stat.txt";
$ps_file = $RESULT_residuals_STAT . ".ps";

unlink($RESULT_residuals_VEC);
unlink($RESULT_residuals_STAT);
unlink($ps_file);
undef(%OPR); #observed predicted residual
#$OPR{ gene_name }{"Observed"}  <- expression
#$OPR{ gene_name }{"Predicted"} <- vector of predicted
#$OPR{ gene_name }{"Residual"}  <- vector of obs - pred

&Get_RESIDUALs( \@resid_files, \%OPR, $RESULT_residuals_STAT, $RESULT_residuals_VEC, $ps_file, $howmany_runs); 
$howmany_genes = &CountHash(\%OPR);
print $howmany_genes . " genes\n";
if( -e  $RESULT_residuals_STAT){
    print "\t" . $RESULT_residuals_STAT . " created\n";
}


exit;




#====================== Sub-functions============    
sub Read_Results{
    my ($ref_resfile, $ref_hash) = @_;
    my ($i, $j, $run, $howmany, $ret, $this_runID, $file, @data);

    foreach $file ( @$ref_resfile ){

	# get run ID
	undef(@data);
	@data = split(/\//, $file);
	pop( @data );
	$this_runID = pop( @data );

	### Read final model (Variable, RC, Tv, Pv)
	undef(@V); undef(@C); undef(@T); undef(@P);
	my($RSE, $DF) = &Get_Regression_Coeff($file, \@V, \@C, \@T, \@P);

	$$ref_hash{ $this_runID}{"RSE"} = $RSE;
	$$ref_hash{ $this_runID}{"DF"}  = $DF;

	#-----------------------------#
	# second order effects
	undef(@data);
	@data = `grep \"^PARKSJ\" $file | grep \"\:\"`;
	chop(@data);

	if(scalar(@data) > 0){
	    foreach $each (@data){
		( $name, $this_pair ) = split(/\s+/, $each);
		$this_pair =~ s/\"//g;
		$$ref_hash{ $this_runID }{ $name } = $this_pair;
	    }
	}
	#-----------------------------#


	#-----------------------------#
	for($i=0; $i<scalar(@V); $i++){
	    if($V[$i] =~ /Intercept/){ next; }

	    # RC for this feature
	    if( $V[$i] =~ /^PARKSJ/){
		# get the name of pair
		my($pair1, $pair2) = split(/\:/, $$ref_hash{ $this_runID }{ $V[$i] });
		$$ref_hash{ $this_runID }{ "Feature" }{ $pair1 } = $C[$i];
		$$ref_hash{ $this_runID }{ "Feature" }{ $pair2 } = $C[$i];
	    }else{
		$$ref_hash{ $this_runID }{"Feature" }{ $V[$i] } = $C[$i];
	    }
	}
	#-----------------------------#

	$line = `grep \"After Try:\" $file`;
	chop($line);
	if(length($line) > 1){
	    undef(@data);
	    @data = split(/\s+/, $line);
	    $$ref_hash{ $this_runID}{"Y"} = $data[5]; # the final number of Y
	    $$ref_hash{ $this_runID}{"X"} = $data[7]; # the final number of X
	}

	$line = `grep \"After Try CV:\" $file`;
	chop($line);
	if(length($line) > 1){
	    undef(@data);
	    @data = split(/\s+/, $line);
	    $$ref_hash{ $this_runID }{ "CV_LC" } = $data[4];
	    $$ref_hash{ $this_runID }{ "CV_TC" } = $data[6];
	}

    } #next RUN

    return 1;
}



sub Get_Regression_Coeff{
    my ($file, $ref_v, $ref_c, $ref_t,$ref_p) = @_;
    my ($line, @data,$vari,@V,@C,@T,@P, $RSE, $DF, $line_number, $cnt);

    # get the final model
    undef(@data);
    @data = `grep -n Coefficients $file`; #multiple hits. use final one.  "line#:Coefficients"
    $found = $data[ scalar(@data)-1];
    undef(@data);
    @data = split(/\:/, $found);
    $line_number = $data[0];
    #------------------------

    $cnt=0;
    open(GRC, $file);
    while($line=<GRC>){
	$cnt++;
	if( $cnt < $line_number){ next; } #just skip

	if($line=~ /^Coefficients/){
	    $line=<GRC>; #colnum header
	    undef(@V); #variable
	    undef(@C); #coeff
	    undef(@T); #tvalue
	    undef(@P); # pv

	    while($line=<GRC>){
		undef(@data);
		chop($line);

		#========== This is the RG is rlm
		#rlm: one line format is "Value   Std. Error t value"
		# lm: #Estimate Std. Error t value Pr(>|t|)#

                #in mlr case or including "Signif. codes..."
		#if($line=~/^\-\-\-/){ last; }

                #--------------------- this is lm case
		# This is lm output case
	        # in the case of empty line 
		if($line=~/^\-\-\-/){ next; }
		if($line=~/^Signif/){ next; }
                if(length($line) < 1){ next; }

		if($line =~ /^Residual\s+standard\s+error\:\s(.*)\son\s(.*)\sdegrees.*/){
		    $RSE = $1;
		    $DF = $2;
		    last;
		}
		#RSE: Residual standard error

		#Name Value   Std. Error t value
		@data = split(/\s+/, $line);
                $vari = $data[0];     
                if($vari =~ /\`(.*)\`/){
	           $vari = $1;
           	}
		$pv = $data[4];
		$pv =~ s/\<//g;
		$pv =~ s/\>//g;
		$pv =~ s/\s+//g;

		push(@V, $vari);        # Name
		push(@C, $data[1]);   # Value
		push(@T, $data[3]);   # t value
		push(@P, $pv);           # empty or P-value
	    }
	}
    }
    close(GRC);
    #==================================================#

    @$ref_v = @V;
    @$ref_c = @C;
    @$ref_t = @T;
    @$ref_p = @P;

    return ($RSE, $DF);
}


#====================================#
#====================================#
#====================================#
#====================================#

sub Get_Feature_Stat{
    my($infile, $outfile, $ps_file, $t_thres )= @_;
    my $r_cmd   = $infile . ".r.cmd";
    my $r_out   = $infile . ".r.out";
    my ($line, @data, $i, $j);

    unlink($ps_file);

    $line = `grep -v "\#" $infile | head -1`;
    chop($line);

    undef(@data);
    @data= split(/\t/, $line);
    
    $R_COLNAMES = ""; # for R script
    for($i=0; $i<scalar(@data); $i++){
	# TC position is so important!!!
	if($data[$i] eq "TC"){
	    # start name of feature from here
	    for($j=$i+1; $j<scalar(@data); $j++){
		$R_COLNAMES .= "\"" . $data[$j] . "\",";
	    }
	}
    }
    $R_COLNAMES =~ s/\,$//;

    ## R script from here
    open(R, ">$r_cmd");
    print R <<EOF;

d     <- read.table("$infile", sep="\\t", header=TRUE, row.names=1)
start <- which( colnames(d)  == "TC") +1
end   <- ncol(d)

# read feature RCs
FEAT <- d[, c(start:end) ]
colnames(FEAT) <- c($R_COLNAMES)
rownames(FEAT) <- rownames(d)

# freq. of Columns
FREQ <- FEAT; FREQ [ FREQ != 0 ] <- 1
FREQ <- apply(FREQ, 2, sum) # vector

# Construct Stat Matrix
mean_w     <- apply( FEAT, 2, mean); names(mean_w) <- colnames(FEAT)
sd_w       <- apply( FEAT, 2, sd);   names(sd_w)   <- colnames(FEAT)

# two-tailed T-test
pv         <- apply( FEAT, 2, function(x){ t.test(x, mu=0)\$p.value } ); names(pv) <- colnames(FEAT)
# bonferroni adjusted pvalue
qv         <- p.adjust(pv, "bonferroni"); names(qv) <- colnames(FEAT)
signif     <- ifelse( qv < $t_thres, "Yes", "no"); names(qv) <- colnames(FEAT)

### RCM: freq man sd mean+sd mean-sd pv qv signif
RCM <- cbind( mean_w, sd_w)
RCM <- cbind( FREQ, RCM, mean_w + sd_w, mean_w - sd_w, pv, qv)
colnames(RCM)  <- c("freq", "mean", "sd", "mean+sd", "mean-sd", "pv", "qv")
rownames(RCM) <- names(mean_w)

# Write table to output file
tmp <- cbind(rownames(RCM), round(RCM,6), signif )
colnames(tmp) <- c("Feature", colnames(RCM), "signif(qv<$t_thres)" )
write.table(tmp, file="$outfile", sep="\\t", quote=FALSE, row.names=FALSE, col.names=TRUE, append=FALSE)


# Plotting
idx  <- signif == "Yes"
if( sum(idx) < 1){ q() } # stop
    
selFEAT <- RCM[idx,]
selFEAT.name <- rownames(RCM)[idx] # !important

y1  <- min( selFEAT[,"mean-sd"] )
y2  <- max( selFEAT[,"mean+sd"] )

postscript("$ps_file", horizontal=T)

# Plot RC and SD (no boxplot)
op <- par( cex=1.2 )
#this must be P-value test
idx <- sort( selFEAT[, "mean"], index.return=TRUE)\$ix
col <- ifelse( selFEAT[,"mean"] >= 0, "brown2", "dodgerblue" )

labY  <- c();
for(i in 1:length(idx)) {
    labY  <- c( labY, selFEAT[ idx[i], "mean+sd"] )
}

#idx is sorted index
plot( selFEAT[idx ,"mean"], type="n", pch=19, col=col[idx], cex=0.7,
#      ylab="RC",
      ylab="Regression Coeff.",
      xlab="Features (FDR<$t_thres)", 
      main  = "mean +/- sd plot",
      ylim=c(y1,y2), axes=F)

arrows(
    c(1:nrow(selFEAT)), selFEAT[ idx, "mean-sd"],
    c(1:nrow(selFEAT)), selFEAT[ idx, "mean+sd"],
    code=3, angle=90, length=0.03, col=col[idx]
)
points( selFEAT[idx ,"mean"], type="p", pch=19, col=col[idx], cex=0.7)

abline(h=0, lwd=1)
axis(2)

fsize <- 0.6
text( c(1:nrow(selFEAT)), labY, #+0.005,
      paste0("  ", rownames(selFEAT)[idx]),
      col=col[idx], srt=90, cex=fsize,
      adj=c(0,0.5),
      xpd=NA ) #c(0,1) )

# End
#============================


#============= Boxplot. using selFEAT.name (qv < thres)
## use FEAT not RCM
RCs <- FEAT[, selFEAT.name] # dim = runs x selFeatures
RCs [RCs == 0.0 ] <- NA #remove zero

# get quantile statics. not plot
bx <- boxplot( RCs, outline=F, na.rm=T, plot=FALSE)
BXstat <- bx\$stats
colnames(BXstat) <- bx\$names
# a row has 5 (box stats)
MV <- BXstat[3,]; names(MV) <- bx\$names

## sort using MV
idx  <- sort( MV, index.return=TRUE)\$ix
RCs2 <- RCs[, idx ] #for each box data
MV   <- MV[idx]     #sorted mean value

## coloring for default
posc <- "darksalmon"
negc <- "skyblue"
color_vector  <- rep(posc,  length(MV)) #for box
color_vector2 <- rep("red", length(MV)) # for borderline
color_vector3 <- rep("red", length(MV)) # for label

# counting and coloring using MV (mean value)
npos <- nneg <- 0
for(i in 1:length(MV)){
  if(MV[i] <= 0) { color_vector[i] <- negc; nneg <- nneg+1; color_vector2[i]<-"blue"; color_vector3[i] <- "blue" }
  if(MV[i] > 0) { npos <- npos+1 }
}

# Labelling
labY  <- c(); aj <- c();
for(i in 1:ncol(RCs2)){
    name <- colnames(RCs2)[i]
    labY  <- c(labY, BXstat[5, name] ) #upper 
  if(MV[i] < 0){
     aj <- c(aj, 0)
  }else{
     aj <- c(aj, 1)
  }
}



fsize <- ifelse(ncol(RCs2)<20, 0.7, 0.6)
bx <- boxplot( RCs2, outline=F, na.rm=T,
               col=color_vector,
               names=NA, axes=F,
               #main=paste("$STAGE (pos=", npos, ", neg=",nneg, ")",sep=""),
               main="$THIS_ID",
               sub=paste0("(Activator=", npos, ", Repressor=",nneg, ", FDR<$t_thres)"),
#              sub="($maximal_genes genes)",
              ylab="Regression Coeff.",
#               ylab="Regulatory Activity",
               border=color_vector2,
               lwd=0.8, cex=fsize)

abline(h=0, lwd=1)
axis(2)

text( c(1:ncol(RCs2)), labY, #+0.005,
      paste0("  ", colnames(RCs2), sep=""),
      col=color_vector3, srt=90, cex=fsize,
      adj=c(0,0.5), xpd=NA) #c(0,1) )


par(op)

dev.off()
q()
EOF
        
close(R);

    # Run it
    system("R CMD BATCH --no-save $r_cmd $r_out");
#    undef(@data);
#    @data = `tail -3 $r_out`;
    #system("cat $r_out");
    #if( !($data[0] =~ /^\>\sproc\.time.*/) ){
#	system("cat $r_out");
#	exit;
#    }
#    unlink($r_cmd);
#    unlink($r_out);
}



sub Write_Table2{
    my( $ref_Features, $ref_GLM, $outfile ) = @_;
    my ($i, $key, $C);

    open(OUT, ">$outfile");
    print OUT "# Feature RCs found in each run\n";
    print OUT "Run\t" . "Gene\t" . "Features\t" . "RSE\t". "DF\t";
    print OUT "LC\t" . "TC";
    for($i=0; $i < scalar( @{$ref_Features} ); $i++){
	print OUT "\t" . $$ref_Features[$i];
    }
    print OUT "\n";

    foreach $key (sort{$a cmp $b} keys %{$ref_GLM} ){
	$RSE   = $GLM{ $key }{"RSE"};
	$DF    = $GLM{ $key }{"DF"};
	$X     = $GLM{ $key }{"X"};
	$Y     = $GLM{ $key }{"Y"};
	$CV_LC = $GLM{$key}{"CV_LC"};
	$CV_TC = $GLM{$key}{"CV_TC"};

	print OUT $key   . "\t"; # run ID
	print OUT $Y     . "\t" . $X . "\t" . $RSE . "\t" . $DF . "\t";;
	print OUT $CV_LC . "\t" . $CV_TC;
	for($i=0; $i<scalar(@{$ref_Features}); $i++){
	    $C = $$ref_GLM{ $key }{"Feature"}{ $$ref_Features[$i] };
	    if($C eq ""){ $C = 0.0; } # not found in this regression model

	    print OUT "\t" . $C;
	}
	print OUT "\n";
    }
    close(OUT);

    return 1;
}


sub Get_RESIDUALs{
    my($ref_files, $ref_OPR, $output_file, $vector_file, $ps_file, $howmany_runs ) = @_;
    my($file, $line, @data, $cnt_genes, $gene);
       
    foreach $file ( @{$ref_files} ){
	print $file . "\n";

	# build vectors
	open(GR, $file);
	while($line = <GR>){
	    if($line =~ /Observed/){ next; }

	    chop($line);
	    undef(@data);
	    @data = split(/\t/, $line);
	    
	    $$ref_OPR{ $data[0] }{"Observed"}   = $data[1];
	    $$ref_OPR{ $data[0] }{"Predicted"} .= $data[2] . " "; #vectors
	    $$ref_OPR{ $data[0] }{"Residual"}  .= abs( $data[3] ) . " "; #abs vectors
	}
	close(GR);
    }

    # check the length of vectors
    foreach $gene (sort {$a cmp $b} keys %{$ref_OPR} ){

	#-------------------------------------------
	$line = $$ref_OPR{ $gene }{"Predicted"};
	$line =~ s/\s+$//;

	undef(@data);
	@data = split(/\s+/, $line);

	if( $howmany_runs != scalar(@data) ){
	    print $gene . "\t" . "the precited vectors was not " . $howmany_runs . " in length\n";
	    print "Aborted\n";
	    exit;
	}
	#-------------------------------------------

	#-------------------------------------------
	$line = $$ref_OPR{ $gene }{"Residual"};
	$line =~ s/\s+$//;

	undef(@data);
	@data = split(/\s+/, $line);
	
	if( $howmany_runs != scalar(@data) ){
	    print $gene . "\t" . "the residual vectors was not " . $howmany_runs . " in length\n";
	    print "Aborted\n";
	    exit;
	}
    }

    # Step 3. create output file
    open(GR, ">$vector_file");
    foreach $gene (sort {$a cmp $b} keys %{$ref_OPR} ){
	#-------------------------------------------
	$line = $$ref_OPR{ $gene }{"Predicted"};
	$line =~ s/\s+$//;
	$$ref_OPR{ $gene }{"Predicted"} = $line;
	#-------------------------------------------

	#-------------------------------------------
	$line = $$ref_OPR{ $gene }{"Residual"};
	$line =~ s/\s+$//;
	$$ref_OPR{ $gene }{"Residual"} = $line;
	#-------------------------------------------

	# gene observed runs PredVector ResidualVector
	print GR $gene . "\t" . $$ref_OPR{$gene}{"Observed"} . "\t" . $howmany_runs . "\t" . $$ref_OPR{$gene}{"Predicted"} . "\t" . $$ref_OPR{$gene}{"Residual"} . "\n";
    }
    close(GR);

    # tfbs_found is used at here
    &R_ResidualZscore($vector_file, $output_file, $ps_file);

    return 1;
}

# Use R 3.6
sub R_ResidualZscore{
    my ($input, $output, $ps_file)  = @_;
    my $r_cmd   = $input . ".r.cmd";
    my $r_out   = $input . ".r.out";

    unlink($ps_file);
    unlink($output);

    $temp = $output . ".outlier_*.txt";
    system("rm -fr $temp");

    open(R, ">$r_cmd");
    print R <<EOF;

d   <- read.table("$input", sep="\\t", header=FALSE, row.names=1)

pred <- c()
for(i in 1:nrow(d)){
    t <- as.numeric( unlist( strsplit(as.character(d[i,3]), " "))) 
    pred <- rbind( pred, t)
}
rownames(pred) <- rownames(d)

resid <- c()
for(i in 1:nrow(d)){
    resid <- rbind( resid, as.numeric( unlist( strsplit( as.character(d[i,4]), " "))) )
}
resid <- abs(resid)
rownames(resid) <- rownames(d)

#---------------
predMean <- apply(pred, 1, mean)
predSD   <- apply(pred, 1, sd)
resiMean <- apply(resid, 1, mean)
resiSD   <- apply(resid, 1, sd)
#---------------


#---------------
mean_resiMean <- mean( resiMean )
sd_resiMean   <- sd( resiMean )
#---------------

zscore <- (resiMean - mean_resiMean ) / sd_resiMean

TBL <- cbind(d[,1], d[,2], predMean, predSD, resiMean, resiSD, zscore)
TBL <- round(TBL, 6)
TBL <- cbind(rownames(d), TBL)

colnames(TBL) <- c( "Gene", "Observed", "runs", "mean_Predicted", "sd_Predicted", "mean_absResidual", "sd_absResidual", "zscore" )
write.table(TBL,file="$output", sep="\\t", quote=FALSE, append=FALSE, row.names=FALSE, col.names=TRUE)

# pick up
# +/- 1.65 90%, +/- 1.96 95%, +/- 2.58 99%
zthres <- c(1.65, 1.96, 2.58)
zid <- c("90per", "95per", "99per")

d <- read.table("$output", sep="\\t", header=TRUE, row.names=1)
minx<-min( d[,c(1,3)] )
maxx<-max( d[,c(1,3)] )

postscript("$ps_file", horizontal=T)

op <- par(pty="s", cex=1.2)
DF   <- data.frame( d[,"mean_Predicted"], d[,"Observed"] )
minx <- min(DF)
max  <- max(DF)
corr <- cor(DF[,1], DF[,2])

library(LSD)
heatscatter(DF[,1],DF[,2], cex=1.2, pch=19,
	    xlab="Predicted (mean)", ylab="Observed",
	    xlim=c(minx,maxx), ylim=c(minx,maxx), main="" )
legend("topleft", paste("R=", round( corr, 2)), bty="n", cex= 1.1)
abline(0,1,col="red")

### add corr to output file
cat( paste0("#R=", round(corr,4)),"\n" , file="$output", append=TRUE)
par(op)


op <- par(mfrow=c(2,3), cex=1.1, pty="s")

hist( d[,"Observed"], xlab="Expression", ylab="Frequency", col="gray", main="")
hist( d[,"mean_absResidual"], xlab="Residual", ylab="Frequency", col="gray",main="" )

hist( d[,"zscore"], xlab="Zscore", ylab="Frequency", col="gray",main="" )
abline(v=zthres, col="red", lwd=1)

for( i in 1:length(zthres)){
    cut <- zthres[i]
    idx <- abs( d[,"zscore"] ) >= cut

    filename <- paste0("$output", ".outlier_" , zid[i] , ".txt")

    C <- as.vector( cor(d[,"mean_Predicted"], d[,"Observed"] ) )
    plot( d[,"mean_Predicted"], d[,"Observed"], col="gray", xlab="Predicted", ylab="Observed",
	  xlim=c(minx,maxx), ylim=c(minx,maxx),
	  main=paste0("zscore >= ", cut),
	  sub=paste0("R=", round(C,3), " ", sum(idx), " genes")
	  )
    
    if(sum(idx)){
	genes <- rownames(d)[idx]

        sink(filename, append=FALSE)
        cat("# Zscore cut: >=", cut, "\\n")
        cat("# Gene Prediction Observation" , "\\n")
        sink()

        write.table(genes, file=filename, quote=FALSE, append=FALSE, row.names=FALSE, col.names=FALSE, sep=" ")

	points( d[genes, "mean_Predicted"], d[genes, "Observed"], col="red", pch=16, cex=1.2)
    }
    abline(0,1,col="red", lwd=2)
}

dev.off()
q()
EOF
        
close(R);

    # Run it
    system("R CMD BATCH --no-save $r_cmd $r_out");
    #undef(@data);
    #@data = `tail -3 $r_out`;
    #system("cat $r_out");
    #if( !($data[0] =~ /^\>\sproc\.time.*/) ){
#	system("cat $r_out");
#	exit;
    #}
    #unlink($r_cmd);
    #unlink($r_out);
}




# $howmany = &CountHash(\%hash)
sub CountHash{
    my $ref_hash = $_[0];
    return keys %$ref_hash;
}


# $dir = &myMKDIR($dir);
sub myMKDIR{
    my $dir = $_[0];

    if(substr($dir, -1) ne "/"){
	$dir .= "/";
    }

    if( !-e $dir){
	print "Created " . $dir . "\n";
	system("mkdir -p $dir");
    }
    return $dir;
}
	

#if(&Chk_File($file) == 0){ exit; }
sub Chk_File{
    my $file = $_[0];
    
    if(!-e $file){
        print "## File: \" " . $file . "\" does not exist\n";
        return 0;
    }

    return 1;
}


#if( !defined($dir = &Chk_Dir($dir)) ){ exit; }
sub Chk_Dir{
    my $dir = $_[0];

    if(!-e $dir){
        print "## Directory: \" " . $dir . "\" does not exist\n";
        return undef($dir);
    }

    if(substr($dir, -1) ne "/"){
	$dir .= "/";
    }

    return $dir;
}

## get unique and sort character
sub Get_Unique_Elements_Sort{
    my @data = @_;
    my (%tmp, @unique, @sorted);

    @unique = grep(!$tmp{$_}++, @data);
    @sorted = sort{$a cmp $b} @unique;

    return @sorted;
}
