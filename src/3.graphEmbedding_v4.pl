#!/usr/bin/env perl

if( @ARGV != 2 ){
    print "%perl " . $0 . " [SCORE_MATRIX_FILE] [INPUT_ROOT_DIR]\n";
    exit;
}

$SCORE_MATRIX = shift(@ARGV);
$INDIR  = shift(@ARGV);

if( &Chk_File($SCORE_MATRIX) == 0 ){ exit; }
if( !defined($INDIR = &Chk_Dir($INDIR)) ){ exit; }

# Where the line programs are located??
# Same directory with this script
undef(@data);
@data = split(/\//, $0);
pop(@data);
$LINE_PROG_DIR = join("/", @data) . "/LINE/linux/";
$LINE_RECONS = $LINE_PROG_DIR . "reconstruct";
$LINE_PROG   = $LINE_PROG_DIR . "line";
$LINE_NORM   = $LINE_PROG_DIR . "normalize";
$LINE_CONC   = $LINE_PROG_DIR . "concatenate";
if( &Chk_File($LINE_RECONS) == 0 ){ exit; }
if( &Chk_File($LINE_PROG) == 0 ){ exit; }
if( &Chk_File($LINE_NORM) == 0 ){ exit; }
if( &Chk_File($LINE_CONC) == 0 ){ exit; }


# Create Directories
$NW_OUTDIR      = &myMKDIR( $INDIR . "Network");
$RESULT_nw_sif  = $NW_OUTDIR . "cytoscape.sif";
$RESULT_nw_edge = $NW_OUTDIR . "cytoscape.edge.txt";
$RESULT_nw_node = $NW_OUTDIR . "cytoscape.node.txt";

# network embedding: source space target degree(=weight)
$LINE_input         = $NW_OUTDIR . "LINE_input.txt";
$RESULT_line_log    = $NW_OUTDIR . "LINE_log.txt";
$RESULT_line_vector = $NW_OUTDIR . "LINE_vectors.txt";
$RESULT_line_tsne   = $NW_OUTDIR . "LINE_tsne.ps";


$RESULT_features_STAT = $INDIR . "Features/feature_stat.txt";
if( &Chk_File($RESULT_features_STAT) == 0 ){ exit; }

#-----------------------------------
print "1. Get score matrix\n";
#-----------------------------------

undef( %EXPRESSION );
undef( @ITEMS );
undef( %SCORES );
&Get_Score_Info( $SCORE_MATRIX, \%EXPRESSION, \@ITEMS, \%SCORES);
$howmany_genes = &CountHash(\%EXPRESSION);
$howmany_items = scalar( @ITEMS );
print "\t*** " . "Genes: "    . $howmany_genes . "\n";
print "\t*** " . "Features: " . $howmany_items . "\n";


#-----------------------------------
print "2. Create Cytoscape Network\n";
#-----------------------------------
## Using only significant feature only
&Create_NW( $RESULT_features_STAT, \%EXPRESSION, \%SCORES, $RESULT_nw_sif, $RESULT_nw_edge, $RESULT_nw_node, $LINE_input);


#-----------------------------------
print "3. Run Graph Embedding\n";
#-----------------------------------
if(!-e $LINE_input){
    print "LINE input: ". $LINE_input . " was not found\n";
    print "Aborted\n";
    exit;
}

&Run_LINE( $LINE_input, $RESULT_line_log, $RESULT_line_vector);

print "\tRunning TSNE........\n";
&Plot_LINE($RESULT_line_vector, $RESULT_nw_node, $RESULT_line_tsne);
print $RESULT_line_tsne . "\n";

exit;


sub Run_LINE{
    my($input, $log, $result) = @_;
    my( $vec1, $vec2, $vec3, $vec4, $command0, $command1, $command2);
    my($command3, $command4, $command5,$command6);

    $vec1 = $result . ".1";
    $vec2 = $result . ".2";
    $vec3 = $result . ".1.norm";
    $vec4 = $result . ".2.norm";

    #LINE Programs defined in the main 
    $command0 = $LINE_RECONS . " -train " . $input . " -output " . $input . ".reconst -depth 2 -threshold 1";
    $command1 = $LINE_PROG   . " -train " . $input . ".reconst " . " -output " . $vec1 . " -binary 1 -size 100 -order 1 -negative 5";
    $command2 = $LINE_PROG   . " -train " . $input . ".reconst " . " -output " . $vec2 . " -binary 1 -size 100 -order 2 -negative 5";
    $command3 = $LINE_NORM   . " -input " . $vec1  . " -output " . $vec3 .  " -binary 1";
    $command4 = $LINE_NORM   . " -input " . $vec2  . " -output " . $vec4 .  " -binary 1";
    $command5 = $LINE_CONC   . " -input1 " . $vec3 . " -input2 " . $vec4 . " -output " . $result . " -binary 0";
    $command6 = "rm -f " . $vec1 . " " . $vec2 . " " . $vec3 . " " . $vec4;

    unlink($log);
    open(OUT, ">$log");
    print OUT $command0 . "\n";
    print OUT $command1 . "\n";
    print OUT $command2 . "\n";
    print OUT $command3 . "\n";
    print OUT $command4 . "\n";
    print OUT $command5 . "\n";
    print OUT $command6 . "\n";
    print OUT "\n";
    close(OUT);

    # Run it
    print "\tRunning LINE ..........." . "\n";
    system($command0 . " >> " . $log);
    system($command1 . " >> " . $log);
    system($command2 . " >> " . $log);
    system($command3 . " >> " . $log);
    system($command4 . " >> " . $log);
    system($command5 . " >> " . $log);
    system($command6);
    print "\tFinished\n";
    print "\t" . $log . "\n";

    return 1;
}





sub Create_NW{
    my($TFBS_file, $ref_Expression, $ref_SCORE, $sif_file, $edge_file, $node_file, $line_input) = @_;
    my($i, $tmp1, $tmp2, @data, $FID);

    # Capture "Yes" feature only 
    # Name Frequncy meanRC
    undef(@data);
    @data = `awk '{if(\$9 == \"Yes\"){ print \$1 \"\\t\" \$2 \"\\t\" \$3} }' $TFBS_file`;
    chop( @data );
    
    #----
    undef(%SIGNIF_RC);
    undef(%FEATURE_FREQ);
    foreach $each (@data){
	# feature_name is ITEM, not transfc
	my($FID, $freq, $rc_mean) = split(/\t/, $each);
	$SIGNIF_RC{ $FID }   = $rc_mean;
	$FEATURE_FREQ{ $FID} = $freq;
    }
    $howmany_signif = &CountHash( \%SIGNIF_RC );
    print "\t*** " . $howmany_signif . " Features found\n";

    if( $howmany_signif < 1){
	print "Aborted: zero sifnificant Features\n";
	exit;
    }
    #----

    undef(@NODES);
    undef(@EDGES);
    undef(@SIF);
    foreach $gene (sort {$a cmp $b} keys %{$ref_SCORE})  {
	$expression = $$ref_Expression{ $gene };

	foreach $feature (sort {$a cmp $b} keys %{$$ref_SCORE{$gene} } ){
	    if( $$ref_SCORE{$gene}{$feature} == 0.0 ) { next; }

	    if( $SIGNIF_RC{ $feature } ne ""){
		# cytoscape network file
		push(@SIF, $gene . " gp " . $feature . "\n");
	
		#edge type expression freq RC TFBSname
		$line  = $gene . " (gp) "   . $feature . "\t" . "Gene-TF" . "\t";
		$line .= $expression . "\t" . $FEATURE_FREQ{$feature} . "\t";
		$line .= $SIGNIF_RC{ $feature } . "\t" . $feature . "\n";
		push(@EDGES, $line);

		#node_name type value
		push(@NODES, $gene    . "\t" . "Gene"    . "\t" . $expression . "\n");
		push(@NODES, $feature . "\t" . "Feature" . "\t" . $SIGNIF_RC{ $feature } . "\n");
	    }

	}
    } # Next gene

    #remove redundancy
    @EDGES = &Get_Unique_Elements_Sort( @EDGES );
    @NODES = &Get_Unique_Elements_Sort( @NODES );
    @SIF   = &Get_Unique_Elements_Sort( @SIF );

    # Create files
    open(OUTSIF, ">$sif_file"  );
  
    open(OUTND,  ">$node_file" );
    print OUTND "Node\tType\tValue(Expres., RC)\n";

    open(OUTED,  ">$edge_file" );
    print OUTED "Edge\tType\tExpression\tFreq\tRC\tTFBS\n";

    print OUTSIF @SIF;
    print OUTND @NODES;
    print OUTED @EDGES;

    close(OUTND);
    close(OUTED);
    close(OUTSIF);

    open(LOUT, ">$line_input");
    open(INED, "$edge_file"  );
    while($line = <INED>){
	if( $line =~ /^Edge/ ){ next; }
	chop($line);

	undef(@data);
	@data = split(/\t/, $line);
	my($node1, $temp, $node2) = split(/\s+/, $data[0]);
	$freq = $data[3];

	# bi-directional
	print LOUT $node1 . " " . $node2 . " " . $freq . "\n";
	print LOUT $node2 . " " . $node1 . " " . $freq . "\n";
    }
    close(INED);
    close(LOUT);

    return 1;
}


sub Plot_LINE{
    my ($input, $info_file, $ps_file) = @_;

    my $r_cmd  = $input . ".R.cmd";
    my $r_out  = $input . ".R.out";
    unlink($ps_file);

    ## R/3.6.0
    open(R, ">$r_cmd");
    print R <<EOF;

library(ggplot2)
library(Rtsne)
#library(ggrepel)

Info <- read.table("$info_file", header=T, row.name=1, sep="\\t")
VEC  <- read.table("$input",     skip=1, row.names=1)
Type <- Info[ rownames(VEC),  "Type"]  #DEG, DEG|TF, others

set.seed(123)

tmp <- Rtsne::Rtsne(VEC, check_duplicates=FALSE, pca=TRUE, perplexity=30, theta=0.5, dims=2)

tsED <- data.frame(
    V1=tmp\$Y[,1],  V2=tmp\$Y[,2],
    color=Type,
    label=rownames(VEC),
    category=Type
)

#grouped <- kmeans(x=tsED, centers=3, nstart=100)
my.theme.squre <- theme (
    plot.title=element_text(hjust=0.5, size=16) 
    ,aspect.ratio=1
    ,axis.text.x=element_text(size=14)
    ,axis.text.\y=element_text(size=14)
    #,legend.position=c(1,0)
    #,legend.justification=c(0,0)
    ,legend.box.margin=margin(c(0,0,0,2)) 
    
    ,axis.title =element_blank(),
    legend.title=element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.text=element_text(size=16)
 )

col <- c(
    "Gene"     ="coral1",
    "Feature"  ="darkgreen"
)

shape <- c(
    "Gene"     = 19,
    "Feature"  = 24
)

p <- ggplot(tsED,  aes(x=V1, \y=V2, color=color, shape=color), size=4)
p <- p + geom_point( aes(color=color, shape=color), size=3 )
p <- p + scale_colour_manual( values = col )
p <- p + scale_shape_manual( values = shape)
#p <- p + geom_text_repel(aes(label=label), size=2, segment.size=0.1, max.overlaps=20)
p <- p + theme_bw() + my.theme.squre
p <- p + guides(color=guide_legend(override.aes=list(size=5)))

postscript( "$ps_file", horizontal=TRUE)

op <- par(cex=1.2)
p
par(op)
dev.off()
q()
EOF
close(R);
	    
    system("R CMD BATCH --no-save $r_cmd $r_out");
    #system("cat $r_out");
    
    undef(@data);
    @data = `tail -3 $r_out`;
    if( !($data[0] =~ /^\>\sproc\.time.*/) ){
	print @data;
    }
#    unlink($r_cmd);
#    unlink($r_out);
	    
    return 1;
}


sub Get_Score_Info{
    my ($matrix, $refTPMs, $refITEMS, $refSCORES) = @_;
    my ($i, $line, $gene, $express, $vector, @data, $ITEM);

    open(GSI, $matrix);
    while($line = <GSI>){
	chop($line);

	# Get scores: {gene}{tfbs} => scores
	if( !($line =~ /^\#/) ){
	    undef(@data);
	    @data = split(/\t/, $line);

	    $gene     = shift(@data);
	    $express  = shift(@data);
	    $vector   = shift(@data); #<--- vector !!!

	    $$refTPMs{ $gene } = $express;
	    $VEC{ $gene }      = $vector; ## anyway keep this vector
	}
    }
    close(GSI);

    # ITEMs == feature names
    $line = `grep \"ITEM:" $matrix`;
    chop($line);
    $line =~ /.*ITEM:\s(.*)/;
    $line = $1;
    @$refITEMS = split(/\s+/, $line);

    foreach $gene (keys %VEC){
	undef(@data);
	@data = split(/\s+/, $VEC{ $gene });

	if(scalar(@data) != scalar(@$refITEMS)){
	    print "\t" . "the lengths of ITEMs and Score vector are different." . "\n";
	    print "\t" . $matrix . "\n";
	    print "\t" . "Aborted\n";
	    exit;
	}


	# $$refSCORES{ $gene }{ $$refITEMS[$i] }
	for($i=0; $i<scalar(@data); $i++){
	    $$refSCORES{ $gene }{ $$refITEMS[ $i] } = $data[$i];
	}
    }

    return 1;
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
