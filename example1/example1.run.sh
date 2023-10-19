
if [ -e first_result ]; then
	rm -fr first_result
fi

if [ -e second_result ]; then
        rm -fr second_result
fi

if [ -e third_result ]; then
	rm -fr third_result
fi

if [ -e fourth_result ]; then
	rm -fr fourth_result
fi

### "../src/GLMGE.R" is required.
PROG="../src/GLMGE_v4.pl"
INPUT="../src/matrix/example.matrix.txt"

#------------ Run regression modeling
echo "first_result"
perl ${PROG} ${INPUT} first_result/result.txt TEST 1 123456 NULL NULL

#------------ Run again the regression modeling
#------------ after removing outliers (1% in zscore)
skip_geneList=first_result/result.txt.outlier_99per.txt
if [ -e ${skip_geneList} ]; then
        cat first_result/result.txt.outlier_99per.txt > "skipped_genes_in_2nd.txt"

	echo "second_result"
	perl ${PROG} ${INPUT} second_result/result.txt TEST 1 123456 skipped_genes_in_2nd.txt NULL
fi

#------------ Do you want to run once again?
skip_geneList=second_result/result.txt.outlier_99per.txt
if [ -e ${skip_geneList} ]; then
        #merge the first and second outliers
        cat first_result/result.txt.outlier_99per.txt second_result/result.txt.outlier_99per.txt > "skipped_genes_in_3rd.txt"

	echo "third_result"
        perl ${PROG} ${INPUT} third_result/result.txt TEST 1 123456 skipped_genes_in_3rd.txt NULL
fi

#------------ One more?
skip_geneList=third_result/result.txt.outlier_99per.txt
if [ -e ${skip_geneList} ]; then
        #merge the first and second and third outliers
        cat first_result/result.txt.outlier_99per.txt second_result/result.txt.outlier_99per.txt third_result/result.txt.outlier_99per.txt > "skipped_genes_in_4th.txt"

	echo "fourth_result"
        perl ${PROG} ${INPUT} fourth_result/result.txt TEST 1 123456 skipped_genes_in_4th.txt NULL
fi

echo "DONE"
