# in 05.Beta_pvalue_check directory
for i in $(ls chr1-22.250227.*0.1M.rankPbeta+.hybrid);do out=$(basename $i | cut -d "." -f 1,2,3,4);cat $i | awk -F " " '{if ($20>=5) print $3}' | tr -d '"' > $out.valid.snp.txt ;done
for i in $(ls 05.Beta_pvalue_check/*valid.snp.txt);do out=$(basename $i | awk -F "_" '{print $1"_"$2}');/input/tools/plink2 --pfile 03.PCA30/$out --extract $i --make-pgen --out 06-1.LD_check/$out 1>06-1.LD_check/$out.log 2>06-1.LD_check/$out.err ;done

