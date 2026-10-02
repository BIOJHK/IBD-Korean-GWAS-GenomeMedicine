for i in $(ls 05.Beta_pvalue_check/*valid.snp.txt);do out=$(basename $i | awk -F "_" '{print $1"_"$2}');/input/tools/plink2 --pfile 06-1.LD_check/$out --r2-phased --out 06-1.LD_check/$out.LD_check 1>06-1.LD_check/$out.LD_check.log 2>06-1.LD_check/$out.LD_check.err ;done
for i in $(ls 05.Beta_pvalue_check/*valid.snp.txt);do out=$(basename $i | awk -F "_" '{print $1"_"$2}');python vcor_LD.change.step062.py 06-1.LD_check/$out*.vcor 06-1.LD_check/$out.LD.txt  ;done

