for i in $(ls 04.firthfall/*logistic.hybrid);do out=$(basename $i | cut -d "." -f 1,2,3,4);Rscript 05.Beta_pvalue_check/Rscript.beta_rankP.plus.R $i 05.Beta_pvalue_check/$out ;done
