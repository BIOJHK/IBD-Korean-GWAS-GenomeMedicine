for i in $(ls *bcf | grep -v "0.1M");do out=$(basename $i | cut -d "." -f 1,2,3,4);bcftools query -l $out.bcf | awk -F " " '{if ($1 ~ /^KDC/) $2 = "1"; else $2 = "2"; print $0}' > $out.pheno.txt ;done
#after this process, nano pheno.txt --> IID DIS typing in first line
