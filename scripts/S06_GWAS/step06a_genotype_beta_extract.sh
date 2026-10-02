
#extract only rsid, not the "." in rsid
for i in $(ls *GT*txt);do out=$(basename $i | awk -F "." '{print "chr1-22.250227."$4"."$3"_0.1M"}');cat $i | awk -F " " '{if ($1!=".") print $0}' > /input/data/train_test_split/06.Genotype_beta_convert/$out.GT.needtocheck.txt ;done
