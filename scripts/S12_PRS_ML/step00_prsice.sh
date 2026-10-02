awk -F " " '{if ($6<=5e-08) print $0}' /input/data/UC526KDC.03.ff.b38con.PHENO1.glm.logistic.cc.assoc > UC526.03.ff.b38con.5e-08.assoc
le UC526.03.ff.b38con.5e-08.assoc
sed -i '1i SNP CHR BP A1 A2 P OR' UC526.03.ff.b38con.5e-08.assoc
head UC526.03.ff.b38con.5e-08.assoc
/input/tools/PRSice.R --dir . --prsice /input/tools/PRSice --base 227loci.CD.assoc --target 03.train_test_split/CD238KDC/CD238combine.grp5.tst.pc30.b38c.ff --thread 10 --stat OR --no-clump --logit-perm --bar-levels 1,0.5,0.2,0.1,1e-02,1e-03,1e-04,1e-05,1e-06,5e-08 --fastscore --print-snp --out 227.lociCD --score sum 1>227.lociCD.log 2>227.lociCD.err
/input/tools/PRSice.R --dir . --prsice /input/tools/PRSice --base 227loci.UC.assoc --target 03.train_test_split/UC526KDC/UC526KDC.grp1.tst.pc30.b38c.ff --thread 10 --stat OR --no-clump --logit-perm --bar-levels 1,0.5,0.2,0.1,1e-02,1e-03,1e-04,1e-05,1e-06,5e-08 --fastscore --print-snp --out 227.lociUC --score sum 1>227.lociUC.log 2>227.lociUC.err


/input/tools/PRSice.R --dir . --prsice /input/tools/PRSice --base check.log --target CD238KDC/CD238combine.grp5.tst.pc30.b38c.ff --thread 10 --stat OR --no-clump --logit-perm --bar-levels 1,0.5,0.2,0.1,1e-02,1e-03,1e-04,1e-05,1e-06,5e-08 --fastscore --print-snp --out 44.lociCD --score sum 1>44.lociCD.log 2>44.lociCD.err
/input/tools/PRSice.R --dir . --prsice /input/tools/PRSice --base UC526KDc.check.5e-08.real.assoc --target UC526KDC/UC526KDC.grp1.tst.pc30.b38c.ff --thread 10 --stat OR --no-clump --logit-perm --bar-levels 1,0.5,0.2,0.1,1e-02,1e-03,1e-04,1e-05,1e-06,5e-08 --fastscore --print-snp --out 44.lociUC --score sum 1>44.lociUC.log 2>44.lociUC.err
