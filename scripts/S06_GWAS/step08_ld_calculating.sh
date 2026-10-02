#susieR
for i in $(ls 05.Beta_pvalue_check/*_0.1M.rankPbeta+.hybrid);do out=$(basename $i | awk -F "us" '{print $1"us"}');Rscript Rscript.LD_create_and_susieR.R 06-1.LD_check/$out*.vcor $i 06-1.LD_check/$out ;done
[1] "data importing process was done"
[1] "ld matrix making process was done"
[1] "cs1 값은 다음과 같습니다:"
[1] 28 29 30 31 32
[1] "susieR process was done"
[1] "data importing process was done"
[1] "ld matrix making process was done"
[1] "cs1 값은 다음과 같습니다:"
[1] 18 19 20 21 22 23
[1] "susieR process was done"
[1] "data importing process was done"
[1] "ld matrix making process was done"
[1] "cs1 값은 다음과 같습니다:"
[1] 4
[1] "susieR process was done"
[1] "data importing process was done"
[1] "ld matrix making process was done"
[1] "cs1 값은 다음과 같습니다:"
[1] 33
[1] "susieR process was done"
[1] "data importing process was done"
[1] "ld matrix making process was done"
Warning message:
In susie(as.matrix(LD), Z) :
  IBSS algorithm did not converge in 100 iterations!
[1] "cs1 값은 다음과 같습니다:"
[1] 50
[1] "susieR process was done"
[1] "data importing process was done"
[1] "ld matrix making process was done"
[1] "cs1 값은 다음과 같습니다:"
 [1]  60  62  63  64  65  66  67  74  76  77  78  79  80  86 106
[1] "susieR process was done"
[1] "data importing process was done"
[1] "ld matrix making process was done"
Warning message:
In susie(as.matrix(LD), Z) :
  IBSS algorithm did not converge in 100 iterations!
[1] "cs1 값은 다음과 같습니다:"
 [1] 128 130 132 133 135 139 141 142 143 144 146 147 150 152 154 155 178 179 181
[20] 184 193 194 195 196 197 198 232 233 300 304 322 332 333 335 337 338 341 342
[39] 343 347 348 349 350 351 352 353 360 361 370 376 377 446 447 448 456 467 490
[58] 493 494 501 524 601
[1] "susieR process was done"
[1] "data importing process was done"
[1] "ld matrix making process was done"
[1] "cs1 값은 다음과 같습니다:"
[1] 830 833
[1] "susieR process was done"

#Finemapping

#개별 데이터 고치기
#개별 데이터 고치기
z.for.z
for i in $(ls 05.Beta_pvalue_check/*_0.1M.rankPbeta+.hybrid);do out=$(basename $i | awk -F "us" '{print $1"us"}');cat $i | awk -F " " '{if ($20>=5) print $3,$1,$2,$4,$5,$9,$21,$14}' > 06-1.LD_check/$out.250304/$out.zfor.z ;done
snp.snp
chr1-22.250227.test.KDCCD_chr7minus_ld.250228.txt to chr1-22.250227.test.KDCCD_chr7minus_ld.250228.ld
for i in $(ls 05.Beta_pvalue_check/*_0.1M.rankPbeta+.hybrid);do out=$(basename $i | awk -F "us" '{print $1"us"}');cp 06-1.LD_check/$out"_susieR.snp.250228.txt" 06-1.LD_check/$out.250304/ ;done
cred
for i in $(ls 05.Beta_pvalue_check/*_0.1M.rankPbeta+.hybrid);do out=$(basename $i | awk -F "us" '{print $1"us"}');cat $i | awk -F " " '{if ($20>=5) print $3,$18}' > 06-1.LD_check/$out.250304/$out.cs1.cred ;done

master making
master ::
z;ld;snp;config;cred;log;n_samples
chr1-22.250227.test.KDCCD_chr7minus.zfor.z;chr1-22.250227.test.KDCCD_chr7minus_ld.250228.ld;chr1-22.250227.test.KDCCD_chr7minus.snp.snp;chr1-22.250227.test.KDCCD_chr7minus.config;chr1-22.250227.test.KDCCD_chr7minus.cred;chr1-22.250227.test.KDCCD_chr7minus.log;50
finemap --sss --in-files chr1-22.250227.test.KDCCD_chr7minus.master

#CaVIAR
for i in $(ls 05.Beta_pvalue_check/*_0.1M.rankPbeta+.hybrid);do out=$(basename $i | awk -F "us" '{print $1"us"}');cat 05.Beta_pvalue_check/chr1-22.250227.test.KDCCD_chr7minus_250227_0.1M.rankPbeta+.hybrid | awk -F " " '{if ($20>=5) print $3,$17}' > 06-1.LD_check/$out.zstat.value.txt ;done
for i in $(ls 05.Beta_pvalue_check/*_0.1M.rankPbeta+.hybrid);do out=$(basename $i | awk -F "us" '{print $1"us"}');CAVIAR -o 06-1.LD_check/$out.CAVIAR.snp -z 06-1.LD_check/$out.zstat.value.txt -l 06-1.LD_check/$out"_ld.250228.txt" 1>06-1.LD_check/$out.CAVIAR.snp.log 2>06-1.LD_check/$out.CAVIAR.snp.err ;done
