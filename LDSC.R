

cd /data/yxh/data/ref/1000G_Phase3_baselineLD_v2.2_ldscores
rename baselineLD.*.l2.ldscore.gz *.l2.ldscore.gz ?.l2.ldscore.gz

setwd("/data/yxh/data/program/LDSC/Alltraits")
filename = list.files("/data/yxh/data/program/LDSC/Alltraits")
write.table(filename,"filename.txt",row.names=F,quote=F)

# LDSC genetic correlation
#-----------------------------------------------------------------------------------

	bash
	source activate ldsc
	for dx in ldsc_Sleep.duration ldsc_Getting.up ldsc_Morning.evening ldsc_Nap.during.day ldsc_Sleeplessness ldsc_Snoring ldsc_Daytime.dozing ; do
	for ex in ldsc_RA; do
	ldsc=/home/yxh/software/ldsc/ldsc.py
	cd /data/yxh/data/program/LDSC/result
	home=/data/yxh/data/program/LDSC/Alltraits
	home2=/data/yxh/data/program/LDSC
	${ldsc} --rg ${home}/${dx}.txt,${home2}/${ex}.txt \
		--ref-ld-chr /data/yxh/data/ref/1000G_Phase3_baselineLD_v2.2_ldscores/ \
		--w-ld-chr   /data/yxh/data/ref/1000G_Phase3_weights_hm3_no_MHC/weights.hm3_noMHC. \
		--frqfile-chr 1000G_Phase3_frq/1000G.EUR.QC. \
		--chisq-max 5000.0 \
		--out ./${dx}_${ex}
	done
	done

