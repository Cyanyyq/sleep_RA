#-------------------------------------------------------------------------------------------------------------------------------------#
#	QC
#-------------------------------------------------------------------------------------------------------------------------------------#
	library(data.table)
	setwd("/data/yxh/data/temp/rawdata/phenotype")
	RA.pheno = data.frame(fread("pheno_RA_time.txt",header=T))
	ind = which(RA.pheno$RA.tw>0)
	ind2 = which(RA.pheno$RA.tw<0)
	RA.diag = rep(0,502422)
	RA.diag[ind] = 1
	RA.diag[ind2] = NA
	RA.time = RA.pheno$RA.tw
	time2 = as.Date("2022-01-01")-as.Date(RA.pheno$Date)
	RA.time[-ind] = time2[-ind]
	RA.classification = data.frame(fread("RA_classification.txt",header=T))

#	RA.exclusion$exclusion = apply(RA.exclusion[,c(3:9)],1,sum)
#	RA.exclusion$exclusion[which(RA.exclusion$exclusion%in%c(1:100))] <-1
#	indd = which(RA.exclusion$exclusion==1)
#	RA.diag[indd] = NA
	RA.survival = cbind(eid=RA.pheno$eid,RA.diag,RA.time,RA.Seropositive=RA.classification[,"M05"])
	setwd("/data/yxh/data/program/PRS/Airpollution/eBMD")
	cx <- data.frame(fread("cx_new_2.txt", sep="\t",head = T))
	Pheno = cbind(RA.survival,cx[,6:45])
	setwd("/data/yxh/data/temp/rawdata/phenotype")
	sleep = data.frame(fread("pheno_sleep.txt",header=T))
	physical = data.frame(fread("pheno_physical_measures.txt",header=T,sep="\t"))
	blood_biochemistry = data.frame(fread("pheno_blood_biochemistry.txt",header=T,sep="\t"))
	sleep.1 = data.frame(sleep,physical[,c(2,10:11)],blood_biochemistry[,c(17,19,10,28)])
	cx.2 = data.frame(Pheno,sleep.1)
	setwd("/data/yxh/data/program/PRS/sleep")
	fwrite(data.frame(cx.2),"cx_new.txt", row.names = F, col.names = T,quote = F, sep = "\t")


	library(data.table)
	setwd("/data/yxh/data/temp/rawdata/phenotype")
	RA.pheno = data.frame(fread("pheno_RA_time.txt",header=T))
	ind = which(RA.pheno$RA.tw>0)
	ind2 = which(RA.pheno$RA.tw<0)
	RA.diag = rep(0,502422)
	RA.diag[ind] = 0
	RA.diag[ind2] = 1
	RA.time = RA.pheno$RA.tw
	time2 = as.Date("2022-01-01")-as.Date(RA.pheno$Date)
	RA.time[-ind] = time2[-ind]
	setwd("/data/yxh/data/temp/rawdata/phenotype")
	RA.NCIC = data.frame(fread("pheno_Non-cancer_illness_code.txt",header=T))
	indd = which(RA.NCIC[,2]%in%c("1464"))
	RA.NCIC[indd,2] = 1; RA.NCIC[-indd,2] = 0
	RF = data.frame(fread("pheno_blood_biochemistry.txt",header=T))[,"Rheumatoid_factor"]
	setwd("/data/yxh/data/program/PRS/Airpollution/eBMD")
	cx <- data.frame(fread("cx_new_2.txt", sep="\t",head = T))
	Pheno = cbind(RA.NCIC,RF,cx[,6:45])
	setwd("/data/yxh/data/temp/rawdata/phenotype")
	sleep = data.frame(fread("pheno_sleep.txt",header=T))
	physical = data.frame(fread("pheno_physical_measures.txt",header=T,sep="\t"))
	blood_biochemistry = data.frame(fread("pheno_blood_biochemistry.txt",header=T,sep="\t"))
	sleep.1 = data.frame(sleep,physical[,c(2,10:11)],blood_biochemistry[,c(17,19,10,28)])
	cx.2 = data.frame(Pheno,sleep.1)
	setwd("/data/yxh/data/program/PRS/sleep")
	fwrite(data.frame(cx.2),"cx_new_cc.txt", row.names = F, col.names = T,quote = F, sep = "\t")



	library(data.table)
	library(MatchIt)
	setwd("/data/yxh/data/program/PRS/sleep")
	cx.PRS = data.frame(fread("cx_new_cc.txt",header=T))
	cx.PRS$NCIC[which(cx.PRS$Sleep.duration%in%c(-3,-1))] = NA
	cx.PRS$NCIC[which(cx.PRS$Getting.up%in%c(-3,-1))] = NA
	cx.PRS$NCIC[which(cx.PRS$Morning.evening%in%c(-3,-1))] = NA
	cx.PRS$NCIC[which(cx.PRS$Nap.during.day%in%c(-3))] = NA
	cx.PRS$NCIC[which(cx.PRS$Sleeplessness%in%c(-3))] = NA
	cx.PRS$NCIC[which(cx.PRS$Snoring%in%c(-1,-3))] = NA
	cx.PRS$NCIC[which(cx.PRS$Daytime.dozing%in%c(-1,-3,3))] = NA
	cx.PRS$batch <- cut(cx.PRS$batch, breaks = c(-12,0,100),labels = c(1,2))
	para = c("batch","centre","sex","TDI","age","race","BMI","smoking_status","alcohol","Diet","PA")
	cx.PRS.1 = na.omit(cx.PRS[,c("eid","NCIC.time","NCIC",para)])
	para2 = colnames(cx.PRS.1)[c(4:14)]
	cx.PRS.1$NCIC = as.logical(cx.PRS.1$NCIC)
	formula <- as.formula(paste('NCIC~',paste(para2, collapse = '+')))
	set.seed(1234)
	match.it <- matchit(formula, data = cx.PRS.1, method="nearest", ratio = 4,reestimate=TRUE,verbose = TRUE)
	df.match <- match.data(match.it)
	ID = df.match$eid
	cx.PRS.2 = cx.PRS[which(cx.PRS$eid%in%ID),]
	fwrite(data.frame(cx.PRS.2),"cx_new_cc_PSM.txt", row.names = F, col.names = T,quote = F, sep = "\t")	

	cx.PRS$NCIC[which(cx.PRS$race!=1)] = NA
	para = c("batch","centre","sex","TDI","age","race","BMI","smoking_status","alcohol","Diet","PA")
	cx.PRS.1 = na.omit(cx.PRS[,c("eid","NCIC.time","NCIC",para)])
	para2 = colnames(cx.PRS.1)[c(4:14)]
	cx.PRS.1$NCIC = as.logical(cx.PRS.1$NCIC)
	formula <- as.formula(paste('NCIC~',paste(para2, collapse = '+')))
	set.seed(1234)
	match.it <- matchit(formula, data = cx.PRS.1, method="nearest", ratio = 4,reestimate=TRUE,verbose = TRUE)
	df.match <- match.data(match.it)
	ID = df.match$eid
	cx.PRS.2 = cx.PRS[which(cx.PRS$eid%in%ID),]
	fwrite(data.frame(cx.PRS.2),"cx_new_cc_PSM_eur.txt", row.names = F, col.names = T,quote = F, sep = "\t")	
	


	
	library(data.table)
	library(nnet)
	library(glmnet)
	library(pROC)	
	setwd("/data/yxh/data/temp/rawdata/phenotype")
	sleep = data.frame(fread("pheno_sleep.txt",header=T))
	physical = data.frame(fread("pheno_physical_measures.txt",header=T,sep="\t"))
	blood_biochemistry = data.frame(fread("pheno_blood_biochemistry.txt",header=T,sep="\t"))
	setwd("/data/yxh/data/temp/rawdata")
	idmatch = read.table("eid_zeng_suda_match_337198.txt",header = T)
	sleep.1 = data.frame(sleep,physical[,c(2:4,10:11)],blood_biochemistry[,c(17,19,10,28)])
	sleep.2 = merge(idmatch,sleep.1,by.x="eid.suda",by.y="eid",all.x=TRUE,sort=FALSE)
	setwd("/data/yxh/data/program/PRS/RA")
	cx.PRS = data.frame(fread("cx_PRS_new2.txt",header=T))
	cx.PRS.2 = data.frame(cx.PRS,sleep.2)

	covariate.3 = c("batch","centre","sex","TDI","age","race","BMI","smoking_status","alcohol","Diet","PA","DBP","SBP","HDL_cholesterol","LDL_direct","Cholesterol","Triglycerides","Cancer")
	setwd("/data/yxh/data/program/PRS/sleep")
	fwrite(data.frame(cx.PRS.2), paste("cx_PRS_new.txt",sep=""),quote=F,sep='\t',col.names=T,row.names=F)


#-------------------------------------------------------------------------------------------------------------------------------------#
#	Main Analysis
#-------------------------------------------------------------------------------------------------------------------------------------#
	library(data.table)
	library(survival)
	setwd("/data/yxh/data/program/PRS/sleep")
	cx.PRS = data.frame(fread("cx_new.txt",header=T))
	cx.PRS$batch <- cut(cx.PRS$batch, breaks = c(-12,0,100),labels = c(1,2))
	cx.PRS$RA.diag[which(cx.PRS$race%in%c(3:5))] = NA

	covariate.1 = c("batch","centre","sex","TDI","age","race")
	covariate.2 = c("batch","centre","sex","TDI","age","race","BMI","smoking_status","alcohol")
	covariate.3 = c("batch","centre","sex","TDI","age","race","BMI","smoking_status","alcohol","Diet","PA","DBP","SBP","HDL_cholesterol","LDL_direct","Cholesterol","Triglycerides","Cancer")
	cx.PRS$Sleep.duration[which(cx.PRS$Sleep.duration%in%c(-3,-1))] = NA
	cx.PRS$Sleep.duration = factor(cut(cx.PRS$Sleep.duration,breaks = c(0,6.5,8.5,21)),labels = c("7h-","7~8h","8h+"))
	cx.PRS$Sleep.duration = factor(cx.PRS$Sleep.duration,levels = c("7~8h","7h-","8h+"),labels = c("7~8h","7h-","8h+"))
	cx.PRS$Getting.up[which(cx.PRS$Getting.up%in%c(-3,-1))] = NA
	cx.PRS$Getting.up = factor(cx.PRS$Getting.up,levels = c(4,3,2,1),labels = c("Very easy","Fairly easy","Not very easy","Not at all easy"))
	cx.PRS$Morning.evening[which(cx.PRS$Morning.evening%in%c(-3,-1))] = NA
	cx.PRS$Morning.evening = factor(cx.PRS$Morning.evening,levels = c(1,2,3,4),labels = c("Definitely morning","Morning more","Evening more","Definitely evening"))
	cx.PRS$Nap.during.day[which(cx.PRS$Nap.during.day%in%c(-3))] = NA
	cx.PRS$Nap.during.day = factor(cx.PRS$Nap.during.day,levels = c(1,2,3),labels = c("Never/rarely","Sometimes","Usually"))
	cx.PRS$Sleeplessness[which(cx.PRS$Sleeplessness%in%c(-3))] = NA
	cx.PRS$Sleeplessness = factor(cx.PRS$Sleeplessness,levels = c(1,2,3),labels = c("Never/rarely","Sometimes","Usually"))
	cx.PRS$Snoring[which(cx.PRS$Snoring%in%c(-1,-3))] = NA
	cx.PRS$Snoring = factor(cx.PRS$Snoring,levels = c(2,1),labels = c("No","Yes"))
	cx.PRS$Daytime.dozing[which(cx.PRS$Daytime.dozing%in%c(-1,-3,3))] = NA
	cx.PRS$Daytime.dozing = factor(cx.PRS$Daytime.dozing,levels = c(0,1,2),labels = c("Never/rarely","Sometimes","Usually"))
	cx.PRS$PSS = as.numeric(cx.PRS$Sleep.duration)+as.numeric(cx.PRS$Getting.up)+as.numeric(cx.PRS$Nap.during.day)+as.numeric(cx.PRS$Sleeplessness)+as.numeric(cx.PRS$Daytime.dozing)

	tmp = summary(cx.PRS$PSS)
	cx.PRS$PSS.grade = cut(cx.PRS$PSS, breaks = c(-100,tmp[2],tmp[5],100),labels = c("Low score","Intermediate score","High score"))
	cx.PRS$PSS_quantile4 = factor(cut(cx.PRS$PSS,
				breaks = c(-100,tmp[2],tmp[3],tmp[5],100),
				labels = c("Score Q1","Score Q2","Score Q3","Score Q4")))
	cx.PRS.1 = cx.PRS
	RES = NULL
	para = c("Sleep.duration","Getting.up","Morning.evening",
		"Nap.during.day","Sleeplessness","Snoring","Daytime.dozing","PSS","PSS.grade")
	for(i in 1:length(para)){
		cx.PRS = na.omit(cx.PRS.1[,c("RA.time","RA.diag",covariate.1,para)])
		sur <- Surv(time = cx.PRS$RA.time, event = cx.PRS$RA.diag)
		formulA = as.formula(paste("sur~",para[i],"+",paste(covariate.1,collapse = "+"),sep=""))
		fit = coxph(formulA, data = cx.PRS,method="breslow")
		res = data.frame(summary(fit)$n,summary(fit)$coefficients)
		RES = rbind(RES,res)
		print(i)
	}
	fwrite(data.frame(RES),"result_sleep_model.1.csv",row.names=T)

	RES = NULL
	for(i in 1:length(para)){
		cx.PRS = na.omit(cx.PRS.1[,c("RA.time","RA.diag",covariate.2,para)])
		sur <- Surv(time = cx.PRS$RA.time, event = cx.PRS$RA.diag)
		formulA = as.formula(paste("sur~",para[i],"+",paste(covariate.2,collapse = "+"),sep=""))
		fit = coxph(formulA, data = cx.PRS,method="breslow")
		res = data.frame(summary(fit)$n,summary(fit)$coefficients)
		RES = rbind(RES,res)
		print(i)
	}
	fwrite(data.frame(RES),"result_sleep_model.2.csv",row.names=T)

	RES = NULL
	for(i in 1:length(para)){
		cx.PRS = na.omit(cx.PRS.1[,c("RA.time","RA.diag",covariate.3,para)])
		sur <- Surv(time = cx.PRS$RA.time, event = cx.PRS$RA.diag)
		formulA = as.formula(paste("sur~",para[i],"+",paste(covariate.3,collapse = "+"),sep=""))
		fit = coxph(formulA, data = cx.PRS,method="breslow")
		res = data.frame(summary(fit)$n,summary(fit)$coefficients)
		RES = rbind(RES,res)
		print(i)
	}
	fwrite(data.frame(RES),"result_sleep_model.3.csv",row.names=T)



#-------------------------------------------------------------------------------------------------------------------------------------#
#	Main Analysis RF+
#-------------------------------------------------------------------------------------------------------------------------------------#
	library(data.table)
	library(survival)
	setwd("/data/yxh/data/program/PRS/sleep")
	cx.PRS = data.frame(fread("cx_new.txt",header=T))
	cx.PRS$batch <- cut(cx.PRS$batch, breaks = c(-12,0,100),labels = c(1,2))
	cx.PRS$RA.Seropositive[which(cx.PRS$race%in%c(3:5))] = NA

	covariate.1 = c("batch","centre","sex","TDI","age","race")
	covariate.2 = c("batch","centre","sex","TDI","age","race","BMI","smoking_status","alcohol")
	covariate.3 = c("batch","centre","sex","TDI","age","race","BMI","smoking_status","alcohol","Diet","PA","DBP","SBP","HDL_cholesterol","LDL_direct","Cholesterol","Triglycerides","Cancer")
	cx.PRS$Sleep.duration[which(cx.PRS$Sleep.duration%in%c(-3,-1))] = NA
	cx.PRS$Sleep.duration = factor(cut(cx.PRS$Sleep.duration,breaks = c(0,6.5,8.5,21)),labels = c("7h-","7~8h","8h+"))
	cx.PRS$Sleep.duration = factor(cx.PRS$Sleep.duration,levels = c("7~8h","7h-","8h+"),labels = c("7~8h","7h-","8h+"))
	cx.PRS$Getting.up[which(cx.PRS$Getting.up%in%c(-3,-1))] = NA
	cx.PRS$Getting.up = factor(cx.PRS$Getting.up,levels = c(4,3,2,1),labels = c("Very easy","Fairly easy","Not very easy","Not at all easy"))
	cx.PRS$Morning.evening[which(cx.PRS$Morning.evening%in%c(-3,-1))] = NA
	cx.PRS$Morning.evening = factor(cx.PRS$Morning.evening,levels = c(1,2,3,4),labels = c("Definitely morning","Morning more","Evening more","Definitely evening"))
	cx.PRS$Nap.during.day[which(cx.PRS$Nap.during.day%in%c(-3))] = NA
	cx.PRS$Nap.during.day = factor(cx.PRS$Nap.during.day,levels = c(1,2,3),labels = c("Never/rarely","Sometimes","Usually"))
	cx.PRS$Sleeplessness[which(cx.PRS$Sleeplessness%in%c(-3))] = NA
	cx.PRS$Sleeplessness = factor(cx.PRS$Sleeplessness,levels = c(1,2,3),labels = c("Never/rarely","Sometimes","Usually"))
	cx.PRS$Snoring[which(cx.PRS$Snoring%in%c(-1,-3))] = NA
	cx.PRS$Snoring = factor(cx.PRS$Snoring,levels = c(2,1),labels = c("No","Yes"))
	cx.PRS$Daytime.dozing[which(cx.PRS$Daytime.dozing%in%c(-1,-3,3))] = NA
	cx.PRS$Daytime.dozing = factor(cx.PRS$Daytime.dozing,levels = c(0,1,2),labels = c("Never/rarely","Sometimes","Usually"))
	cx.PRS$PSS = as.numeric(cx.PRS$Sleep.duration)+as.numeric(cx.PRS$Getting.up)+as.numeric(cx.PRS$Nap.during.day)+as.numeric(cx.PRS$Sleeplessness)+as.numeric(cx.PRS$Daytime.dozing)

	ind = which(cx.PRS$RA.diag==1&cx.PRS$RA.Seropositive==0)
	cx.PRS$RA.Seropositive[ind] = NA
	tmp = summary(cx.PRS$PSS)
	cx.PRS$PSS.grade = cut(cx.PRS$PSS, breaks = c(-100,tmp[2],tmp[5],100),labels = c("Low score","Intermediate score","High score"))
	cx.PRS$PSS_quantile4 = factor(cut(cx.PRS$PSS,
				breaks = c(-100,tmp[2],tmp[3],tmp[5],100),
				labels = c("Score Q1","Score Q2","Score Q3","Score Q4")))
	cx.PRS.1 = cx.PRS
	RES = NULL
	para = c("Sleep.duration","Getting.up","Morning.evening",
		"Nap.during.day","Sleeplessness","Snoring","Daytime.dozing","PSS","PSS.grade")
	for(i in 1:length(para)){
		cx.PRS = na.omit(cx.PRS.1[,c("RA.time","RA.Seropositive",covariate.1,para)])
		sur <- Surv(time = cx.PRS$RA.time, event = cx.PRS$RA.Seropositive)
		formulA = as.formula(paste("sur~",para[i],"+",paste(covariate.1,collapse = "+"),sep=""))
		fit = coxph(formulA, data = cx.PRS,method="breslow")
		res = data.frame(summary(fit)$n,summary(fit)$coefficients)
		RES = rbind(RES,res)
		print(i)
	}
	fwrite(data.frame(RES),"result_sleep_model.1_RF+.csv",row.names=T)

	RES = NULL
	for(i in 1:length(para)){
		cx.PRS = na.omit(cx.PRS.1[,c("RA.time","RA.Seropositive",covariate.2,para)])
		sur <- Surv(time = cx.PRS$RA.time, event = cx.PRS$RA.Seropositive)
		formulA = as.formula(paste("sur~",para[i],"+",paste(covariate.2,collapse = "+"),sep=""))
		fit = coxph(formulA, data = cx.PRS,method="breslow")
		res = data.frame(summary(fit)$n,summary(fit)$coefficients)
		RES = rbind(RES,res)
		print(i)
	}
	fwrite(data.frame(RES),"result_sleep_model.2_RF+.csv",row.names=T)

	RES = NULL
	for(i in 1:length(para)){
		cx.PRS = na.omit(cx.PRS.1[,c("RA.time","RA.Seropositive",covariate.3,para)])
		sur <- Surv(time = cx.PRS$RA.time, event = cx.PRS$RA.Seropositive)
		formulA = as.formula(paste("sur~",para[i],"+",paste(covariate.3,collapse = "+"),sep=""))
		fit = coxph(formulA, data = cx.PRS,method="breslow")
		res = data.frame(summary(fit)$n,summary(fit)$coefficients)
		RES = rbind(RES,res)
		print(i)
	}
	fwrite(data.frame(RES),"result_sleep_model.3_RF+.csv",row.names=T)


#-------------------------------------------------------------------------------------------------------------------------------------#
#	PRS Analysis
#-------------------------------------------------------------------------------------------------------------------------------------#
	library(data.table)
	library(survival)
	setwd("/data/yxh/data/program/PRS/sleep")
	cx.PRS = data.frame(fread("cx_PRS_new.txt",header=T))
	cx.PRS$RA.diag[which(cx.PRS$race%in%c(3:5))] = NA
	covariate.1 = c("batch","centre","sex","TDI","age")
	covariate.2 = c("batch","centre","sex","TDI","age","BMI","smoke","alcohol")
	covariate.3 = c("batch","centre","sex","TDI","age","BMI","smoke","alcohol","diet_healthy","Physical_activity","DBP","SBP","HDL_cholesterol","LDL_direct","Cholesterol","Triglycerides")
	cx.PRS$Sleep.duration[which(cx.PRS$Sleep.duration%in%c(-3,-1))] = NA
	cx.PRS$Sleep.duration = factor(cut(cx.PRS$Sleep.duration,breaks = c(0,6.5,8.5,21)),labels = c("7h-","7~8h","8h+"))
	cx.PRS$Sleep.duration = factor(cx.PRS$Sleep.duration,levels = c("7~8h","7h-","8h+"),labels = c("7~8h","7h-","8h+"))
	cx.PRS$Getting.up[which(cx.PRS$Getting.up%in%c(-3,-1))] = NA
	cx.PRS$Getting.up = factor(cx.PRS$Getting.up,levels = c(4,3,2,1),labels = c("Very easy","Fairly easy","Not very easy","Not at all easy"))
	cx.PRS$Morning.evening[which(cx.PRS$Morning.evening%in%c(-3,-1))] = NA
	cx.PRS$Morning.evening = factor(cx.PRS$Morning.evening,levels = c(1,2,3,4),labels = c("Definitely morning","Morning more","Evening more","Definitely evening"))
	cx.PRS$Nap.during.day[which(cx.PRS$Nap.during.day%in%c(-3))] = NA
	cx.PRS$Nap.during.day = factor(cx.PRS$Nap.during.day,levels = c(1,2,3),labels = c("Never/rarely","Sometimes","Usually"))
	cx.PRS$Sleeplessness[which(cx.PRS$Sleeplessness%in%c(-3))] = NA
	cx.PRS$Sleeplessness = factor(cx.PRS$Sleeplessness,levels = c(1,2,3),labels = c("Never/rarely","Sometimes","Usually"))
	cx.PRS$Snoring[which(cx.PRS$Snoring%in%c(-1,-3))] = NA
	cx.PRS$Snoring = factor(cx.PRS$Snoring,levels = c(2,1),labels = c("No","Yes"))
	cx.PRS$Daytime.dozing[which(cx.PRS$Daytime.dozing%in%c(-1,-3,3))] = NA
	cx.PRS$Daytime.dozing = factor(cx.PRS$Daytime.dozing,levels = c(0,1,2),labels = c("Never/rarely","Sometimes","Usually"))
	cx.PRS$PSS = as.numeric(cx.PRS$Sleep.duration)+as.numeric(cx.PRS$Getting.up)+as.numeric(cx.PRS$Nap.during.day)+as.numeric(cx.PRS$Sleeplessness)+as.numeric(cx.PRS$Daytime.dozing)
	tmp = summary(cx.PRS$PSS)
	cx.PRS$PSS.grade = cut(cx.PRS$PSS, breaks = c(-100,tmp[2],tmp[5],100),labels = c("Low score","Intermediate score","High score"))


	n = dim(cx.PRS)[1]
	set.seed(12345)
	index = sample(1:n,0.2*n,rep=F)
	indd.eid = cx.PRS$eid[-index]
	indd.eid2 = cx.PRS$eid[index]
	cx.PRS.train = cx.PRS[which(cx.PRS$eid%in%indd.eid2),]
	cx.PRS.test = cx.PRS[which(cx.PRS$eid%in%indd.eid),]
	tmp = summary(cx.PRS.test$GRS_withMHC_5E6)
	cx.PRS.test$PRS = cut(cx.PRS.test$GRS_withMHC_5E6, breaks = c(-100,tmp[2],tmp[5],100),labels = c("Low risk","Intermediate risk","High risk"))
	
	cx.PRS.test$group = paste(cx.PRS.test$PRS,cx.PRS.test$PSS.grade,sep="_")
	cx.PRS.test$group[which(cx.PRS.test$group%in%c("High risk_NA","Intermediate risk_NA","Low risk_NA"))] = NA
	cx.PRS.test$group = factor(cx.PRS.test$group,levels = c("Low risk_Low score","Low risk_Intermediate score","Low risk_High score",
					"Intermediate risk_Low score","Intermediate risk_Intermediate score","Intermediate risk_High score",
					"High risk_Low score","High risk_Intermediate score","High risk_High score"))
	cx.PRS.test.male = cx.PRS.test[which(cx.PRS.test$sex==1),]
	cx.PRS.test.female = cx.PRS.test[which(cx.PRS.test$sex==0),]

	#	RA~PRS+PSS
	#-------------------------------------------------------------------------------------------------------------------------------#
		covariate = c("batch","centre","age","sex","TDI","PC.1","PC.2","PC.3","PC.4","PC.5","PC.6","PC.7","PC.8","PC.9","PC.10")
		sur <- Surv(time = cx.PRS.test$RA.time, event = cx.PRS.test$RA.diag)	
		formulA = as.formula(paste("sur~group+",paste(covariate,collapse = "+"),sep=""))
		fit = coxph(formulA, data = cx.PRS.test,method="breslow")
		res = summary(fit)$coefficients
		sur <- Surv(time = cx.PRS.test.male$RA.time, event = cx.PRS.test.male$RA.diag)	
		fit.male = coxph(formulA, data = cx.PRS.test.male,method="breslow")
		res.male = summary(fit.male)$coefficients
		sur <- Surv(time = cx.PRS.test.female$RA.time, event = cx.PRS.test.female$RA.diag)	
		fit.female = coxph(formulA, data = cx.PRS.test.female,method="breslow")
		res.female = summary(fit.female)$coefficients
		RES = rbind(cbind(summary(fit)$n,res),cbind(summary(fit.male)$n,res.male),cbind(summary(fit.female)$n,res.female))
		setwd("/data/yxh/data/program/PRS/sleep")
		write.csv(RES, file = "cox_result_PSS+PRS.csv")

	#	RA~PRS*PSS
	#-------------------------------------------------------------------------------------------------------------------------------#
		covariate = c("batch","centre","age","sex","TDI","PC.1","PC.2","PC.3","PC.4","PC.5","PC.6","PC.7","PC.8","PC.9","PC.10")
		sur <- Surv(time = cx.PRS.test$RA.time, event = cx.PRS.test$RA.diag)	
		formulA = as.formula(paste("sur~GRS_withMHC_5E6*PSS+",paste(covariate,collapse = "+"),sep=""))
		fit = coxph(formulA, data = cx.PRS.test,method="breslow")
		res = summary(fit)$coefficients
		sur <- Surv(time = cx.PRS.test.male$RA.time, event = cx.PRS.test.male$RA.diag)	
		fit.male = coxph(formulA, data = cx.PRS.test.male,method="breslow")
		res.male = summary(fit.male)$coefficients
		sur <- Surv(time = cx.PRS.test.female$RA.time, event = cx.PRS.test.female$RA.diag)	
		fit.female = coxph(formulA, data = cx.PRS.test.female,method="breslow")
		res.female = summary(fit.female)$coefficients
		RES = rbind(cbind(summary(fit)$n,res),cbind(summary(fit.male)$n,res.male),cbind(summary(fit.female)$n,res.female))
		setwd("/data/yxh/data/program/PRS/sleep")
		write.csv(RES, file = "cox_result_PSS:PRS.csv")





#-------------------------------------------------------------------------------------------------------------------------------------#
#	Main Analysis Case control
#-------------------------------------------------------------------------------------------------------------------------------------#
	library(data.table)
	library(survival)
	setwd("/data/yxh/data/program/PRS/sleep")
	cx.PRS = data.frame(fread("cx_new_cc_PSM.txt",header=T))

	covariate.1 = c("batch","sex","TDI","scale(age)","race")
	covariate.2 = c("batch","sex","TDI","age","race","BMI","smoking_status","alcohol")
	covariate.3 = c("batch","sex","TDI","age","race","BMI","smoking_status","alcohol","Diet","PA","DBP","SBP","HDL_cholesterol","LDL_direct","Cholesterol","Triglycerides","Cancer")
	cx.PRS$Getting.up = factor(cx.PRS$Getting.up,levels = c(4,3,2,1),labels = c("Very easy","Fairly easy","Not very easy","Not at all easy"))
	cx.PRS$Morning.evening = factor(cx.PRS$Morning.evening,levels = c(1,2,3,4),labels = c("Definitely morning","Morning more","Evening more","Definitely evening"))
	cx.PRS$Nap.during.day = factor(cx.PRS$Nap.during.day,levels = c(1,2,3),labels = c("Never/rarely","Sometimes","Usually"))
	cx.PRS$Sleeplessness = factor(cx.PRS$Sleeplessness,levels = c(1,2,3),labels = c("Never/rarely","Sometimes","Usually"))
	cx.PRS$Snoring = factor(cx.PRS$Snoring,levels = c(2,1),labels = c("No","Yes"))
	cx.PRS$Daytime.dozing = factor(cx.PRS$Daytime.dozing,levels = c(0,1,2),labels = c("Never/rarely","Sometimes","Usually"))
	cx.PRS$Sleep.duration.1 = factor(cut(cx.PRS$Sleep.duration,breaks = c(0,6.5,8.5,21)),labels = c("7h-","7~8h","8h+"))
	cx.PRS$Sleep.duration.1 = factor(cx.PRS$Sleep.duration.1,levels = c("7~8h","7h-","8h+"),labels = c("7~8h","no 7~8h","no 7~8h"))
	cx.PRS$PSS = as.numeric(cx.PRS$Sleep.duration.1)+as.numeric(cx.PRS$Getting.up)+as.numeric(cx.PRS$Nap.during.day)+as.numeric(cx.PRS$Sleeplessness)+as.numeric(cx.PRS$Daytime.dozing)
	tmp = summary(cx.PRS$PSS)
	cx.PRS$PSS.grade = cut(cx.PRS$PSS, breaks = c(-100,tmp[2],tmp[5],100),labels = c("Low score","Intermediate score","High score"))
	cx.PRS$PSS_quantile4 = factor(cut(cx.PRS$PSS,
				breaks = c(-100,tmp[2],tmp[3],tmp[5],100),
				labels = c("Score Q1","Score Q2","Score Q3","Score Q4")))


	library(MASS)
	para = c("Getting.up","Morning.evening","Nap.during.day",
		"Sleeplessness","Daytime.dozing","PSS.grade")
	covariate.1 = c("batch","sex","TDI","scale(age)","race")
	RES = NULL
	for(i in 1:length(para)){
		formulA = as.formula(paste(para[i],"~NCIC+",paste(covariate.1,collapse = "+"),sep=""))
		fit <- polr(formulA, data = cx.PRS, Hess=TRUE)
		ctable <- coef(summary(fit))
		p <- pnorm(abs(ctable[, "t value"]), lower.tail = FALSE)*2
		ctable <- data.frame(ctable, "p value" = p)
		ctable$OR = exp(ctable$Value)
		ctable$LOW = exp(ctable$Value-1.96*ctable$Std..Error)
		ctable$UP = exp(ctable$Value+1.96*ctable$Std..Error)
		res = data.frame(para[i],ctable[1,])
		RES = rbind(RES,res)
		print(i)
	}
	RES.1 = RES[,c("para.i.","OR","LOW","UP","p.value")]
	colnames(RES.1) = c("para","estimate","LOW","UP","P"); rownames(RES.1)=NULL
	
	RES = NULL
	covariate.1 = c("batch","sex","TDI","age","race")
	para = c("Snoring")
	for(i in 1:length(para)){
		formulA = as.formula(paste(para[i],"~NCIC+",paste(covariate.1,collapse = "+"),sep=""))
		fit = glm(formulA, data = cx.PRS, family = "binomial")
		ctable <- data.frame(summary(fit)$coefficients)
		ctable$LOW = exp(ctable$Estimate-1.96*ctable$Std..Error)
		ctable$UP = exp(ctable$Estimate+1.96*ctable$Std..Error)
		res = data.frame(para[i],ctable["NCIC",c("Estimate","LOW","UP","Pr...z..")])
		RES = rbind(RES,res)
		print(i)
	}
	RES.2 = RES
	colnames(RES.2) = c("para","estimate","LOW","UP","P"); rownames(RES.2)=NULL


	RES = NULL
	covariate.1 = c("batch","sex","TDI","age","race")
	para = c("Sleep.duration")
	for(i in 1:length(para)){
		formulA = as.formula(paste(para[i],"~NCIC+",paste(covariate.1,collapse = "+"),sep=""))
		fit = lm(formulA, data = cx.PRS)
		ctable <- data.frame(summary(fit)$coefficients)
		ctable$LOW = exp(ctable$Estimate-1.96*ctable$Std..Error)
		ctable$UP = exp(ctable$Estimate+1.96*ctable$Std..Error)
		res = data.frame(para[i],ctable["NCIC",c("Estimate","LOW","UP","Pr...t..")])
		RES = rbind(RES,res)
		print(i)
	}
	RES.3 = RES
	colnames(RES.3) = c("para","estimate","LOW","UP","P"); rownames(RES.3)=NULL
	RES = rbind(RES.1,RES.2,RES.3)
	fwrite(data.frame(RES),"result_CC.csv",row.names=T)
