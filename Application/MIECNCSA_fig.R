source("/Users/kararudolph/Documents/JHSPHpostdoc/measurementErrorImputation/MIEC/MI-EC algorithm.r")
source("/Users/kararudolph/Documents/JHSPHpostdoc/measurementErrorImputation/MIEC/FunctionsMIECcombUrbanSubsetATT.R")
library(mitools)
library(sandwich)
library(ggplot2)

setwd("/Users/kararudolph/Documents/JHSPHpostdoc/measurementErrorImputation")
validationcomp<-read.csv("validationcomp.csv", header=TRUE)
urb<-validationcomp[validationcomp$urbancat==1,]
validationcomp<-urb
## What if the correlations are 0.9, 0.7, and 0.4 and the validation sample is smaller--10% of main sample: N=416.
## using the subset of urban, there are 400 in the calibration sample and 1526 in the main sample
set.seed(12309)
calibration<-validationcomp[, c("cmatAge", "cmatAgeworse", "cmatAgeevenworse", "mmatAge", "SampleID")]
smallcalib<-calibration[sample(nrow(calibration), 400), ]
var(smallcalib$cmatAge)
var(smallcalib$mmatAge)
var(validationcomp$cmatAge)
var(validationcomp$mmatAge)
cor(smallcalib$cmatAge, smallcalib$mmatAge)
cor(smallcalib$cmatAgeworse, smallcalib$mmatAge)
cor(smallcalib$cmatAgeevenworse, smallcalib$mmatAge)
validationcomp$substance<-ifelse(validationcomp$d_drug12_NIMH2==1 | validationcomp$d_alcohol12_NIMH2==1, 1, 0)

mainsample<-validationcomp[!validationcomp$SampleID %in% smallcalib$SampleID, c("cmatAge", "cmatAgeworse", "cmatAgeevenworse", "d_drug12_NIMH2", "cp_CdOddh12_NIMH2", "tertscore", "SEXF", "age_cent", "urbancat", "suburb", "black", "latino", "other", "midwest", "south", "west", "cinc", "mmatAge", "fath", "moth", "nohighschool", "highschool", "smcollege", "pc_pa_minor", "pc_pa_severe", "pc_psych_minor", "d_eating12_NIMH2", "d_alcohol12_NIMH2", "imgen", "d_mdddys12_NIMH2" , "d_anxiety12_NIMH2", "any", "internal", "substance")]

table(mainsample$internal)
table(mainsample$d_drug12_NIMH2)
table(I(mainsample$d_drug12_NIMH2==1 | mainsample$d_alcohol12_NIMH2==1))
table(mainsample$substance)
#urbanvalid<-validationcomp[validationcomp$urbancat==1,]
#urbancalibration<-urbanvalid[, c("cmatAge", "cmatAgeworse", "cmatAgeevenworse", "mmatAge", "SampleID")]
#urbansmallcalib<-urbancalibration[sample(nrow(urbancalibration), 479), ]
#mainsampleurban<-urbanvalid[!urbanvalid$SampleID %in% urbansmallcalib$SampleID, c("cmatAge", "cmatAgeworse", "cmatAgeevenworse", "d_drug12_NIMH2", "cp_CdOddh12_NIMH2", "tertscore", "SEXF", "age_cent", "urbancat", "suburb", "black", "latino", "other", "midwest", "south", "west", "cinc", "mmatAge", "fath", "moth", "nohighschool", "highschool", "smcollege", "pc_pa_minor", "pc_pa_severe", "pc_psych_minor", "d_eating12_NIMH2", "d_alcohol12_NIMH2", "imgen", "d_mdddys12_NIMH2" , "d_anxiety12_NIMH2", "any", "internal")]

confounders<-c("SEXF", "age_cent", "moth", "fath", "black", "latino", "other", "midwest", "south", "west", "cinc")

#mainsample<-validationcomp[!validationcomp$SampleID %in% smallcalib$SampleID, c("cmatAge", "cmatAgeworse", "cmatAgeevenworse", "d_drug12_NIMH2", "tertscore", "SEXF", "age_cent", "urbancat", "suburb", "black", "latino", "other", "midwest", "south", "west", "cinc", "mmatAge", "nohighschool", "highschool", "smcollege")]

#confounders<-c("SEXF", "age_cent", "urbancat", "suburb", "black", "latino", "other", "midwest", "south", "west", "cinc", "nohighschool", "highschool", "smcollege")


m<-12
n<-3

n_calib<-nrow(smallcalib)
n_main<-nrow(mainsample)

T.outcomes<-c("tertscore", "substance")
Z.covariates<-confounders

q<-length(T.outcomes)
r<-length(Z.covariates)

set.seed(2394022)


## get treatment effect estimate for 1) truth, 2) 0.9 correlation W, 3) 0.7 correlation W, 4) 0.4 correlation W
MIEC_data_point9 <- MIEC(mainsample[,c("cmatAge", T.outcomes, Z.covariates)], smallcalib[,c("cmatAge", "mmatAge")], n_calib,n_main,M=m,N=n,K=q,S=r)

lovplotgen<-function(data){
  c<-data.frame(data[,c(1,3:ncol(data))])
  means<-ddply(c, .(tertscore), colMeans)
  vars<-ddply(c, .(tertscore), colwise(var))
  stdmeandif<-(unlist(means[2,]) - unlist(means[1,]))/sqrt(unlist(vars[2,]))
  return(c(stdmeandif))
}

pre<-lovplotgen(MIEC_data_point4)

stdmeandif<-matrix(rep(NA, 36*12),ncol=12)
wt<-list(rep(NA))

ddply(c, .(tertscore), function(x)
  )
wtlovplotgen<-function(data){
 # data<-MIEC_data_point9
  for(i in 14:ncol(data)){
 #   i<-14
    txmodel <- glm(data[,1] ~ data[,i]+ data[,3] +  data[,4] +  data[,5] +  data[,6] +  data[,7] +  data[,8] +  data[,9] +  data[,10] +  data[,11] +  data[,12] +  data[,13] , family = "binomial") 
    ps<- predict(txmodel, type="response")
    wt[[i]]<-ifelse(data[,1] ==1, 1, ps/(1-ps))
  c<-data.frame(data)
  c$wt<-wt[[i]]
  des<-svydesign(ids = ~1, weights = ~wt, data = c)
  means<-svyby(make.formula(names(c)[c(3:13,i)]) ,~tertscore, des, svymean)
  #vars<-svyby(make.formula(names(c)[c(3:13,i)]) ,~tertscore, des, svyvar)
  deff<-sum(c$wt[c$tertscore==1])/((sum(c$wt[c$tertscore==1]))^2 - sum(c$wt[c$tertscore==1]^2))
  for(j in 2:13){
    varout[j-1]<-deff*sum((c[c$tertscore==1,j+1]-means[2,j])^2)
  }
  
  stdmeandif[i-13,]<-(unlist(means[2,c(2:13)]) - unlist(means[1,c(2:13)]))/varout
}
  return(stdmeandif)
}

post<-wtlovplotgen(MIEC_data_point4)
post2<-cbind(post[,11],post[,5], post[,4], post[,6], post[,2], post[,8], post[,9], post[,3], post[,1], post[,7], post[,10], post[,12])
post<-post2
pre<-abs(pre)
nam2<-c("household income", "black", "live with father", "latino","age", "midwest", "south", "live with mother", "gender", "other race",   "west", "maternal age at birth")
orderedpear<-pre[2:12][order(pre[2:12], decreasing=TRUE)]
orderedpear1<-matrix(c(rep(orderedpear, 36)), byrow=TRUE, nrow=36)
orderedpear2<-cbind(orderedpear1, pre[13:48])
orderedpear3<-apply(orderedpear2, c(1,2), function(x) as.numeric(x)*100)
peach1<-apply(post, c(1,2), function(x) abs(as.numeric(x)*100))

pdf("average.pdf")
dotchart(colMeans(peach1), pch="", labels=nam2, cex=0.75) 
mtext("Standardized Difference (%)", side=1, line=2)
points(colMeans(orderedpear3), seq(1:12), 
pch=21, col="blue", cex=1.2) 
points(colMeans(peach1),seq(1:12), 
pch=16, col="red", cex=1.2) 
abline(v=0, lty=1) 
abline(v=10, lty=2, lwd=1, col="grey") 
abline(v=20, lty=2, lwd=1, col="black") 
legend("topright", legend = c("Unweighted", "Weighted"), col=c("blue", "red"), text.col=c("blue", "red"), pch=c(21,16))
dev.off()

pdf("onetonine.pdf")
par(mfrow=c(3,3))
for(i in 1:9){
dotchart(peach1[i,], pch="", labels=nam2, cex=0.75) 
mtext("Standardized Difference (%)", side=1, line=2)
points(orderedpear3[i,], seq(1:length(orderedpear3[i,])), 
pch=21, col="blue", cex=1.2) 
points(peach1[i,],seq(1:length(peach1[i,])), 
pch=16, col="red", cex=1.2) 
abline(v=0, lty=1) 
abline(v=10, lty=2, lwd=1, col="grey") 
abline(v=20, lty=2, lwd=1, col="black") 
legend("topright", legend = c("Unweighted", "Weighted"), col=c("blue", "red"), text.col=c("blue", "red"), pch=c(21,16))
}
dev.off()
pdf("tento18.pdf")
par(mfrow=c(3,3))
for(i in 10:18){
dotchart(peach1[i,], pch="", labels=nam2, cex=0.75) 
mtext("Standardized Difference (%)", side=1, line=2)
points(orderedpear3[i,], seq(1:length(orderedpear3[i,])), 
pch=21, col="blue", cex=1.2) 
points(peach1[i,],seq(1:length(peach1[i,])), 
pch=16, col="red", cex=1.2) 
abline(v=0, lty=1) 
abline(v=10, lty=2, lwd=1, col="grey") 
abline(v=20, lty=2, lwd=1, col="black") 
legend("topright", legend = c("Unweighted", "Weighted"), col=c("blue", "red"), text.col=c("blue", "red"), pch=c(21,16))
}
dev.off()
pdf("nineteento27.pdf")
par(mfrow=c(3,3))
for(i in 19:27){
dotchart(peach1[i,], pch="", labels=nam2, cex=0.75) 
mtext("Standardized Difference (%)", side=1, line=2)
points(orderedpear3[i,], seq(1:length(orderedpear3[i,])), 
pch=21, col="blue", cex=1.2) 
points(peach1[i,],seq(1:length(peach1[i,])), 
pch=16, col="red", cex=1.2) 
abline(v=0, lty=1) 
abline(v=10, lty=2, lwd=1, col="grey") 
abline(v=20, lty=2, lwd=1, col="black") 
legend("topright", legend = c("Unweighted", "Weighted"), col=c("blue", "red"), text.col=c("blue", "red"), pch=c(21,16))
}
dev.off()
pdf("twentyeightto36.pdf")
par(mfrow=c(3,3))
for(i in 28:36){
dotchart(peach1[i,], pch="", labels=nam2, cex=0.75) 
mtext("Standardized Difference (%)", side=1, line=2)
points(orderedpear3[i,], seq(1:length(orderedpear3[i,])), 
pch=21, col="blue", cex=1.2) 
points(peach1[i,],seq(1:length(peach1[i,])), 
pch=16, col="red", cex=1.2) 
abline(v=0, lty=1) 
abline(v=10, lty=2, lwd=1, col="grey") 
abline(v=20, lty=2, lwd=1, col="black") 
legend("topright", legend = c("Unweighted", "Weighted"), col=c("blue", "red"), text.col=c("blue", "red"), pch=c(21,16))
}
dev.off()

dotchart(peach1[1,], pch="", labels=nam2, cex=0.75) 
mtext("Standardized Difference (%)", side=1, line=2)
points(orderedpear3[1,], seq(1:length(orderedpear3[1,])), 
pch=21, col="blue", cex=1.2) 
points(peach1[1,],seq(1:length(peach1[1,])), 
pch=16, col="red", cex=1.2) 
abline(v=0, lty=1) 
abline(v=10, lty=2, lwd=1, col="grey") 
abline(v=20, lty=2, lwd=1, col="black") 
legend("topright", legend = c("Unweighted", "Weighted"), col=c("blue", "red"), text.col=c("blue", "red"), pch=c(21,16)) 


dat<-data.frame(SMD=orderedpear1, Variable=nam2, Sample=c(rep("Survey Sub-Sample", 3)))
dat2<-data.frame(SMD=peach1, Variable=nam2, Sample=c(rep("Survey Sample", 3)))
dat3<-data.frame(rbind(dat, dat2))
dat3$Variable<-factor(dat3$Variable, c("Y", "W2", "W1"))

grapha<-ggplot(dat3, aes(x=SMD, y=Variable, shape=Sample)) + geom_point(size=4)+  scale_shape_manual(values=c(1,2))+ theme_bw(base_size = 14, base_family = "Helvetica") + xlab("Standarized Mean Difference, %")  + geom_vline(aes(xintercept=0), linetype="dashed")+ theme_classic() +theme(axis.text = element_text(size=14), axis.title=element_text(size=14))+ theme(axis.text.y=element_text(face="italic")) +  scale_x_continuous(limits = c(-5,150)) + theme(legend.background = element_rect(fill = "white", color = "black",    linetype="solid"), legend.key = 
element_rect(fill = 'white', color = "white", size = 0.1), 
              legend.justification=c(1,1), legend.position=c(1,1), legend.text=element_text(size=14), legend.title = element_blank())


#Second, generate data from Scen2
#f<-svysub[,c('z1', 'z2', 't', 'y')]

ctxmeans<-daply(c, .(t), colwise(mean))
ctxvars<-daply(c, .(t), colwise(var))

MIEC_data_point7 <- MIEC(mainsample[,c("cmatAgeworse", T.outcomes, Z.covariates)], smallcalib[,c("cmatAge", "mmatAge")], n_calib,n_main,M=m,N=n,K=q,S=r)
MIEC_data_point4 <- MIEC(mainsample[,c("cmatAgeevenworse", T.outcomes, Z.covariates)], smallcalib[,c("cmatAge", "mmatAge")], n_calib,n_main,M=m,N=n,K=q,S=r)

T.outcomes<-c("tertscore", "internal")
MIEC_data_point9internal <- MIEC(mainsample[,c("cmatAge", T.outcomes, Z.covariates)], smallcalib[,c("cmatAge", "mmatAge")], n_calib,n_main,M=m,N=n,K=q,S=r)
MIEC_data_point7internal <- MIEC(mainsample[,c("cmatAgeworse", T.outcomes, Z.covariates)], smallcalib[,c("cmatAge", "mmatAge")], n_calib,n_main,M=m,N=n,K=q,S=r)
MIEC_data_point4internal <- MIEC(mainsample[,c("cmatAgeevenworse", T.outcomes, Z.covariates)], smallcalib[,c("cmatAge", "mmatAge")], n_calib,n_main,M=m,N=n,K=q,S=r)

point9<-multiple_imputation_EC_congenial(mainsample[,c("cmatAge", T.outcomes, Z.covariates)], smallcalib[,c("cmatAge", "mmatAge")], outcomevector=mainsample$substance)
point7<-multiple_imputation_EC_congenial(mainsample[,c("cmatAgeworse", T.outcomes, Z.covariates)], smallcalib[,c("cmatAgeworse", "mmatAge")], outcomevector=mainsample$substance)
point4<-multiple_imputation_EC_congenial(mainsample[,c("cmatAgeevenworse", T.outcomes, Z.covariates)], smallcalib[,c("cmatAgeevenworse", "mmatAge")], outcomevector=mainsample$substance)

T.outcomes<-c("tertscore", "internal")
point9internal<-multiple_imputation_EC_congenial(mainsample[,c("cmatAge", T.outcomes, Z.covariates)], smallcalib[,c("cmatAge", "mmatAge")], outcomevector=mainsample$internal)
point7internal<-multiple_imputation_EC_congenial(mainsample[,c("cmatAgeworse", T.outcomes, Z.covariates)], smallcalib[,c("cmatAgeworse", "mmatAge")], outcomevector=mainsample$internal)
point4internal<-multiple_imputation_EC_congenial(mainsample[,c("cmatAgeevenworse", T.outcomes, Z.covariates)], smallcalib[,c("cmatAgeevenworse", "mmatAge")], outcomevector=mainsample$internal)

#point9internalcong<-multiple_imputation_EC_congenial(mainsample[,c("cmatAge", T.outcomes, Z.covariates)], smallcalib[,c("cmatAge", "mmatAge")], outcomevector=mainsample$internal)
#point7internalcong<-multiple_imputation_EC_congenial(mainsample[,c("cmatAgeworse", T.outcomes, Z.covariates)], smallcalib[,c("cmatAgeworse", "mmatAge")], outcomevector=mainsample$internal)
#point4internalcong<-multiple_imputation_EC_congenial(mainsample[,c("cmatAgeevenworse", T.outcomes, Z.covariates)], smallcalib[,c("cmatAgeevenworse", "mmatAge")], outcomevector=mainsample$internal)

# Truth
#truemodel<- 'tertscore ~ poly(mmatAge,2) + SEXF + age_cent + urbancat + suburb + black + latino + other + midwest + south + west + cinc'
txmodel<- 'tertscore ~ mmatAge + SEXF + age_cent + moth + fath + black + latino + other + midwest + south + west + cinc'
outcomemodel<- 'substance~ tertscore'

truth<-IPTWest(txmodel, mainsample, outcomemodel)

outcomemodelinternal<- 'internal ~ tertscore'
truthinternal<-IPTWest(txmodel, mainsample, outcomemodelinternal)

# Naive
naivepoint9model<- 'tertscore ~ cmatAge + SEXF + age_cent + moth + fath + black + latino + other + midwest + south + west + cinc'
naivepoint7model<- 'tertscore ~ cmatAgeworse + SEXF + age_cent + moth + fath + black + latino + other + midwest + south + west + cinc'
naivepoint4model<- 'tertscore ~ cmatAgeevenworse + SEXF + age_cent + moth + fath  + black + latino + other + midwest + south + west + cinc'
naivepoint9<-IPTWest(naivepoint9model, mainsample, outcomemodel)
naivepoint7<-IPTWest(naivepoint7model, mainsample, outcomemodel)
naivepoint4<-IPTWest(naivepoint4model, mainsample, outcomemodel)

naivepoint9internal<-IPTWest(naivepoint9model, mainsample, outcomemodelinternal)
naivepoint7internal<-IPTWest(naivepoint7model, mainsample, outcomemodelinternal)
naivepoint4internal<-IPTWest(naivepoint4model, mainsample, outcomemodelinternal)

#MIEC without Y in the imputation model
T.outcomes<-c("tertscore")
q<-length(T.outcomes)
point9noy<-multiple_imputation_EC_noY(mainsample[,c("cmatAge", T.outcomes, Z.covariates)], smallcalib[,c("cmatAge", "mmatAge")], outcomevector=mainsample$substance)
point7noy<-multiple_imputation_EC_noY(mainsample[,c("cmatAgeworse", T.outcomes, Z.covariates)], smallcalib[,c("cmatAgeworse", "mmatAge")], outcomevector=mainsample$substance)
point4noy<-multiple_imputation_EC_noY(mainsample[,c("cmatAgeevenworse", T.outcomes, Z.covariates)], smallcalib[,c("cmatAgeevenworse", "mmatAge")], outcomevector=mainsample$substance)

point9noyinternal<-multiple_imputation_EC_noY(mainsample[,c("cmatAge", T.outcomes, Z.covariates)], smallcalib[,c("cmatAge", "mmatAge")], outcomevector=mainsample$internal)
point7noyinternal<-multiple_imputation_EC_noY(mainsample[,c("cmatAgeworse", T.outcomes, Z.covariates)], smallcalib[,c("cmatAgeworse", "mmatAge")], outcomevector=mainsample$internal)
point4noyinternal<-multiple_imputation_EC_noY(mainsample[,c("cmatAgeevenworse", T.outcomes, Z.covariates)], smallcalib[,c("cmatAgeevenworse", "mmatAge")], outcomevector=mainsample$internal)

#MIEC without T or Y in the imputation model
T.outcomes<-NULL
q<-length(T.outcomes)
point9noty<-multiple_imputation_EC_noTY(mainsample[,c("cmatAge", T.outcomes, Z.covariates)], smallcalib[,c("cmatAge", "mmatAge")], treatmentvector=mainsample$tertscore,  outcomevector=mainsample$substance)
point7noty<-multiple_imputation_EC_noTY(mainsample[,c("cmatAgeworse", T.outcomes, Z.covariates)], smallcalib[,c("cmatAgeworse", "mmatAge")], treatmentvector=mainsample$tertscore,outcomevector=mainsample$substance)
point4noty<-multiple_imputation_EC_noTY(mainsample[,c("cmatAgeevenworse", T.outcomes, Z.covariates)], smallcalib[,c("cmatAgeevenworse", "mmatAge")], treatmentvector=mainsample$tertscore,outcomevector=mainsample$substance)

point9notyinternal<-multiple_imputation_EC_noTY(mainsample[,c("cmatAge", T.outcomes, Z.covariates)], smallcalib[,c("cmatAge", "mmatAge")], treatmentvector=mainsample$tertscore,outcomevector=mainsample$internal)
point7notyinternal<-multiple_imputation_EC_noTY(mainsample[,c("cmatAgeworse", T.outcomes, Z.covariates)], smallcalib[,c("cmatAgeworse", "mmatAge")], treatmentvector=mainsample$tertscore, outcomevector=mainsample$internal)
point4notyinternal<-multiple_imputation_EC_noTY(mainsample[,c("cmatAgeevenworse", T.outcomes, Z.covariates)], smallcalib[,c("cmatAgeevenworse", "mmatAge")], treatmentvector=mainsample$tertscore, outcomevector=mainsample$internal)

#Regression Calibration
model_rc_point9<-lm(mmatAge ~ cmatAge, data=smallcalib)
model_rc_point7<-lm(mmatAge ~ cmatAgeworse, data=smallcalib)
model_rc_point4<-lm(mmatAge ~ cmatAgeevenworse, data=smallcalib)

mainsample$xpredpoint9<-predict(model_rc_point9, newdata=mainsample)
mainsample$xpredpoint7<-predict(model_rc_point7, newdata=mainsample)
mainsample$xpredpoint4<-predict(model_rc_point4, newdata=mainsample)

rc_point9_model<- 'tertscore ~ xpredpoint9 + SEXF + age_cent + moth + fath + black + latino + other + midwest + south + west + cinc'
rc_point7_model<- 'tertscore ~ xpredpoint7 + SEXF + age_cent +moth + fath + black + latino + other + midwest + south + west + cinc'
rc_point4_model<- 'tertscore ~ xpredpoint4 + SEXF + age_cent + moth + fath + black + latino + other + midwest + south + west + cinc'

rcpoint9<-IPTWest(rc_point9_model, mainsample, outcomemodel)
rcpoint7<-IPTWest(rc_point7_model, mainsample, outcomemodel)
rcpoint4<-IPTWest(rc_point4_model, mainsample, outcomemodel)

rcpoint9internal<-IPTWest(rc_point9_model, mainsample, outcomemodelinternal)
rcpoint7internal<-IPTWest(rc_point7_model, mainsample, outcomemodelinternal)
rcpoint4internal<-IPTWest(rc_point4_model, mainsample, outcomemodelinternal)

#changed not to include regression calibration 23 Jul 2014 
#changed to  not include y and t in the uncongeniel version
dfc<-data.frame(rbind(naivepoint9, naivepoint7, naivepoint4, point9noty, point7noty, point4noty, point9, point7, point4, truth, naivepoint9internal, naivepoint7internal, naivepoint4internal, point9notyinternal, point7notyinternal, point4notyinternal, point9internal, point7internal, point4internal, truthinternal) )
colnames(dfc)<-c("coef", "se", "CI_low", "CI_upp")
dfc$correlation<-rep(c(rep(c(cor(smallcalib$cmatAge, smallcalib$mmatAge), cor(smallcalib$cmatAgeworse, smallcalib$mmatAge), cor(smallcalib$cmatAgeevenworse, smallcalib$mmatAge)), 3),  1),2)
dfc$method<-rep(factor(c(rep("Naive (W)", 3), rep("Uncongenial MI-EC", 3),  rep("Congenial MI-EC", 3), "Truth (X)")),2)
dfc$outcome<-c(rep("substance", 10), rep("internal",10))
dfc$outcome<-factor(dfc$outcome, levels=c("substance", "internal"), labels=c("Substance Abuse/Dependence", "Anxiety or Depression"))
dfc$corrfact<-factor(dfc$correlation, levels=c(cor(smallcalib$cmatAgeevenworse, smallcalib$mmatAge),cor(smallcalib$cmatAgeworse, smallcalib$mmatAge), cor(smallcalib$cmatAge, smallcalib$mmatAge), 1), labels=c("0.30", ".72", "0.94", "true"))
dfc$method<-factor(dfc$method, levels=c("Truth (X)","Naive (W)", "Congenial MI-EC", "Uncongenial MI-EC"), labels=c("Truth (X)","Naive (W)", "Congenial MI-EC", "Uncongenial MI-EC"))

setwd("/Users/kararudolph/Documents/JHSPHpostdoc/measurementErrorImputation/MIEC")
#png("IllustExResnofact.png")
pd <- position_dodge(.1)
nofact<-ggplot(dfc, aes(x=correlation, y=coef, group=method, linetype=method)) + facet_wrap(~outcome) +
    geom_errorbar(width=.1, aes(ymin=CI_low, ymax=CI_upp), position=pd) +
    geom_point(position=pd) + ylab("ATT") + xlab("correlation between X and W") + xlim(.2, 1.1) + scale_linetype_manual(values=c("Uncongenial MI-EC"="longdash","Congenial MI-EC"= "solid", "Naive (W)"="dotdash", "Truth (X)"="dotted")) + scale_colour_grey()+ theme_bw()
#png("IllustExResnofactcol.png")
nofactcol<-ggplot(dfc, aes(x=correlation, y=coef, group=method, linetype=method, color=method)) + facet_wrap(~outcome) +
    geom_errorbar(width=.1, aes(ymin=CI_low, ymax=CI_upp), position=pd) +
    geom_point(position=pd) + ylab("ATT") + xlab("correlation between X and W") + xlim(.2, 1.1) + scale_linetype_manual(values=c("Uncongenial MI-EC"="longdash","Congenial MI-EC"= "solid", "Naive (W)"="dotdash", "Truth (X)"="dotted")) + theme_bw()

ggsave(nofact, file="IllustExResnofact.png")
ggsave(nofactcol, file="IllustExResnofactcol.png")

pd <- position_dodge(.5)
png("IllustExRes.png")
ggplot(dfc, aes(x=corrfact, y=coef, group=method, linetype=method)) + facet_wrap(~outcome) +
    geom_errorbar(width=.1, aes(ymin=CI_low, ymax=CI_upp), position=pd) +
    geom_point(position=pd) + ylab("ATT") + xlab("correlation between X and W") + theme_bw()
dev.off()