setwd("/lustre/scratch115/realdata/mdt0/projects/ddd/users/hcm/DDD/chrX_analyses/code_and_data_to_release")
#### Note that this script takes a while to run. You can shorted it by changing the number of options for pi (the TADA parameter for fraction of true disease genes) that it loops over in this line: for(pi in c(0.01,0.05,0.1,0.15,0.2,0.25)){

# run the TADA code to load functions into this script
source("TADA.v.1.2.R")

# load data for TADA
load("input_for_TADA.males.case_control_and_de_novos.RData")

#specify the sample sizes in the analysis.
n.family = 5138 #this is the number of trios
n.case = 7136 #this is the number of male cases
n.ctrl = 8551 #this is the number of fathers

### convert the mutation rates back to a rate per individual
tada.data$mut.cls1=tada.data$mut.cls1/n.family
tada.data$mut.cls2=tada.data$mut.cls2/n.family

### remove any genes which don't have a mutation rate
tada.data=tada.data[!is.na(tada.data$mut.cls1),]

### set a list of sample sizes
#n = data.frame(dn=n.family, ca=n.case, cn=n.ctrl)
#need to set the sample size for the de novo analysis to # trios/2, because TADA is assuming samples are diploid, whereas here we're dealing with haploid males
n = data.frame(dn=n.family/2, ca=n.case, cn=n.ctrl)
sample.counts <- list(cls1=n, cls2=n)
sample.counts$cls1["dn"]=round(sample.counts$cls1["dn"])
sample.counts$cls2["dn"]=round(sample.counts$cls2["dn"])

# set up the count data so that it can be used by the main TADA function; cls1=LoF; cls2=functional (missense/inframe)
# The data of one mutational class is stored in a data frame with three columns for de novo, case and control counts
cls1.counts=data.frame(dn=tada.data$dn.cls1, ca=tada.data$case.cls1, cn=tada.data$ctrl.cls1)
cls2.counts=data.frame(dn=tada.data$dn.cls2, ca=tada.data$case.cls2, cn=tada.data$ctrl.cls2)
rownames(cls1.counts)=tada.data$gene.id
rownames(cls2.counts)=tada.data$gene.id
tada.counts=list(cls1=cls1.counts,cls2=cls2.counts)

#set mutation rates
mu=data.frame(cls1=tada.data$mut.cls1,cls2=tada.data$mut.cls2)

#We first set up TADA in the full, instead of denovo, mode. Note that it is possible to run TADA in the de novo mode for some classes, but the full mode for other ones.
denovo.only=data.frame(cls1=FALSE,cls2=FALSE)

### function to calculate the relative risk given the burden of a particular class of variant, and a prior assumption about the fraction of risk genes
relative.risk=function(burden,fraction.risk.genes){
  return(1 + ((burden - 1)/fraction.risk.genes))
}

# define a set of known X-linked disease genes based on DDG2P
ddg2p=read.csv('DDG2P_25_1_2019.csv',header=T,as.is=T)
hemi.genes=unique(ddg2p[ddg2p$DDD.category %in% c("confirmed","probable") & ddg2p$allelic.requirement %in% "hemizygous","gene.symbol"])
xld.genes=unique(ddg2p[ddg2p$DDD.category %in% c("confirmed","probable") & ddg2p$allelic.requirement %in% "x-linked dominant","gene.symbol"])
xlod.genes=unique(ddg2p[ddg2p$DDD.category %in% c("confirmed","probable") & ddg2p$allelic.requirement %in% "x-linked over-dominance","gene.symbol"])
disease.genes=unique(c(hemi.genes,xld.genes,xlod.genes))

# set up a list for the results
overall.results=list()
# set up a matrix to save the parameters used for TADA with different values of pi (the fraction of causal genes)
saved.params=NULL
for(pi in c(0.05)){
#for(pi in c(0.01,0.05,0.1,0.15,0.2,0.25)){
cat('assuming fraction of risk genes = ',pi,'\n')

#We start by setting the parameters for TADA, following the instructions in http://www.compgen.pitt.edu/TADA/TADA_guide.html
rr.lof.dnm=relative.risk(sum(tada.data$dn.cls1)/sum(tada.data$mut.cls1 * n.family),pi)
rr.func.dnm=relative.risk(sum(tada.data$dn.cls2)/sum(tada.data$mut.cls2 * n.family),pi)

#parameters for LoFs
mutation="cls1"
counts <- tada.counts[[mutation]]
N <- sample.counts[[mutation]]
q.ca <- sum(counts[,"ca"])/N[,"ca"]
q.cn <- sum(counts[,"cn"])/N[,"cn"]
rr.lof.cc=relative.risk(burden = q.ca/q.cn,fraction.risk.genes = pi)
if(rr.lof.dnm>4){
  beta.dn.lof=1
}else {
  beta.dn.lof = 10/rr.lof.dnm
}
if(rr.lof.cc>4){
  beta.cc.lof=1
}else {
  beta.cc.lof = 10/rr.lof.dnm
}
counts <- tada.counts[[mutation]]
q.mean.lof <- (sum(counts[rownames(counts) %in% disease.genes,"ca"]) + sum(counts[rownames(counts) %in% disease.genes,"cn"])) / (N[,"ca"] + N[,"cn"])

nu <- 100
rho.lof <- nu * q.mean.lof
rho1.lof <- rho.lof
nu1 <- nu
rho0.lof <- rho.lof
nu0 <- nu

#parameters for missense/inframe variants
mutation="cls2"
counts <- tada.counts[[mutation]]
N <- sample.counts[[mutation]]
q.ca <- sum(counts[,"ca"])/N[,"ca"]
q.cn <- sum(counts[,"cn"])/N[,"cn"]
rr.func.cc=relative.risk(burden = q.ca/q.cn,fraction.risk.genes = pi)
if(rr.func.dnm>4){
  beta.dn.func=1
}else {
  beta.dn.func = 10/rr.func.dnm
}
if(rr.func.cc>4){
  beta.cc.func=1
}else {
  beta.cc.func = 10/rr.func.cc
}
counts <- tada.counts[[mutation]]
q.mean.func <- (sum(counts[rownames(counts) %in% disease.genes,"ca"]) + sum(counts[rownames(counts) %in% disease.genes,"cn"])) / (N[,"ca"] + N[,"cn"])

nu <- 100
rho.func <- nu * q.mean.func
rho1.func <- rho.func
nu1 <- nu
rho0.func <- rho.func
nu0 <- nu

#set up the input parameters, for TADA to use
cls1= data.frame(gamma.mean.dn=rr.lof.dnm, beta.dn=beta.dn.lof, gamma.mean.CC=rr.lof.cc, beta.CC=beta.cc.lof, rho1=rho.lof, nu1=nu1, rho0=rho.lof, nu0=nu1)
cls2= data.frame(gamma.mean.dn=rr.func.dnm, beta.dn=beta.dn.func, gamma.mean.CC=rr.func.cc, beta.CC=beta.cc.func, rho1=rho.func, nu1=nu1, rho0=rho.func, nu0=nu1)
hyperpar=list(cls1=cls1,cls2=cls2)

#save the input parameters for this value of pi in a data frame
saved.params=rbind(saved.params,c(rr.lof.dnm,beta.dn.lof,rr.lof.cc,beta.cc.lof,rho.lof,rr.func.dnm,beta.dn.func,rr.func.cc,beta.cc.func,rho.func))

# call TADA, calculate Bayes factors for all genes
re.TADA <- do.call(cbind.data.frame, TADA(tada.counts=tada.counts, sample.counts=sample.counts, mu=mu, hyperpar=hyperpar, denovo.only=denovo.only))
# calculate q-values for all genes
re.TADA$qval.lof=Bayesian.FDR(re.TADA$BF.cls1, pi0 = 1-pi)
re.TADA$qval=Bayesian.FDR(re.TADA$BF.total, pi0 = 1-pi)
#  obtain p-values using a sampling approach to obtain the null distribution with function TADAnull
re.TADA.null=do.call(cbind.data.frame, TADAnull(tada.counts=tada.counts, sample.counts=sample.counts, mu=mu, hyperpar=hyperpar, denovo.only=denovo.only, nrep=5000))
re.TADA$pval.lof=bayesFactor.pvalue(re.TADA$BF.cls1,re.TADA.null$BFnull.cls1)
re.TADA$pval=bayesFactor.pvalue(re.TADA$BF.total,re.TADA.null$BFnull.total)

#sort by p-value
re.TADA=re.TADA[order(re.TADA$pval),]

#combine with the input data
combined=cbind(tada.data,re.TADA[rownames(tada.data),])

#save input data and results in a list
overall.results[[paste("frac_risk_genes_",pi,sep="")]] = combined

}

colnames(saved.params) =c("rr.dnm.lof","beta.dnm.lof","rr.cc.lof","beta.cc.lof","rho.lof","rr.dnm.func","beta.dnm.func","rr.func.cc","beta.cc.func","rho.func")
saved.params=as.data.frame(saved.params)
saved.params$frac.risk.genes=as.numeric(gsub("frac_risk_genes_","",names(overall.results)))
#write.table(saved.params,"Supp_Table_7.parameters_run_with_TADA.txt",quote=F,sep="\t",row.names=T)

for(i in 1:length(overall.results)){
  write.table(overall.results[[i]],paste("output_from_TADA.",names(overall.results)[i],".txt",sep=""),row.names=F,quote=F,sep="\t")
}
