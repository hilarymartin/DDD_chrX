setwd("/lustre/scratch115/projects/ddd/users/hcm/DDD/chrX_analyses/code_and_data_to_release/")
#read in the PRS data with principal components
master.coreex=read.delim("PRS_data_with_annotations.txt",header=T,as.is=T)
#define the PRSs we will test
prs.cols=c( "ddd_prs.p_1.r2_0.1.8_jan_2018.SD","iq_prs.p_0.05.r2_0.1.26_dec_2017.SD","clozuk_prs.p_0.05.r2_0.1.8_jan_2018.SD","edu_2018_prs.ssgaccorrected.p_1.r2_0.1.19_june_2018.SD")

x="male.suspected.xlinked.excl.xlinked.diag" #this column is a flag indicating whether the patient was suspected to have an X-linked disorder due to family history, excluding those with known X-linked diagnoses.
# We are analysing a set of male probands without an X-linked molecular diagnosis from DDD who have genotype data on the CoreExome chip, from Niemi et al., Nature, 2018
# We do a linear regression of each PRSs on whether or not the patient was suspected by a clinician to have an X-linked disorder due to the family history.
results=NULL
prs.saved=c()
for(prs in prs.cols){
    prs.saved=c(prs.saved,prs)
    master.coreex$x=master.coreex[,x]
    master.coreex$y=master.coreex[,prs]
    my.lm=summary(lm(y~x+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10,data=master.coreex))
    results=rbind(results,c(my.lm$coefficients[2,]))
  }

results=as.data.frame(results)
results$PRS=prs.saved

write.table(results,"Supp_Table_6.PRS_comparisons.txt",quote=F,sep="\t",row.names=F)



