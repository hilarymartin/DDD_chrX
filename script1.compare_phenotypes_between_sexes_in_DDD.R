setwd("/lustre/scratch115/projects/ddd/users/hcm/DDD/chrX_analyses/code_and_data_to_release")
library(epiDisplay)
#define a set of quantitative traits to look at - age at assessment, milestones, anthropometric measurements, number of organ systems affected
quant.phenos=c("decimal_age_at_assessment","talking.cleaned.reached","walking.cleaned.reached","birthweight_sd","birth_ofc_sd","height_sd","weight_sd","ofc_sd","n.nonredundant.organ.systems")

#read the phenotype data
# The DDD phenotype data is only accessible via EGA under managed access. Hence, I am not sharing the phenotype file with the code. To see the format of this phenotype file, see fake_phenotype_data.example_file_to_illustrate_format.txt as an example.
pheno=read.delim("NOT_FOR_RELEASE/phenotype_data.NOT_FOR_RELEASE.txt",header=T,as.is=T)

# for dichotomous traits, let's just consider those that affect at least 1% of individuals in the cohort
traits=colnames(pheno)[2:48]

# Test for a difference in the frequency of different phenotypic abnormalities between sexes, controlling for age at assessment
logreg.results=NULL
for(t in traits){
  pheno$t=pheno[,t]
  my.logreg= glm(pheno$t~pheno$gender + pheno$decimal_age_at_assessment,family=binomial(link='logit'))
  percentage.female=100*sum(pheno$t[pheno$gender=="Female"])/sum(pheno$gender=="Female")
  percentage.male=100*sum(pheno$t[pheno$gender=="Male"])/sum(pheno$gender=="Male")
  logreg.results=rbind(logreg.results,c(percentage.female,percentage.male,summary(my.logreg)$coefficients[2,c(1,4)]))
}
logreg.results=as.data.frame(logreg.results)
colnames(logreg.results)[1:2]=c("percentage.of.females","percentage.of.males")
rownames(logreg.results)=traits
logreg.results$OR=exp(logreg.results$Estimate)
write.table(logreg.results,"Supp_Table_1.results_from_logistic_regression_comparing_presence_absence_of_phenotypes_between_sexes_controlling_for_age_at_assessment.txt",quote=F,sep="\t")

# plot distribution of number of affected organ systems between sexes
pdf("Supp_Figure_1.comparing_number_of_nonredundant_affected_organ_systems_between_sexes.pdf",height=5,width=5)
barplot(prop.table(table(pheno$n.nonredundant.organ.systems[pheno$n.nonredundant.organ.systems>0 & pheno$gender=="Male"])),
        xlab="number of affected organ systems",border="blue",main="",col=NA)
barplot(prop.table(table(pheno$n.nonredundant.organ.systems[pheno$n.nonredundant.organ.systems>0 & pheno$gender=="Female"])),add=T,border="red",col=NA)
legend("topright",c("males","females"),col=c("blue","red"),lty=1)
dev.off()

#### Now let's consider the quantitative phenotypes - anthropometric measurements and growth milestones
# write out a table of the number of samples with missing data on each phenotyp
write.table(t(apply(pheno[,quant.phenos],2,function(x){table(is.na(x))})),"number_of_probands_excluded_from_linear_regressions_on_quantitative_phenotypes_due_to_missing_data_or_not_yet_reached.txt",quote=F,sep="\t")

# Compare quantitative phenotypes between sexes
lm.results=NULL
for(q in quant.phenos){
my.lm=summary(lm(pheno[,q]~pheno$gender))
lm.results=rbind(lm.results,my.lm$coefficients[2,])
}
rownames(lm.results)=quant.phenos
write.table(lm.results,"results_of_comparing_milestones_and_growth_metrics_between_sexes.milestones_on_for_those_achieved_already.txt",quote=F,sep="\t")

# Compare quantitative phenotypes between sexes, controlling for age at assesment
lm.results2=NULL
for(q in quant.phenos[2:length(quant.phenos)]){
my.lm=summary(lm(pheno[,q]~pheno$gender + pheno$decimal_age_at_assessment))
lm.results2=rbind(lm.results2,my.lm$coefficients[2,])
}
rownames(lm.results2)=quant.phenos[2:length(quant.phenos)]
write.table(lm.results2,"Supp_Table_2.results_of_comparing_milestones_and_growth_metrics_between_sexes.milestones_on_for_those_achieved_already.controlling_for_age_at_assessment.txt",quote=F,sep="\t")




