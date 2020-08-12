setwd("/lustre/scratch115/projects/ddd/users/hcm/DDD/chrX_analyses/code_and_data_to_release")
library(data.table)
# load output of burden analysis on DNMs and from case/control analysis
burden.dnms=read.delim("burden_analysis_of_de_novos.different_chrX_gene_sets.txt",header=T,as.is=T)
burden.case.control=read.delim("burden_analysis_on_males.from_case_control_analysis.txt",header=T,as.is=T)

#load sets of X-linked genes
load("gene_sets_to_use_in_burden_analysis.DDG2P_14_1_2020.consensus_with_OMIM.RData")

### Let's start by calculating some key statistics that are mentioned in the text
# What proportion of the excess is in known genes in females?
burden.dnms.f=burden.dnms[burden.dnms$geneset %in% c("XL","all") & burden.dnms$csq %in% c("LoF.and.functional") & burden.dnms$sex %in% c("F"),]
prop.known.females=burden.dnms.f[burden.dnms.f$geneset=="XL","excess.over.N"]/burden.dnms.f[burden.dnms.f$geneset=="all","excess.over.N"]
prop.known.females
#[1] 0.9012927

# What proportion of the excess is in known genes in males?
burden.all.m=burden.case.control[c("LoF.and.functional.0_hemi_gnomAD_all","LoF.and.functional.0_hemi_gnomAD_XL"),]
prop.known.males=burden.all.m["LoF.and.functional.0_hemi_gnomAD_XL","excess.vars.over.N"]/burden.all.m["LoF.and.functional.0_hemi_gnomAD_all","excess.vars.over.N"]
prop.known.males
#[1] 0.6372547

#  What proportion of the excess is in known genes in both sexes combined? i.e. including DNMs for females and DNMs+inherited variants for males
proportion.known.genes=(burden.all.m["LoF.and.functional.0_hemi_gnomAD_XL","excess.vars.over.N"]+burden.dnms.f[burden.dnms.f$geneset=="XL","excess.over.N"])/
  (burden.all.m["LoF.and.functional.0_hemi_gnomAD_all","excess.vars.over.N"]+burden.dnms.f[burden.dnms.f$geneset=="all","excess.over.N"])
proportion.known.genes
# [1] 0.7786889

# What proportion of the cohort can be explained by X-linked causes not in known genes?
(1-proportion.known.genes) * (burden.all.m["LoF.and.functional.0_hemi_gnomAD_all","excess.vars.over.N"]*burden.all.m["LoF.and.functional.0_hemi_gnomAD_all","n.probands"]+
                                burden.dnms.f[burden.dnms.f$geneset=="all","excess.over.N"]*burden.dnms.f[burden.dnms.f$geneset=="all","N"])/
  (burden.all.m["LoF.and.functional.0_hemi_gnomAD_all","n.probands"]+burden.dnms.f[burden.dnms.f$geneset=="all","N"])
#[1] 0.01416263

### In males, what proportion is de novo versus inherited?
burden.dnms.m=burden.dnms[burden.dnms$geneset %in% c("XL","all") & burden.dnms$csq %in% c("LoF.and.functional") & burden.dnms$sex %in% c("M"),]
prop.dnm.in.males=burden.dnms.m[burden.dnms.m$geneset=="all","excess"]/burden.case.control["trios_only.LoF.and.functional.0_hemi_gnomAD_all","excess.vars"]
prop.dnm.in.males
#[1] 0.4240502

# Testing Haldane's theory: In trio males, what proportion of the excess in X-linked recessive genes is de novo versus inherited?  
# We focus on SNVs since the assumption that μ_male=3.5μ_female pertains to SNVs 
obs.dnms=burden.dnms[burden.dnms$csq=="LoF.and.missense.SNV" & burden.dnms$sex=="M" & burden.dnms$geneset=="Hemi.only","excess"]
obs.total=burden.case.control["trios_only.LoF.and.missense.SNV.0_hemi_gnomAD_Hemi.only","excess.vars"]
# approximate estimate of the proportion of the excess in X-linked recessive genes is de novo versus inherited
prop.dnm.in.males.XLR=obs.dnms/obs.total
prop.dnm.in.males.XLR
#[1] 0.3857306

## testing Haldane's equation - fraction de novo in X-linked recessive genes versus expectation
expectation.under.Haldane=1/5.5
prop.test(round(obs.dnms),round(obs.total),p=expectation.under.Haldane)
#X-squared = 19.384, df = 1, p-value = 1.069e-05
#sample estimates:
#  p 
#0.3888889

### Under Sherman's theory, and assuming m=1 and v=3.5mu, where m=reproductive loss in males, mu=mutation rate in egs, v=mutation rate in sperm,
## the expected proportion de novo is given by (f+1)/(5.5-3.5f), where f is reproductive loss in females
## This implies that f=(5.5p-1)/(3.5p+1)
## In UKBB, we estimate that # kids in carrier women is 1.31 and in non-carrier women is 1.76. This gives a reproductive loss in females of:
f=1-1.31/1.76
f
#[1] 0.2556818
## This gives an expected proportion de novo of
(f+1)/(5.5-3.5*f)
#[1] 0.2726712

## If we consider the lower bound of the fertility ratio from UKBB, we get the following expected proportion of de novos:
f=1-0.503
(f+1)/(5.5-3.5*f)
#[1] 0.3980854

#### Compare attributable fraction between males who were versus were not suspected to have X-linked inheritance by clincians, based on family history
attr.frac.xlinked=burden.case.control[c("probands_suspected_Xlinked.LoF.and.functional.0_hemi_gnomAD_all","LoF.and.functional.0_hemi_gnomAD_all"),]
attr.frac.xlinked[,c("excess.vars.over.N","att.frac.lt","att.frac.ut")]

#### calculate confidence intervals on these key statistics from bootstrapping
#load bootstrap results
load("results_of_bootstrapping_for_DNM_burden.RData")
load("results_of_bootstrapping_in_case_control_analysis.all.RData")

### proportion of the excess that is in known DD-associated genes in males
frac.known.m=boot.counts.list.new[["trios.only_FALSE_LoF.and.functional"]][,"fraction.excess.known"]

### proportion of the excess that is in known DD-associated genes in females
frac.known.f=bootstraps.by.sex[["F"]]$fraction.burden.XL

#### fraction of excess in males that is de novo (LoF+functional)
fraction.denovo.all=bootstraps.by.sex[["M"]]$excess.LoF.and.functional.all/boot.counts.list.new[["trios.only_TRUE_LoF.and.functional" ]][,"excess.all"]

### fraction of the excess in males that is de novo, focusing on SNVs in XLR genes
fraction.denovo.Hemi.only=bootstraps.by.sex[["M"]]$excess.LoF.and.missense.SNV.Hemi.only/boot.counts.list.new[[ "trios.only_TRUE_LoF.and.missense.SNV"]][,"excess.Hemi.only"]

### Summarise the observed values and 95% confidence intervals determined by bootstrapping
summary.metrics=cbind(c(prop.known.males,quantile(frac.known.m,c(0.025,0.975))),
c(prop.known.females,quantile(frac.known.f,c(0.025,0.975))),
c(prop.dnm.in.males,quantile(fraction.denovo.all,c(0.025,0.975))),
c(prop.dnm.in.males.XLR,quantile(fraction.denovo.Hemi.only,c(0.025,0.975))))
rownames(summary.metrics)=c("observed","2.5th_percentile_from_boostrapping","97.4th_percentile_from_boostrapping")
colnames(summary.metrics)=c("fraction.in.known.genes.males","fraction.in.known.genes.females","fraction.denovo.in.males.all.genes","fraction.SNVs.denovo.in.males.XLR.genes")
write.table(summary.metrics,"confidence_intervals_from_bootstrapping.txt",quote=F,sep="\t")


####### Figure 1A - Fraction of males and females attributable to rare inherited and de novo coding variants on the X chromosome
cols2=c("#98D7E5","#D37986","#E577E2")#)"#ffffbf","#fc8d59" ,"#d7191c")
names(cols2)=c("lof","functional","LoF.and.functional")
# merge results from the case/control analysis and the de novo analysis
burden.all.cc=burden.case.control[c( "synonymous.0_hemi_gnomAD_all" ,"functional.0_hemi_gnomAD_all" ,"LoF.0_hemi_gnomAD_all","LoF.and.functional.0_hemi_gnomAD_all" ),
                                  c("class","ratio.rate.vars","burden.lt","burden.ut","excess.vars.over.N","att.frac.lt","att.frac.ut","pois.p.ut","ppv.lt","ppv.ut")]
burden.all.cc$sex="M"
burden.all.cc$type="case-control"
burden.all.cc$csq=c("syn","functional","lof","LoF.and.functional")
colnames(burden.all.cc)=c("class","ratio","burden.lt","burden.ut","excess.over.N","att.frac.lt","att.frac.ut","ppois.ut","ppv.lt","ppv.ut","sex","type","csq")

burden.all.dnms=burden.dnms[burden.dnms$geneset=="all" & burden.dnms$csq %in% c("syn","lof","functional","LoF.and.functional") & burden.dnms$sex %in% c("M","F"),
  c("sex","csq","ratio","burden.lt","burden.ut","excess.over.N","att.frac.lt","att.frac.ut","ppois.ut")]
burden.all.dnms$type="de novo"

attr.frac.combined=rbind(burden.all.dnms[burden.all.dnms$csq %in% c("lof","functional","LoF.and.functional"),],burden.all.cc[burden.all.cc$csq %in% c("functional","lof","LoF.and.functional"),colnames(burden.all.dnms)])
attr.frac.combined=attr.frac.combined[order(attr.frac.combined$type,attr.frac.combined$sex,attr.frac.combined$csq),]
attr.frac.combined=attr.frac.combined[c(1:3,7:9,4:6),]
att.fraction=attr.frac.combined
att.fraction$cols=cols2[att.fraction$csq]
# now plot the attributable fraction
pdf("Figure_1A.attributable_fraction.indicating_male_DNM.pdf",height=5,width=6,useDingbats=F)
att.fraction1=att.fraction[!(att.fraction$sex=="M" & att.fraction$type=="de novo"),]
att.fraction2=att.fraction[(att.fraction$sex=="M" & att.fraction$type=="de novo")|(att.fraction$sex=="F"),]

x=barplot(att.fraction1$excess.over.N,col=att.fraction1$cols,names.arg=att.fraction1$sex,legend.text=c("missense/inframe","PTV","PTV+missense/inframe"),args.legend=list(x="topleft",bty="n"),ylab="attributable fraction",
  ylim=c(0,max(att.fraction$att.frac.ut)))
for(i in 1:nrow(att.fraction1)){
  segments(x0=x[i],x1=x[i],y0=att.fraction1[i,"att.frac.lt"],y1=att.fraction1[i,"att.frac.ut"],lwd=2)
}

x2=barplot(att.fraction2$excess.over.N,add=T,angle=45,density=15,col="black")
for(i in 1:nrow(att.fraction2)){
  segments(x0=x2[i],x1=x2[i],y0=att.fraction2[i,"att.frac.lt"],y1=att.fraction2[i,"att.frac.ut"],lwd=1.5)
}
dev.off()


###### Figure 1B: Relative fraction of PTV versus missense/inframe variants amongst ClinVar likely pathogenic or pathogenic variants in X-linked DDG2P genes, versus the fraction inferred in the burden analysis in DDD.
# load DDG2P genes
ddg2p=read.csv("DDG2P_25_1_2019.csv",header=T,as.is=T)

# load in Ensembl transcripts for DDG2P genes
ddg2p.ens=read.delim("DDG2P_25_1_2019.Ensembl_transcript_IDs.Ensembl_V100.txt",header=T,as.is=T)

# define transcripts corresponding to X-linked DD-associated genes
xlinked=ddg2p.ens[gsub("HGNC:","",ddg2p.ens$HGNC.ID) %in% unique(ddg2p[ddg2p$DDD.category %in% c("confirmed","probable") & ddg2p$allelic.requirement %in% c("x-linked dominant","x-linked over-dominance","hemizygous"),
  "hgnc.id"]),]

# load ClinVar variants
annot=as.data.frame(fread("clinvar_annotation.chrX.txt",header=T,stringsAsFactors=F))
x.info=as.data.frame(fread("clinvar.chrX.txt",header=T,stringsAsFactors=F))

# add on annotations to X-linked ClinVar variants
x.info=merge(x.info,annot,by.x="clinvar_code",by.y="clinvar_code",all.x=T)

#restrict to Clinvar Variants in X-linked DDG2P genes
#x.info=x.info[x.info$ensembl_transcript %in% xlinked$Ensembl.Transcript.ID,]
x.info=x.info[x.info$ensembl_transcript %in% xlinked$Transcript.stable.ID,]

# restrict to pathogenic or likely pathogenic ClinVar variants
x.info.path=x.info[x.info$pathogenicity %in% c("Pathogenic","Pathogenic, Likely pathogenic","Likely pathogenic"),]

# split pathogenic/likely pathogenic X-linked ClinVar variants into PTVs versus missense/inframe
lof.csq=c("splice_acceptor_variant","splice_donor_variant","stop_gained","frameshift_variant")
miss.csq=c("missense_variant","inframe_insertion","inframe_deletion")
x.info.path.summary=unique(x.info.path[,c("clinvar_code","vep_consequence")])
x.path.lof=x.info.path.summary[x.info.path.summary$vep_consequence %in% lof.csq,]
x.path.miss=x.info.path.summary[x.info.path.summary$vep_consequence %in% miss.csq,]
# remove variants in the missense/inframe group that are also annotated as PTV in another transcript
x.path.miss=x.path.miss[!x.path.miss[,1] %in% x.path.lof[,1],]

# Count pathogenic/likely pathogenic X-linked ClinVar variants that are missense/inframe versus PTV
clinvar.estimates=c( length(unique(x.path.miss[,1])), length(unique(x.path.lof[,1])))
names(clinvar.estimates)=c("missense/inframe","ptv")

# Determine the approximate number of males and females in DDD explained by LoF or missense/inframe variants overall in X-linked DD-associated genes, using the excess from the burden analysis
lof.cont=burden.dnms[burden.dnms$geneset %in% c("XL") & burden.dnms$csq %in% c("lof") & burden.dnms$sex %in% c("F"),"excess"]+burden.case.control[c( "LoF.0_hemi_gnomAD_XL"),"excess.vars"]
func.cont=burden.dnms[burden.dnms$geneset %in% c("XL") & burden.dnms$csq %in% c("functional") & burden.dnms$sex %in% c("F"),"excess"]+burden.case.control[c( "functional.0_hemi_gnomAD_XL"),"excess.vars"]
ddd.estimates=c(func.cont,lof.cont)
names(ddd.estimates)=c("missense/inframe","ptv")

# Compare the contribution of PTV versus missense/inframe variants between the ClinVar variants and the DDD burden analysis
prop.table(ddd.estimates)
prop.table(clinvar.estimates)

# test for a significant difference
compare.clinvar.to.ddd=fisher.test(cbind(ddd.estimates,clinvar.estimates))
compare.clinvar.to.ddd
compare.clinvar.to.ddd$p.value
#[1] 5.620052e-21

# Estimate the fraction of truly pathogenic missense/inframe variants that are not being classed as pathogenic
1-(clinvar.estimates[1]/clinvar.estimates[2])/(ddd.estimates[1]/ddd.estimates[2])
#       0.5866561 

pdf("Figure_1B.proportion_of_functional_versus_LoF_variants_in_ClinVar_versus_DDD_burden.pdf",height=5,width=4)
barplot(cbind(prop.table(clinvar.estimates),prop.table(ddd.estimates)),beside=F,names.arg=c("ClinVar pathogenic","DDD burden"),legend.text=c("missense/inframe","PTV"),col=c("#D37986","#98D7E5"),ylab="proportion of variants")
dev.off()

###### Figure 1C: Estimated attributable fraction versus positive predictive value for DNMs and inherited variants in males in X-linked DDG2P genes.

mynames=c("PTV","missense","PTV","missense","missense MPC<1","missense 1<MPC<2","missense 1<MPC<2, CADD>25","missense MPC>2","missense MPC>2, CADD>25")
cols=c("#98D7E5","#D37986","#98D7E5","#D37986","#FFF7EC","#FDD49E","#FDBB84","#EF6F3C","#EA3A1C","#d7301f")

burden.dnms$ppv=burden.dnms$excess/burden.dnms$obs
dnm.ppv=burden.dnms[burden.dnms$geneset=="XL" & burden.dnms$csq %in% c("lof","mis","functional"),
  c("sex","excess","obs","csq","exp","ppv.lt","ppv.ut",   "excess.over.N"  ,"att.frac.lt", "att.frac.ut"  )]
dnm.ppv$ppv=as.numeric(dnm.ppv$excess)/as.numeric(dnm.ppv$obs)
dnm.ppv=dnm.ppv[dnm.ppv$sex=="combined" &dnm.ppv$csq %in% c('lof','mis'),]

dnm.ppv=dnm.ppv[order(dnm.ppv$csq),]

inher.ppv=burden.case.control[c("inherited_only.trios_only.LoF.0_hemi_gnomAD_XL",
  "inherited_only.trios_only.missense.0_hemi_gnomAD_XL","inherited_only.trios_only.missense.MPC.lt.1.0_hemi_gnomAD_XL",
  "inherited_only.trios_only.missense.MPC.1.to.2.0_hemi_gnomAD_XL",
  "inherited_only.trios_only.missense.MPC.1.to.2.CADDgt25.0_hemi_gnomAD_XL",
  "inherited_only.trios_only.missense.MPC.gt.2.0_hemi_gnomAD_XL","inherited_only.trios_only.missense.MPC.gt.2.CADDgt25.0_hemi_gnomAD_XL"  ),]

pdf("Figure_1C.ppv_versus_attributable_fraction.pdf",height=6,width=8,useDingbats=F)
plot(c(dnm.ppv$ppv,inher.ppv$prop.vars.pathogenic),c(dnm.ppv$excess.over.N,inher.ppv$excess.vars.over.N),bg=cols,xlab="positive predictive value",ylab="attributable fraction",pch=c(rep(21,2),rep(23,7)),cex=1,col="black")
cis=rbind(dnm.ppv[,c("ppv.lt","ppv.ut")],inher.ppv[,c("ppv.lt","ppv.ut")])
x=c(dnm.ppv$ppv,inher.ppv$prop.vars.pathogenic)
y=c(dnm.ppv$excess.over.N,inher.ppv$excess.vars.over.N)
for(i in 1:nrow(cis)){
  segments(y0=y[i],y1=y[i],x0=cis[i,1],x1=cis[i,2],lwd=1.5,col="grey")
}
cis=rbind(dnm.ppv[,c("att.frac.lt","att.frac.ut")],inher.ppv[,c("att.frac.lt","att.frac.ut")])
for(i in 1:nrow(cis)){
  segments(x0=x[i],x1=x[i],y0=cis[i,1],y1=cis[i,2],lwd=1.5,col="grey")
}

points(c(dnm.ppv$ppv,inher.ppv$prop.vars.pathogenic),c(dnm.ppv$excess.over.N,inher.ppv$excess.vars.over.N),bg=cols,pch=c(rep(21,2),rep(23,7)),cex=1.5,col="black")
i=1
cis=rbind(dnm.ppv[,c("ppv.lt","ppv.ut")],inher.ppv[,c("ppv.lt","ppv.ut")])

legend.text=c("DNM","inherited","PTV","missense","missense MPC<1","missense 1<MPC<2","missense 1<MPC<2 & CADD>25","missense MPC>2","missense MPC>2 & CADD>25")
mycols=cols[3:length(cols)]
legend('top',legend.text,col=rep("black",9),pch=c(21,23,rep(23,7)),pt.bg=c("white","white",mycols))
dev.off()


#### Figure 2A: bivariate plot of de novo burden in males versus females, for the genome-wide significant genes
tada.genes=read.delim("output_from_TADA.frac_risk_genes_0.05.used_in_paper.txt",header=T,as.is=T)
dnm.genes=read.delim("de_novo_enrichment_tests_per_gene.with_sex_bias_test.txt",header=T,as.is=T)

bonferroni.cutoff=0.05/(6*804 + 2*(19685-804))
bonferroni.cutoff
#[1] 1.174095e-06

# define the genes that pass genome-wide significance
tada.significant.genes=tada.genes[tada.genes$pval<bonferroni.cutoff|tada.genes$pval.lof <bonferroni.cutoff,"gene.id"]
dnm.significant.genes.females=dnm.genes[dnm.genes$min.pval.F < bonferroni.cutoff,"gene"]
dnm.significant.genes.sexes.combined=dnm.genes[dnm.genes$min.pval.mega < bonferroni.cutoff,"gene"]
all.significant.genes=unique(c(tada.significant.genes,dnm.significant.genes.females,dnm.significant.genes.sexes.combined))
# how many pass genome-wide significance with the old method?
bonferroni.cutoff.old=0.05/(2*19685)
old.significant.genes=dnm.genes[dnm.genes$min.pval.mega < bonferroni.cutoff.old,"gene"]
length(old.significant.genes)
#[1] 19
# how many pass genome-wide significance with the new set of tests?
length(all.significant.genes)
#[1] 23

dnm.genes=dnm.genes[dnm.genes$gene %in% all.significant.genes,]
#determine ratio of observed to expected variants per gene (i.e. burden)
dnm.genes$ratio.M=dnm.genes$obs.lof.and.func.M/dnm.genes$exp.lof.and.func.M
dnm.genes$ratio.F=dnm.genes$obs.lof.and.func.F/dnm.genes$exp.lof.and.func.F
#for genes with 0 observed de novos in one sex or the other, need to add in a pseudocount since we will plot on a log scale and we want these genes to be included
dnm.genes$ratio.F2=dnm.genes$ratio.F
dnm.genes$ratio.F2[dnm.genes$ratio.F2==0]=2
dnm.genes$ratio.M2=dnm.genes$ratio.M
dnm.genes$ratio.M2[dnm.genes$ratio.M2==0]=2
#colour genes by the inheritance pattern reported in the literature
dnm.genes$col="blue"
dnm.genes$col[dnm.genes$gene %in% gene.sets$Hemi.only]="#78c679"
dnm.genes$col[dnm.genes$gene %in% gene.sets$XLD.only]="orange"

# now plot burden of de novos in females versus males, for genome-wide significant genes
pdf("Figure_2A.significant_X_genes_de_novo_enrichment_in_females_versus_males.pdf",height=6,width=9,useDingbats=F)
plot(dnm.genes$ratio.M2,dnm.genes$ratio.F2,cex=0.5,pch=19,col=dnm.genes$col,log="xy",xaxt="null",ylab="observed/expected damaging de novos in females",xlab="observed/expected damaging de novos in males",yaxt="null")
axis(1,c(3,5,10,20,50,100))
axis(2,c(3,5,10,20,50,100,200))
abline(a=0,b=1)
text(x=dnm.genes$ratio.M2,y=1.1*dnm.genes$ratio.F2,dnm.genes$gene,col=dnm.genes$col)
dev.off()


### Figure 2B: Burden of damaging (PTV + missense/inframe) DNMs for males and females in the indicated gene sets.
# start by creating a data frame with the right sets of genes
burden.combined=rbind(burden.all.dnms[burden.all.dnms$csq!="LoF.and.functional",],burden.all.cc[burden.all.cc$csq!="LoF.and.functional",colnames(burden.all.dnms)])

burden.case.control$sex="M"
burden.case.control$type="case-control"
cols.to.keep=c("class","ratio.rate.vars","burden.lt","burden.ut","excess.vars.over.N","att.frac.lt","att.frac.ut","pois.p.ut","ppv.lt","ppv.ut","sex","type")
burden.all.cc.Hemi=burden.case.control[c( "synonymous.0_hemi_gnomAD_Hemi.only" ,"functional.0_hemi_gnomAD_Hemi.only" ,"LoF.0_hemi_gnomAD_Hemi.only","LoF.and.functional.0_hemi_gnomAD_Hemi.only" ),cols.to.keep]
burden.all.cc.Hemi$csq=c("syn","functional","lof","LoF.and.functional")
burden.all.cc.XLD=burden.case.control[c( "synonymous.0_hemi_gnomAD_XLD.only" ,"functional.0_hemi_gnomAD_XLD.only" ,"LoF.0_hemi_gnomAD_XLD.only","LoF.and.functional.0_hemi_gnomAD_XLD.only"),cols.to.keep]
burden.all.cc.XLD$csq=c("syn","functional","lof","LoF.and.functional")
colnames(burden.all.cc.XLD)=c("class","ratio","burden.lt","burden.ut","excess.over.N","att.frac.lt","att.frac.ut","ppois.ut","ppv.lt","ppv.ut","sex","type","csq")
colnames(burden.all.cc.Hemi)=c("class","ratio","burden.lt","burden.ut","excess.over.N","att.frac.lt","att.frac.ut","ppois.ut","ppv.lt","ppv.ut","sex","type","csq")

cols.to.keep2=c("sex","csq","ratio","burden.lt","burden.ut","excess.over.N","att.frac.lt","att.frac.ut","ppois.ut","type")
burden.dnms$type="de novo"
burden.all.dnms.Hemi=burden.dnms[burden.dnms$geneset=="Hemi.only" & burden.dnms$csq %in% c("syn","lof","functional","LoF.and.functional") & burden.dnms$sex %in% c("M","F"),cols.to.keep2]
burden.all.dnms.XLD=burden.dnms[burden.dnms$geneset=="XLD.only" & burden.dnms$csq %in% c("syn","lof","functional","LoF.and.functional") & burden.dnms$sex %in% c("M","F"),cols.to.keep2]

burden.combined.Hemi=rbind(burden.all.dnms.Hemi[burden.all.dnms.Hemi$csq!="LoF.and.functional",],burden.all.cc.Hemi[burden.all.cc.Hemi$csq!="LoF.and.functional",colnames(burden.all.dnms.Hemi)])
burden.combined.XLD=rbind(burden.all.dnms.XLD[burden.all.dnms.XLD$csq!="LoF.and.functional",],burden.all.cc.XLD[burden.all.cc.XLD$csq!="LoF.and.functional",colnames(burden.all.dnms.XLD)])

burden.combinedB=rbind(burden.all.dnms[burden.all.dnms$csq=="LoF.and.functional",],burden.all.cc[burden.all.cc$csq=="LoF.and.functional",colnames(burden.all.dnms)])
burden.combined.HemiB=rbind(burden.all.dnms.Hemi[burden.all.dnms.Hemi$csq=="LoF.and.functional",],burden.all.cc.Hemi[burden.all.cc.Hemi$csq=="LoF.and.functional",colnames(burden.all.dnms.Hemi)])
burden.combined.XLDB=rbind(burden.all.dnms.XLD[burden.all.dnms.XLD$csq=="LoF.and.functional",],burden.all.cc.XLD[burden.all.cc.XLD$csq=="LoF.and.functional",colnames(burden.all.dnms.XLD)])

burden1=rbind(burden.combinedB[burden.combinedB$type=="de novo",])
burden1$geneset="all genes"
burden2=rbind(burden.combined.XLDB[burden.combined.XLDB$type=="de novo",])
burden2$geneset="X-linked dominant genes"
burden3=rbind(burden.combined.HemiB[burden.combined.HemiB$type=="de novo",])
burden3$geneset="X-linked recessive genes"
burden=rbind(burden1,burden2,burden3)
cols=c("#A59E98","#FFA500","#78C679")
names(cols)=c("all genes","X-linked dominant genes","X-linked recessive genes")
burden$cols=cols[burden$geneset]
burden=burden[order(burden$sex),]

#now make the plot of burden of de novos in males and females in different gene sets
pdf("Figure_2B.total_burden_in_all_versus_DDG2P_Xlinked_genes.with_error_bars.pdf",height=6,width=6,useDingbats=F)
par(mfrow=c(1,1))
x=barplot(burden$ratio,col=burden$cols,legend.text=names(cols),names.arg=c(NA,"females",NA,NA,"males",NA),ylab="observed/expected de novos",args.legend=list(x="topright",bty="n"),ylim=c(0,max(burden$ratio)+10))
abline(h=1)
for(i in 1:nrow(burden)){
    segments(x0=x[i],x1=x[i],y0=burden[i,"burden.lt"],y1=burden[i,"burden.ut"],lwd=1.5)
}
text(x=x,y=burden$ratio+2,paste("p=",formatC(burden$ppois.ut,1),"\nAttFr=",formatC(100*burden$excess.over.N,2),"%",sep=""))
dev.off()

#### Supplementary Figure 1: Histogram of the number of affected organ systems in male versus female probands in DDD.
# This figure is made in script1.compare_phenotypes_between_sexes_in_DDD.R

#### Supplementary Figure 2: Histograms indicating the distribution of bootstrapped values for metrics reported in the main text, from 10,000 iterations.

# Plot the results from boostrapping indicating the observed values and confidence intervals. Note that we are excluding some outliers here for each of visualisation
pdf("Supp_Figure_2.results_from_bootstrapping.all_metrics.pdf",height=10,width=10)
par(mfrow=c(2,2))
hist(frac.known.m,breaks=100,xlab="values from bootstraps",     main="Proportion of excess in males that is in known genes",col="black")
abline(v=c(prop.known.males,quantile(frac.known.m,c(0.025,0.975))),col=c("red","purple","purple"),lwd=2)
legend("topright",c( "observed","95% confidence limits"),lty=1,col=c("red","purple"),lwd=2)

hist(frac.known.f,breaks=100,main="Proportion of excess in females that is in known genes",col="black",xlab="values from bootstraps")
abline(v=c( prop.known.females,quantile(frac.known.f,c(0.025,0.975))),col=c("red","purple","purple"),lwd=2)
legend("topright",c( "observed","95% confidence limits"),lty=1,col=c("red","purple"),lwd=2)

hist(fraction.denovo.all[fraction.denovo.all> quantile(fraction.denovo.all,0.004) & fraction.denovo.all < quantile(fraction.denovo.all,0.995)],breaks=100,main="Proportion of excess in males that is de novo",
     xlim=quantile(fraction.denovo.all,c(0.004,0.995)),
     col="black",xlab="values from bootstraps")
abline(v=c(prop.dnm.in.males,quantile(fraction.denovo.all,c(0.025,0.975))),col=c("red","purple","purple"),lwd=2)
legend("topright",c( "observed","95% confidence limits"),lty=1,col=c("red","purple"),lwd=2)

hist(fraction.denovo.Hemi.only[fraction.denovo.Hemi.only>quantile(fraction.denovo.Hemi.only,0.001) & fraction.denovo.Hemi.only<quantile(fraction.denovo.Hemi.only,0.99)],
     breaks=100,xlim=quantile(fraction.denovo.Hemi.only,c(0.001,0.99)),main="Proportion of excess in X-linked recessive genes in males that is de novo",col="black",xlab="values from bootstraps")
abline(v=c(expectation.under.Haldane,prop.dnm.in.males.XLR,quantile(fraction.denovo.Hemi.only,c(0.025,0.975))),col=c("orange","red","purple","purple"),lwd=2)
legend("topright",c("expected under Haldane's rule", "observed","95% confidence limits"),lty=1,col=c("orange","red","purple"),lwd=2)

dev.off()

#### Supplementary Figure 3 - Positive predictive values estimated for inherited variants in trio probands with different filters.
# Let's first construct data frames with the necessary results
# Results when only variants with 0 hemizygous in gnomAD are retained:
myrows=c("inherited_only.trios_only.LoF.0_hemi_gnomAD","inherited_only.trios_only.missense.0_hemi_gnomAD",         "inherited_only.trios_only.missense.MPC.lt.1.0_hemi_gnomAD",
         "inherited_only.trios_only.missense.MPC.1.to.2.0_hemi_gnomAD",  "inherited_only.trios_only.missense.MPC.1.to.2.CADDgt25.0_hemi_gnomAD", "inherited_only.trios_only.missense.MPC.gt.2.0_hemi_gnomAD",
         "inherited_only.trios_only.missense.MPC.gt.2.CADDgt25.0_hemi_gnomAD")
compare.ppv.inherited=NULL
for(maf.cutoff in c(0.001,0.0001,0.00005)){
  if(maf.cutoff<0.001){
    inher.ppv2=burden.case.control[paste(myrows,".MAF_lt_",maf.cutoff,"_XL",sep=""),]
  } else {
    inher.ppv2=burden.case.control[paste(myrows,"_XL",sep=""),]
  }

inher.ppv2=inher.ppv2[,c("n.vars.dads","n.vars.probands","prop.vars.pathogenic","ppv.lt","ppv.ut")]
colnames(inher.ppv2)=paste(gsub("prop.vars.pathogenic","ppv",colnames(inher.ppv2)),".MAF_lt_",maf.cutoff,sep="")
  if(is.null(compare.ppv.inherited)){
    compare.ppv.inherited=inher.ppv2
  }else{
    compare.ppv.inherited=cbind(compare.ppv.inherited,inher.ppv2)
  }
}
# Results when variants are only filtered by overall MAF, ignoring the number of hemizygotes in gnomAD
myrows=c("inherited_only.trios_only.LoF",          "inherited_only.trios_only.missense",          "inherited_only.trios_only.missense.MPC.lt.1",
         "inherited_only.trios_only.missense.MPC.1.to.2",  "inherited_only.trios_only.missense.MPC.1.to.2.CADDgt25","inherited_only.trios_only.missense.MPC.gt.2",  "inherited_only.trios_only.missense.MPC.gt.2.CADDgt25")
compare.ppv.inherited2=NULL
for(maf.cutoff in c(0.001,0.0001,0.00005)){
  if(maf.cutoff<0.001){
    inher.ppv2=burden.case.control[paste(myrows,".MAF_lt_",maf.cutoff,"_XL",sep=""),]
  } else {
    inher.ppv2=burden.case.control[paste(myrows,"_XL",sep=""),]
  }
  inher.ppv2=inher.ppv2[,c("n.vars.dads","n.vars.probands","prop.vars.pathogenic","ppv.lt","ppv.ut")]
  colnames(inher.ppv2)=paste(gsub("prop.vars.pathogenic","ppv",colnames(inher.ppv2)),".MAF_lt_",maf.cutoff,sep="")
  if(is.null(compare.ppv.inherited2)){
    compare.ppv.inherited2=inher.ppv2
  }else{
    compare.ppv.inherited2=cbind(compare.ppv.inherited2,inher.ppv2)
  }
}
# set up a legend
mynames=c("PTV","missense","missense MPC<1","missense 1<MPC<2","missense 1<MPC<2, CADD>25","missense MPC>2","missense MPC>2, CADD>25")
cols=c("#98D7E5","#D37986","#FFF7EC","#FDD49E","#FDBB84","#EF6F3C","#EA3A1C")

# Plot the positive predictive value for inherited variants in males with various filters
pdf("Supp_Figure_3.ppv_for_inherited_variants.pdf",height=10,width=8,useDingbats=F)
par(mfrow=c(2,1))

x=barplot(c(compare.ppv.inherited2$ppv.MAF_lt_0.001,compare.ppv.inherited2[,"ppv.MAF_lt_1e-04"],            compare.ppv.inherited2[,"ppv.MAF_lt_5e-05"]),col=c(cols,cols,cols),
          ylab="positive predictive value",xlab="",ylim=c(-0.1,1.2),main="no filter on # hemizygotes in gnomAD")#legend=mynames,args.legend=list("top"),
axis(1,at=x[c(4,11,18),],labels=c("MAF<0.001","MAF<0.0001","MAF<0.00005"),tick=FALSE)
cis=cbind(c(compare.ppv.inherited2[,c("ppv.lt.MAF_lt_0.001")],            compare.ppv.inherited2[,c("ppv.lt.MAF_lt_1e-04")],            compare.ppv.inherited2[,c("ppv.lt.MAF_lt_5e-05")]),
          c(compare.ppv.inherited2[,c("ppv.ut.MAF_lt_0.001")],           compare.ppv.inherited2[,c("ppv.ut.MAF_lt_1e-04")],             compare.ppv.inherited2[,c("ppv.ut.MAF_lt_5e-05")]))

text(x=x,y=cis[,2]+0.05,paste(round(100*c(compare.ppv.inherited2$ppv.MAF_lt_0.001,compare.ppv.inherited2[,"ppv.MAF_lt_1e-04"],
                                          compare.ppv.inherited2[,"ppv.MAF_lt_5e-05"])),"%",sep=""))

for(i in 1:nrow(cis)){
  segments(x0=x[i],x1=x[i],y0=cis[i,1],y1=cis[i,2],lwd=1.5)
}


x=barplot(c(compare.ppv.inherited$ppv.MAF_lt_0.001,compare.ppv.inherited[,"ppv.MAF_lt_1e-04"],
            compare.ppv.inherited[,"ppv.MAF_lt_5e-05"]),col=c(cols,cols,cols),
          ylab="positive predictive value",xlab="",ylim=c(-0.1,1.2),main="0 hemizygotes in gnomAD")

cis=cbind(c(compare.ppv.inherited[,c("ppv.lt.MAF_lt_0.001")],           compare.ppv.inherited[,c("ppv.lt.MAF_lt_1e-04")],           compare.ppv.inherited[,c("ppv.lt.MAF_lt_5e-05")]),
          c(compare.ppv.inherited[,c("ppv.ut.MAF_lt_0.001")],             compare.ppv.inherited[,c("ppv.ut.MAF_lt_1e-04")],             compare.ppv.inherited[,c("ppv.ut.MAF_lt_5e-05")]))

text(x=x,y=cis[,2]+0.05,paste(round(100*c(compare.ppv.inherited$ppv.MAF_lt_0.001,compare.ppv.inherited[,"ppv.MAF_lt_1e-04"],
                                          compare.ppv.inherited[,"ppv.MAF_lt_5e-05"])),"%",sep=""))
for(i in 1:nrow(cis)){
  segments(x0=x[i],x1=x[i],y0=cis[i,1],y1=cis[i,2],lwd=1.5)
}
axis(1,at=x[c(4,11,18),],labels=c("MAF<0.001","MAF<0.0001","MAF<0.00005"),tick=FALSE)
dev.off()


######################################################
#### Supplementary Figure 5: Positive predictive value for different classes of variants in X-linked recessive versus X-linked dominant genes. 
my.geneset="Hemi.only"
inher.ppv=burden.case.control[paste(c("inherited_only.trios_only.LoF.0_hemi_gnomAD", "inherited_only.trios_only.missense.0_hemi_gnomAD",
                                "inherited_only.trios_only.missense.MPC.gt.2.CADDgt25.0_hemi_gnomAD"  ),my.geneset,sep="_"),]

dnm.ppv=burden.dnms[burden.dnms$geneset==my.geneset & burden.dnms$csq %in% c("lof","mis","functional"),
                    c("sex","excess","obs","csq","exp","ppv.lt","ppv.ut")]
dnm.ppv$ppv=as.numeric(dnm.ppv$excess)/as.numeric(dnm.ppv$obs)
dnm.ppv=dnm.ppv[dnm.ppv$sex!="combined" &dnm.ppv$csq %in% c('lof','mis'),]

my.geneset2="XLD.only"
dnm.ppv2=burden.dnms[burden.dnms$geneset==my.geneset2 & burden.dnms$csq %in% c("lof","mis","functional"),
                    c("sex","excess","obs","csq","exp","ppv.lt","ppv.ut")]

dnm.ppv2$ppv=as.numeric(dnm.ppv2$excess)/as.numeric(dnm.ppv2$obs)
dnm.ppv2=dnm.ppv2[dnm.ppv2$sex!="combined" &dnm.ppv2$csq %in% c('lof','mis'),]

dnm.ppv=dnm.ppv[4:1,]
dnm.ppv2=dnm.ppv2[4:1,]

mynames=c("de novo PTV in females","de novo missense in females","de novo PTV in males","de novo missense in males","inherited PTV","inherited missense","inherited missense MPC>2, CADD>25")
cols=c("#7F99CC","#EF5273","#98D7E5","#D37986","#98D7E5","#D37986","#EAB9C3")

pdf("Supp_Figure_5.ppv.split_by_XLR_versus_XLD.pdf",height=8,width=8,useDingbats=F)
x=barplot(c(dnm.ppv$ppv,inher.ppv$prop.vars.pathogenic,dnm.ppv2$ppv),col=c(cols,cols),ylab="positive predictive value",xlab="",ylim=c(-0.2,1.6))
legend("top",mynames,fill=cols)
#legend("top",mynames,angle=c(rep(45,4),NA,NA,NA),density=c(rep(15,4),NA,NA,NA),col=c(rep("black",4),NA,NA,NA))
legend("top",mynames,angle=45,density=15,fill=c(rep("black",4),cols[5:7]))

cis=rbind(dnm.ppv[,c("ppv.lt","ppv.ut")],inher.ppv[,c("ppv.lt","ppv.ut")],dnm.ppv2[,c("ppv.lt","ppv.ut")])
axis(1,x[c(4,10),1],c("X-linked recessive only","X-linked dominant only"))

text(x=x,y=cis[,2]+0.05,paste(round(100*c(dnm.ppv$ppv,inher.ppv$prop.vars.pathogenic,dnm.ppv2$ppv)),"%",sep=""))
barplot(c(dnm.ppv$ppv,0,0,0,dnm.ppv2$ppv),add=T,angle=45,density=15,col="black")

for(i in 1:nrow(cis)){
  segments(x0=x[i],x1=x[i],y0=cis[i,1],y1=cis[i,2],lwd=1.5)
}
dev.off()

#### Supplementary Table 1: Results from logistic regression comparing the prevalence of different phenotypic features amongst males versus females in DDD, controlling for age at assessment.
#This is the file Supp_Table_1.results_from_logistic_regression_comparing_presence_absence_of_phenotypes_between_sexes_controlling_for_age_at_assessment.txt output from script1.compare_phenotypes_between_sexes_in_DDD.R


#### Supplementary Table 2: Results from comparing growth metrics and developmental milestones between males and females in DDD using linear regression.
#This is the file Supp_Table_2.results_of_comparing_milestones_and_growth_metrics_between_sexes.milestones_on_for_those_achieved_already.controlling_for_age_at_ assessment.txt output from script1.compare_phenotypes_between_sexes_in_DDD.R


#### Supplementary Table 3: De novo mutations that passed our filtering and were used in the burden analysis and gene-based tests.
#This is from de_novos_in_both_sexes.post_filtering.0_hemis_gnomAD.MAX_AF_0.001.RData which is the input to script4.prepare_data_for_TADA.R


#### Supplementary Table 4: Variant counts and results from gene-based tests using TADA to analyse de novo and inherited variants in males, or de novo enrichment tests in males and females
# This is made by merging the tables output_from_TADA.frac_risk_genes_0.05.txt output from script5.run_TADA_to_do_per_gene_tests_on_males.R and de_novo_enrichment_tests_per_gene.with_sex_bias_test.txt output from script3.burden_analysis_and_per_gene_tests_on_de_novos.R


tada.output=read.delim("output_from_TADA.frac_risk_genes_0.05.txt",header=T,as.is=T)

tada.output=tada.output[,c("gene.id","case.cls1","case.cls2","ctrl.cls1","ctrl.cls2","BF.cls1","BF.cls2","BF.total","qval.lof","qval","pval.lof","pval")]
colnames(tada.output)=c("gene","n.lof.males.excl.dnms","n.func.males.excl.dnms","n.lof.fathers","n.func.fathers","BF.lof","BF.func","BF.combined","qval.lof","qval.combined","pval.lof","pval.combined")

colnames(tada.output)[6:ncol(tada.output)]=paste("TADA.",colnames(tada.output)[6:ncol(tada.output)],sep="")
                       
dnm.tests=read.delim("de_novo_enrichment_tests_per_gene.with_sex_bias_test.txt",header=T,as.is=T)
dnm.tests=dnm.tests[, c("gene","obs.lof.F","obs.functional.F","obs.lof.M","obs.functional.M","ppois.lof.mega","ppois.lof.and.func.mega" , "ppois.lof.F", "ppois.lof.and.func.F", "ppois.lof.M","ppois.lof.and.func.M"  )]

colnames(dnm.tests)=gsub("obs","n.dnms",colnames(dnm.tests))
colnames(dnm.tests)=gsub("ppois","dnm.pval",colnames(dnm.tests))

gene.results=merge(tada.output,dnm.tests,by.x="gene",by.y="gene")
gene.results=gene.results[,c(1,13:16,2:12,17:ncol(gene.results))]

colnames(gene.results)=gsub("lof","PTV",colnames(gene.results))
colnames(gene.results)=gsub("functional","missenseinframe",colnames(gene.results))
colnames(gene.results)=gsub("func","missenseinframe",colnames(gene.results))

write.table(gene.results,"Supplementary_Table_4.tests_per_gene.txt",quote=F,sep="\t",row.names=F)

#### Supplementary Table 5: Results of tests for a depletion of DNMs in males compared to females, likely due to male lethality.
# This is in the file de_novo_enrichment_tests_per_gene.with_sex_bias_test.txt output from script3.burden_analysis_and_per_gene_tests_on_de_novos.R 


#### Supplementary Table 6: Results from logistic regressions testing the difference in polygenic scores for the indicated traits between 216 males who were suspected by clinicians to have X-linked inheritance versus 3439 who were not.
#This is the file Supp_Table_6.PRS_comparisons.txt from script7.compare_PRS_between_patient_subsets_for_chrX_paper.R


#### Supplementary Table 7 : Parameters used for TADA.
# This is made using the file Supp_Table_7.parameters_run_with_TADA.txt from script5.run_TADA_to_do_per_gene_tests_on_males.R


