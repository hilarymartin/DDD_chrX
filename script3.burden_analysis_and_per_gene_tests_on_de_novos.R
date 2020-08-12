setwd("/lustre/scratch115/projects/ddd/users/hcm/DDD/chrX_analyses/code_and_data_to_release/")
#### Note that this script is quite slow and will take several hours because of the boostrapping part. If you don't want to run the bootstrapping part to get confidence intervals, just change the if(TRUE){ line below to if(FALSE){

n.iterations=100 # number of iterations for bootstrapping
library(plyr)

### Load lists of genes - to use in burden analysis
# DDG2P/OMIM consensus classifications --> gene.sets
load("gene_sets_to_use_in_burden_analysis.DDG2P_14_1_2020.consensus_with_OMIM.RData")
# list of possible DDG2P X-linked genes  - just to use in annotating output --> gene.sets2
load("gene_sets_to_use_in_annotation.DDG2P_possible_14_1_2020.RData")

### Load lists of deidentified IDs of males and females to use in bootstrapping
load("deidentified_IDs_of_male_and_female_trio_probands_to_use_in_bootstrapping_for_de_novo_analysis.RData")

### Load observed de novos
load("de_novos_in_both_sexes.post_filtering.0_hemis_gnomAD.MAX_AF_0.001.RData")

### Load mutation rates for chrX genes, adjusted for coverage using the coverage data from the 250 males and 250 females in DDD
load("mutation_rates_per_gene_using_longest_Ensembl_transcript.chrX_nonPAR_genes.adjusted_for_sex_specific_coverage_in_DDD.RData")
males.x.mut.rates = adjusted.rates.per.longest.gene[["M"]]
females.x.mut.rates = adjusted.rates.per.longest.gene[["F"]]

### Define sample sizes - number of trios with male and female probands
male.n=5138
female.n=3908

#### Calculate the expected number of mutations in males and females accounting for ploidy
x.female.transmissions = male.n + female.n
x.male.transmissions = female.n
alpha = 3.4
male.factor = 2 / (1 + (1 / alpha))
female.factor = 2 / (1 + alpha)
x.male.transmissions = (x.male.transmissions * male.factor)
x.female.transmissions = (x.female.transmissions * female.factor)

x.female.transmissions.to.males= male.n* female.factor
x.male.transmissions.to.males=0

x.female.transmissions.to.females= female.n* female.factor
x.male.transmissions.to.females=female.n* male.factor

# Note that we have excluded PAR genes already from chrX so will now treat all chrX genes as haploid in males
adj.expected.denovos.in.males=(x.female.transmissions.to.males) *males.x.mut.rates[,grep("p_",colnames(males.x.mut.rates))]
rownames(adj.expected.denovos.in.males) = c(males.x.mut.rates[,"gene"])
adj.expected.denovos.in.males=as.data.frame(adj.expected.denovos.in.males)
adj.expected.denovos.in.males$cds.length=c(males.x.mut.rates[,"bp"])

adj.expected.denovos.in.females=(x.female.transmissions.to.females+ x.male.transmissions.to.females) *females.x.mut.rates[,grep("p_",colnames(females.x.mut.rates))]
adj.expected.denovos.in.females=as.data.frame(adj.expected.denovos.in.females)
adj.expected.denovos.in.females$cds.length=c(females.x.mut.rates[,"bp"])
rownames(adj.expected.denovos.in.females) = c(females.x.mut.rates[,"gene"])

### Add expected  mutation rates for frameshifts and inframe indels
# The calculation is based on the probability of a nonsense variant in the gene:
#      P_fs = P_non*scaling_factor
#      P_fs.per.gene = (P_fs*bp) / BPtot
#      P_inframe.indel.per.gene = P_fs.per.gene/9
# Where P_fs is the overall probability of a frameshift mutation, P_non is the overall probability of a nonsense mutation, P_fs.per.gene is the probability of a frameshift mutation in a given gene,
# P_inframe.indel.per.gene is the probability of an inframe indel per gene, bp is the number of coding base pairs in a gene, and BPtot is the total number of coding base pairs.
# We will calculate the scaling factor based on the frameshift:nonsense ratio in the DDD 4K dataset:
nonsense.n = 411
frameshift.n = 610
scaling.factor = frameshift.n/nonsense.n
total.framshift.rate.males=sum(adj.expected.denovos.in.males$p_non) * scaling.factor
total.framshift.rate.females=sum(adj.expected.denovos.in.females$p_non) * scaling.factor
total.cds.length= sum(males.x.mut.rates$bp)

adj.expected.denovos.in.males$p_frameshift=total.framshift.rate.males * adj.expected.denovos.in.males$cds.length/total.cds.length
adj.expected.denovos.in.males$p_missense.indel=adj.expected.denovos.in.males$p_frameshift/9
adj.expected.denovos.in.females$p_frameshift=total.framshift.rate.females * adj.expected.denovos.in.females$cds.length/total.cds.length
adj.expected.denovos.in.females$p_missense.indel=adj.expected.denovos.in.females$p_frameshift/9

# Sum together rates for consensus splice site, frameshift and nonsense mutations to get the overall LoF rate
# Sum together rate for missense and inframe indels to get the overall rate for 'functional' mutations
adj.expected.denovos.in.males$lof=adj.expected.denovos.in.males$p_css + adj.expected.denovos.in.males$p_frameshift + adj.expected.denovos.in.males$p_non
adj.expected.denovos.in.males$functional=adj.expected.denovos.in.males$p_mis + adj.expected.denovos.in.males$p_missense.indel
adj.expected.denovos.in.females$lof=adj.expected.denovos.in.females$p_css + adj.expected.denovos.in.females$p_frameshift + adj.expected.denovos.in.females$p_non
adj.expected.denovos.in.females$functional=adj.expected.denovos.in.females$p_mis + adj.expected.denovos.in.females$p_missense.indel

##### Do the burden analysis for all chrX genes and for different sets of genes
total.gene.sets=NULL
for(i in 1:length(gene.sets)){
  cat(names(gene.sets)[i],"\n")
  #define the gene set to be used here
  genes.to.keep=gene.sets[[i]]
  # calculate the expected and observed number of de novos in this gene set by summing up rates across genes and counting the observed mutations in different classes
  total.these.genes=cbind(c(colSums(adj.expected.denovos.in.males[rownames(adj.expected.denovos.in.males) %in% genes.to.keep,c("p_syn","p_mis","p_non","p_css","p_frameshift","p_missense.indel","lof","functional")]),
    colSums(adj.expected.denovos.in.females[rownames(adj.expected.denovos.in.females) %in% genes.to.keep,c("p_syn","p_mis","p_non","p_css","p_frameshift","p_missense.indel","lof","functional")])),
    c(colSums(d.sep2[d.sep2$sex=="M" & d.sep2$symbol %in% genes.to.keep,c("is.syn","is.mis","is.nonsense","is.css","is.frameshift","is.inframe","is.lof","is.functional")]),
      colSums(d.sep2[d.sep2$sex=="F" & d.sep2$symbol %in% genes.to.keep,c("is.syn","is.mis","is.nonsense","is.css","is.frameshift","is.inframe","is.lof","is.functional")])))
    colnames(total.these.genes)=c("exp","obs")
  rownames(total.these.genes)=NULL
  total.these.genes=as.data.frame(total.these.genes)
  total.these.genes$sex=c(rep("M",8),rep("F",8))
  total.these.genes$csq=c(gsub("is.","",c("is.syn","is.mis","is.nonsense","is.css","is.frameshift","is.inframe","is.lof","is.functional")),gsub("is.","",
    c("is.syn","is.mis","is.nonsense","is.css","is.frameshift","is.inframe","is.lof","is.functional")))

#combine LoF and missense
  for(sex in c("M","F")){
    tmp=as.data.frame(total.these.genes[total.these.genes$sex %in% sex & total.these.genes$csq=="lof",1:2]+total.these.genes[total.these.genes$sex %in% sex & total.these.genes$csq=="functional",1:2])
    colnames(tmp)=colnames(total.these.genes)[1:2]
    tmp$sex=sex
    tmp$csq="LoF.and.functional"
    total.these.genes=rbind(total.these.genes,tmp)

    tmp=as.data.frame(total.these.genes[total.these.genes$sex %in% sex & total.these.genes$csq=="nonsense",1:2]+total.these.genes[total.these.genes$sex %in% sex & total.these.genes$csq=="mis",1:2]+
      total.these.genes[total.these.genes$sex %in% sex & total.these.genes$csq=="css",1:2])
    colnames(tmp)=colnames(total.these.genes)[1:2]
    tmp$sex=sex
    tmp$csq="LoF.and.missense.SNV"
    total.these.genes=rbind(total.these.genes,tmp)
  }
  #combine male and female
  for(csq in c("mis","lof","functional","LoF.and.functional","LoF.and.missense.SNV")){
    tmp=as.data.frame(t(colSums(total.these.genes[total.these.genes$csq %in% csq,1:2])))
    colnames(tmp)=colnames(total.these.genes)[1:2]
    tmp$sex="combined"
    tmp$csq=csq
    total.these.genes=rbind(total.these.genes,tmp)
  }
  #calculate the observed:expected ratio for de novos
  total.these.genes$ratio=total.these.genes$obs/total.these.genes$exp
  # calculate the excess = observed-expected
  total.these.genes$excess=total.these.genes$obs-total.these.genes$exp
  # calculate a Poisson p-value for enrichment
  total.these.genes$ppois.ut=ppois(total.these.genes$obs-1,total.these.genes$exp,lower.tail=F)
  # define the sample sizes                                        
  total.these.genes$N=male.n
  total.these.genes$N[total.these.genes$sex=="F"] = female.n
  total.these.genes$N[total.these.genes$sex=="combined"] = female.n + male.n
  # calculate attributable fraction = excess/N
  total.these.genes$excess.over.N = total.these.genes$excess/total.these.genes$N
  # label these rows with the gene set being used
  total.these.genes$geneset=names(gene.sets)[i]
  # save in a table with results from all different gene sets
  total.gene.sets=rbind(total.gene.sets,total.these.genes)
}

#now calculate confidence intervals for the attributable fraction ([observed-expected]/N), burden (observed/expected) and positive predictive value ([observed-expected]/observed)
total.gene.sets$att.frac.lt=NA
total.gene.sets$att.frac.ut=NA
total.gene.sets$burden.lt=NA
total.gene.sets$burden.ut=NA
total.gene.sets$ppv.lt=NA
total.gene.sets$ppv.ut=NA

for(i in 1:nrow(total.gene.sets)){
  #att.fraction = O/n-E/n
  exp.rate=total.gene.sets[i,"exp"]/total.gene.sets[i,"N"]
  ci.obs.rate=poisson.test(total.gene.sets[i,"obs"],T=total.gene.sets[i,"N"],r=exp.rate)$conf.int
  total.gene.sets$att.frac.lt[i] = ci.obs.rate[1]-exp.rate
  total.gene.sets$att.frac.ut[i] = ci.obs.rate[2]-exp.rate
  #burden=O/E
  total.gene.sets$burden.lt[i] = ci.obs.rate[1]/exp.rate
  total.gene.sets$burden.ut[i] = ci.obs.rate[2]/exp.rate
  #ppv=(O-E)/O = 1 - E/O = 1-(E/n)/(O/n)
  total.gene.sets$ppv.lt[i] = 1-exp.rate/ci.obs.rate[1]
  total.gene.sets$ppv.ut[i] = 1-exp.rate/ci.obs.rate[2]
  }

write.table(total.gene.sets,"burden_analysis_of_de_novos.different_chrX_gene_sets.txt",quote=F,sep="\t",row.names=F)


#######################################################
#### Bootstrap probands to get confidence intervals on key parameters
# This bit is slow - if you don't want to run it, just change the next line to: if(FALSE){
if(TRUE){
bootstrapped.trios=list()
bootstrapped.trios[["M"]]=NULL
bootstrapped.trios[["F"]]=NULL
bootstraps.by.sex=list()
all.obs=split(d.sep2,d.sep2$randomised.ID)
#### iterate over females and males
for(sex in c("M","F")){
  bootstrapped.totals=NULL
  # bootstrap sets of probands n.iterations times
  for( b in 1:n.iterations){
      if(b %/% 50 - b/50 ==0){
        cat(sex,"\t",b,"\n")
      }
     # save bootstrapped indivs for either all genes or known X-linked genes
      to.choose=probands.to.choose[[sex]]
      probs.to.use=to.choose[sample(length(to.choose),length(to.choose),replace=T)]
      bootstrapped.trios[[sex]]=cbind(bootstrapped.trios[[sex]],probs.to.use)

    #pull out the rows corresponding to these
      obs=d.sep2[unlist(sapply(sort(probs.to.use),function(x){which(d.sep2$randomised.ID %in% x)})),]
    
      total.combined=NULL
      for(i in c(1:2,5)){
        genes.to.keep=gene.sets[[i]]
    # calculate overall burden of LoF + functional, then calculate all ratio
        if(sex=="M"){
          rates=adj.expected.denovos.in.males
        } else {
          rates=adj.expected.denovos.in.females
        }
        total.these.genes=cbind(c(colSums(rates[rownames(rates) %in% genes.to.keep,c("p_syn","lof","functional","p_non", "p_mis","p_css"  )]),
          c(colSums(obs[obs$symbol %in% genes.to.keep,c("is.syn","is.lof","is.functional","is.nonsense","is.mis","is.css")]))))
  #combine LoF+functional(missense/inframe) variants, or just LoF+missense SNVs   
        total.these.genes=rbind(total.these.genes,sum(total.these.genes[c("is.lof","is.functional"),]))
        total.these.genes=rbind(total.these.genes,sum(total.these.genes[c("is.nonsense","is.css","is.mis"),]))
        total.these.genes=rbind(total.these.genes,sum(total.these.genes[c("lof","functional"),]))
        total.these.genes=rbind(total.these.genes,sum(total.these.genes[c("p_non","p_css","p_mis"),]))
        
        rownames(total.these.genes)[(nrow(total.these.genes)-3):nrow(total.these.genes)]=c("LoF.and.functional","LoF.and.missense.SNV","p_LoF.and.functional","p_LoF.and.missense.SNV")
        rownames(total.these.genes)=paste(rownames(total.these.genes),names(gene.sets)[[i]],sep=".")
        total.combined=rbind(total.combined,total.these.genes)
      }
      bootstrapped.totals=cbind(bootstrapped.totals,total.combined)
    }
  bootstrapped.totals=as.data.frame(t(bootstrapped.totals))
  #calculate excess for different gene sets, for either all LoF+functional(missense/inframe) variants, or just LoF+missense SNVs
  bootstrapped.totals$excess.LoF.and.functional.Hemi.only = bootstrapped.totals$LoF.and.functional.Hemi.only - bootstrapped.totals$p_LoF.and.functional.Hemi.only
  bootstrapped.totals$excess.LoF.and.functional.XL = bootstrapped.totals$LoF.and.functional.XL - bootstrapped.totals$p_LoF.and.functional.XL
  bootstrapped.totals$excess.LoF.and.functional.all = bootstrapped.totals$LoF.and.functional.all - bootstrapped.totals$p_LoF.and.functional.all
  bootstrapped.totals$excess.LoF.and.missense.SNV.Hemi.only = bootstrapped.totals$LoF.and.missense.SNV.Hemi.only - bootstrapped.totals$p_LoF.and.missense.SNV.Hemi.only
  bootstrapped.totals$excess.LoF.and.missense.SNV.XL = bootstrapped.totals$LoF.and.missense.SNV.XL - bootstrapped.totals$p_LoF.and.missense.SNV.XL
  bootstrapped.totals$excess.LoF.and.missense.SNV.all = bootstrapped.totals$LoF.and.missense.SNV.all - bootstrapped.totals$p_LoF.and.missense.SNV.all
  #calculate the fraction of the excess (observed-expected) that is in known X-linked DD-associated genes
  bootstrapped.totals$fraction.burden.XL=bootstrapped.totals$excess.LoF.and.functional.XL/bootstrapped.totals$excess.LoF.and.functional.all
  bootstraps.by.sex[[sex]]=bootstrapped.totals
}

#combine male and female
bootstraps.by.sex[["combined"]]=bootstraps.by.sex[["M"]]+bootstraps.by.sex[["F"]]
bootstraps.by.sex[["combined"]]$fraction.burden.XL=bootstraps.by.sex[["combined"]]$excess.LoF.and.functional.XL/bootstraps.by.sex[["combined"]]$excess.LoF.and.functional.all

#calculate confidence interval on the fraction of the burden that is in known X-linked DD-associated genes, using the bootstrapped results   
conf.ints.fraction.DNM.burden.XL=NULL
for(sex in names(bootstraps.by.sex)){
  conf.ints.fraction.DNM.burden.XL=rbind(conf.ints.fraction.DNM.burden.XL, quantile(bootstraps.by.sex[[sex]][,"fraction.burden.XL"],c(0.025,0.05,0.95,0.975)))
}
rownames(conf.ints.fraction.DNM.burden.XL)=names(bootstraps.by.sex)
write.table(conf.ints.fraction.DNM.burden.XL,"bootstrapped_confidence_intervals_on_fraction_DNMs_in_known_genes.txt",quote=F,sep="\t",row.names=F)

#save bootstrapped trios
save(bootstrapped.trios,file="trios_selected_for_bootstrapping.RData")

#save bootstrap results --> these will be used in another script
save(bootstraps.by.sex,file="results_of_bootstrapping_for_DNM_burden.RData")
}
#######################################################
#### Do gene-based tests on de novos, for males alone, females alone and then both sexes combined
for(sex in c("M","F","MF")){
  # pull out observed de novos
  if(sex %in% c("M","F")){
    obs=d.sep2[d.sep2$sex==sex,]
  }else {
    obs=d.sep2
  }
  #pull out expected de novos
  if(sex=="F"){
    exp.data=adj.expected.denovos.in.females
  }else if(sex=="M") {
    exp.data=adj.expected.denovos.in.males
  } else {
  #combine expectation for males and females
    exp.data=adj.expected.denovos.in.females[,c("lof","functional")] + adj.expected.denovos.in.males[rownames(adj.expected.denovos.in.females),c("lof","functional")]
  }
  #for LoF and functional (missense/inframe) DNMs separately
  for(csq in c("lof","functional")){
    #count the de novos in this class
    gene.count=as.data.frame.table(table(obs[obs[,paste('is.',csq,sep="")],"symbol"])) 
    gene.count=gene.count[gene.count[,1]!=".",]
    ### merge observed with expected
    exp.data$newobs=0
    exp.data[as.character(gene.count[,1]),"newobs"] = gene.count[,2]
    colnames(exp.data)[ncol(exp.data)] = paste("obs.",csq,sep="")
  }
  # calculate total observed and expected for LoF+missense/inframe DNMs combined                                        
  exp.data$obs.lof.and.func = exp.data$obs.lof + exp.data$obs.functional
  exp.data$exp.lof.and.func= exp.data$lof + exp.data$functional

  # conduct Poisson test of enriched for LoFs
  exp.data$ppois.lof=ppois(exp.data$obs.lof - 1,exp.data$lof,lower.tail=F)
  #  LoF+missense/inframe
  exp.data$ppois.lof.and.func=ppois(exp.data$obs.lof.and.func - 1,exp.data$exp.lof.and.func,lower.tail=F)
  # select lowest p-value out of these two tests
  exp.data$min.pval=apply(exp.data[,c("ppois.lof","ppois.lof.and.func")],1,min,na.rm=T)
  # define a new data frame with the observed and expected mutations per gene
  if(sex=="F"){
    adj.expected.denovos.in.females=exp.data
  }else if(sex=="M") {
    adj.expected.denovos.in.males=exp.data
  }else{
    adj.expected.denovos.in.megaanalysis=exp.data
  }
}
#add column names and gene symbols to the data frames with the observed and expected mutations per gene  
colnames(adj.expected.denovos.in.males)=paste(colnames(adj.expected.denovos.in.males),"M",sep=".")
colnames(adj.expected.denovos.in.females)=paste(colnames(adj.expected.denovos.in.females),"F",sep=".")
colnames(adj.expected.denovos.in.megaanalysis)=paste(colnames(adj.expected.denovos.in.megaanalysis),"mega",sep=".")
adj.expected.denovos.in.males$gene=rownames(adj.expected.denovos.in.males)
adj.expected.denovos.in.females$gene=rownames(adj.expected.denovos.in.females)
adj.expected.denovos.in.megaanalysis$gene=rownames(adj.expected.denovos.in.megaanalysis)

#combine results for males, females and the two sexes together
combined=merge(adj.expected.denovos.in.males,adj.expected.denovos.in.females,by.x="gene",by.y='gene',all=T)
combined=merge(combined,adj.expected.denovos.in.megaanalysis,by.x="gene",by.y='gene',all=T)
combined$min.pval.F[!is.finite(combined$min.pval.F)]= NA
combined$min.pval.M[!is.finite(combined$min.pval.M)]= NA
combined$min.pval.mega[!is.finite(combined$min.pval.mega)]= NA

# Annotate table of p-values with gene classifications according to the consensus of DDG2P and OMIM
combined$ddg2p.prob.conf=NA
combined$ddg2p.prob.conf[combined$gene %in% gene.sets[["Hemi.only"]]] = "hemizygous only"
combined$ddg2p.prob.conf[combined$gene %in% gene.sets[["XLD.only"]]] = "XLD only"
combined$ddg2p.prob.conf[combined$gene %in% gene.sets[["Hemi.and.XLD"]]] = "hemizygous and XLD"
combined$ddg2p.prob.conf[combined$gene %in% gene.sets[["XLOD"]]] = "XL over-dominant"

# Annotate table of p-values with gene classifications with DDG2P "possible enes
combined$ddg2p.possible=NA
combined$ddg2p.possible[combined$gene %in% gene.sets2[["Hemi.only"]]] = "hemizygous only"
combined$ddg2p.possible[combined$gene %in% gene.sets2[["XLD.only"]]] = "XLD only"
combined$ddg2p.possible[combined$gene %in% gene.sets2[["Hemi.and.XLD"]]] = "hemizygous and XLD"
combined$ddg2p.possible[combined$gene %in% gene.sets2[["XLOD"]]] = "XL over-dominant"

### Test if there's a depletion of the fraction of DNMs that are in males compared to what we would expect given the mutation rates, to try to identify X-linked semi-dominant genes
# First we need to calculate the fraction of expected de novos that we expect in males per gene
# Note that, because the rates of expected mutations are coverage-corrected, this differs slightly between genes due to differing coverage levels
exp.fraction.male.per.gene=(adj.expected.denovos.in.males$lof.M+adj.expected.denovos.in.males$functional.M)/(adj.expected.denovos.in.males$lof.M+adj.expected.denovos.in.males$functional.M
                                                                                                            + adj.expected.denovos.in.females$lof.F+adj.expected.denovos.in.females$functional.F)
names(exp.fraction.male.per.gene)=rownames(adj.expected.denovos.in.males)

combined$p.depletion.M.DNMs = NA
combined$exp.fraction.female.DNMs = NA
combined$obs.fraction.female.DNMs = NA
for(i in 1:nrow(combined)){
  if(!is.na(combined[i,"obs.lof.and.func.F"] ) & !is.na(combined[i,"obs.lof.and.func.M"])){
    #if there is at least 1 LoF or missense/inframe DNM in either males or females
    if((combined[i,"obs.lof.and.func.F"] + combined[i,"obs.lof.and.func.M"])>0){
      # compare the fraction of observed DNMs that are in males to the expected fraction using a lower-tailed binomial test
      combined$p.depletion.M.DNMs[i] = prop.test(combined[i,"obs.lof.and.func.M"],combined[i,"obs.lof.and.func.F"] + combined[i,"obs.lof.and.func.M"],exp.fraction.male.per.gene[combined$gene[i]],alternative="less")$p.value
      # save the observed and expected fraction of DNMs that are in females
      combined$exp.fraction.female.DNMs[i]=1-exp.fraction.male.per.gene[combined$gene[i]]
      combined$obs.fraction.female.DNMs[i]=combined[i,"obs.lof.and.func.F"]/(combined[i,"obs.lof.and.func.F"] + combined[i,"obs.lof.and.func.M"])
    }
  }
}


write.table(combined,"de_novo_enrichment_tests_per_gene.with_sex_bias_test.txt",quote=F,sep="\t",row.names=F )


