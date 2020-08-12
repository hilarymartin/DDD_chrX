setwd("/lustre/scratch115/projects/ddd/users/hcm/DDD/chrX_analyses/code_and_data_to_release/")
### Note: This script will probably take >48 hours to run because it involves bootstrapping with a 'for' loop. If you want to speed it up by skipping the bootstrapping, change the if(TRUE){ line below to if(FALSE){

n.iterations=100 #Number of iterations for boostrapping. 

library(data.table)
library(ratesci)

### load gene sets
load("gene_sets_to_use_in_burden_analysis.DDG2P_14_1_2020.consensus_with_OMIM.RData")

load("annotated_variants_for_case_control_analysis.RData")   
d.indicating.denovos=d

### define sample sizes
n.dads.all=8551
n.probands.all=7136 
n.trio.probands.all=5138
n.probands.candidate.xlinked=273
n.trio.probands.candidate.xlinked=180

counts=NULL
classes=c() #we will loop over variant classes
sets=c() #we will loop over gene sets 
my.n.probands=c() #we will loop over different sets of probands

#define consequence classes with different filters     
my.classes=c("synonymous","functional","functional.CADDgt25","missense","missense.CADDgt25","LoF","all.LoF.SNV","missense.MPC.lt.1","missense.MPC.gt.1","missense.MPC.1.to.2","missense.MPC.gt.2",
  "missense.MPC.gt.1.CADDgt25","missense.MPC.1.to.2.CADDgt25","missense.MPC.gt.2.CADDgt25","LoF.and.functional","LoF.and.missense.SNV")
#iterate over MAF cutoffs
for(maf.cutoff in c(0.001,0.0001,0.00005)){
  #either retain all variants in the case/control analysis that passed filtering (remove.denovos=FALSE) or remove variants that passed the de novo filtering, to focus on the inherited variants (remove.denovos=TRUE)
  for(remove.denovos in c(FALSE,TRUE)){
    #we will do the analysis either on all male probands (candidate.xlinked=FALSE) or just on the male probands who were suspected by clinicians to have an X-linked inheritance pattern (candidate.xlinked=TRUE)
    candidate.xlinked.options=c(TRUE,FALSE)
    if(maf.cutoff<0.001){
      candidate.xlinked.options=c(FALSE) # we only want to do the burden analysis on the male probands who were suspected to have X-linked inheritance using the MAF<0.001 filter
    }
    for(candidate.xlinked in candidate.xlinked.options){
     #we will do the analysis either on all male probands (trios.only=FALSE) or only on the male probands who are part of trios (trios.only=TRUE) - we want the latter for analyses involving de novo calls
      trios.only.options=c(FALSE,TRUE)
      if(candidate.xlinked){
        trios.only.options=c(FALSE) # we don't care about the trios-only analysis for the male probands who were suspected to have X-linked inheritance 
      }
      if(maf.cutoff<0.001){
        trios.only.options=c(TRUE) # for the analysis with a lower MAF cutoff, we'll focus on trio probands only, since the point is to focus just on inherited variants here
      }
      for(trios.only in trios.only.options){
        #we want to try the burden analysis either restricting to variants that had 0 hemizygotes in gnomAD (remove.gnomad.hemi=TRUE) or keeping them in (remove.gnomad.hemi=FALSE)
            for(remove.gnomad.hemi in c(FALSE,TRUE)){
              #we will iterate over consequence classes with different filters
              for(myclass in my.classes){
                #we will iterate over gene sets
                for(myset in names(gene.sets)){
#                  cat("maf cutoff: ",maf.cutoff,"; removing de novos=", remove.denovos,"; considering suspected X-linked probands=",candidate.xlinked, "; considering trios only=",trios.only,
#                      "; removing variants with hemizygotes in gnomad=",remove.gnomad.hemi,"; considering ",myclass," in geneset ",myset,"\n")
                  #create a label for the consequence class being tested here
                  if(remove.gnomad.hemi){
                    tmp.class=paste(myclass,"0_hemi_gnomAD",sep=".")
                  }else {
                    tmp.class=myclass
                  }
                  if(maf.cutoff<0.001){
                    tmp.class=paste(tmp.class,".MAF_lt_",maf.cutoff,sep="")
                  }
                  #store the gene set being used
                  sets=c(sets,myset)
                  #define the sample size for probands
                  n.probands=n.probands.all
                  n.trio.probands=n.trio.probands.all
                  if(candidate.xlinked){
                    #if we are focusing on probands with suspected X-linked inheritance, define the sample sizes and add a flag to the label
                    n.probands=n.probands.candidate.xlinked
                    n.trio.probands=n.trio.probands.candidate.xlinked
                    tmp.class=paste("probands_suspected_Xlinked.",tmp.class,sep="")
                  }
                  #store the number of probands being used in the analysis
                  if(trios.only){
                      tmp.class=paste("trios_only.",tmp.class,sep="")
                      my.n.probands=c(my.n.probands,n.trio.probands)
                  }else {
                    my.n.probands=c(my.n.probands,n.probands)
                  }
                  #pull out the set of variants to count for this analysis - are we removing de novos to focus on inherited variants, or keeping all variants?
                  if(remove.denovos){
                    #### d.indicating.denovos$de.novo.passed ==> indicates whether the variant passed our stringent de novo filtering
                    #### d.indicating.denovos$absent.from.mother ==> indicates whether the ALT allele was totally absent from the mother; is NA if mother was not sequenced
                    d=d.indicating.denovos[(!d.indicating.denovos$de.novo.passed & !d.indicating.denovos$absent.from.mother) | is.na(d.indicating.denovos$absent.from.mother),]
                    #add a flag to the label for this analysis
                    tmp.class=paste("inherited_only.",tmp.class,sep="")
                  }else{
                    d=d.indicating.denovos
                  }
                  #store the consequence class being used + other labels
                  classes=c(classes,tmp.class)

                  #subset variants from data frame d                
                  if(myclass %in% c("synonymous","functional","LoF")){
                    e=d[d$class %in% myclass & d$gene %in% gene.sets[[myset]],]
                  }
                  if(myclass %in% c("LoF.and.functional")){
                    e=d[d$class %in% c("LoF","functional") & d$gene %in% gene.sets[[myset]],]
                  }
                  if(myclass %in% "functional.CADDgt25"){
                    e=d[d$class %in% "functional" & !is.na(d$CADD) & d$CADD>25 & d$gene %in% gene.sets[[myset]],]
                  }
                  if(myclass %in% "missense"){
                    e=d[1:nrow(d) %in% grep("missense",d$consequence) & d$gene %in% gene.sets[[myset]],]
                  }
                  if(myclass %in% "missense.CADDgt25"){
                    e=d[1:nrow(d) %in% grep("missense",d$consequence) & !is.na(d$CADD) & d$CADD>25 & d$gene %in% gene.sets[[myset]],]
                  }
                  if(myclass %in% "missense.MPC.lt.1"){
                    e=d[1:nrow(d) %in% grep("missense",d$consequence) & !is.na(d$MPC) & d$MPC<1 & d$gene %in% gene.sets[[myset]],]
                  }
                  if(myclass %in% "missense.MPC.gt.1"){
                    e=d[1:nrow(d) %in% grep("missense",d$consequence) & !is.na(d$MPC) & d$MPC>=1 & d$gene %in% gene.sets[[myset]],]
                  }
                  if(myclass %in% "missense.MPC.1.to.2"){
                    e=d[1:nrow(d) %in% grep("missense",d$consequence) & !is.na(d$MPC) & d$MPC>=1 & d$MPC<2 & d$gene %in% gene.sets[[myset]],]
                  }
                  if(myclass %in% "missense.MPC.gt.2"){
                    e=d[1:nrow(d) %in% grep("missense",d$consequence) & !is.na(d$MPC) & d$MPC>=2 & d$gene %in% gene.sets[[myset]],]
                  }
                  if(myclass %in% "missense.MPC.gt.1.CADDgt25"){
                    e=d[1:nrow(d) %in% grep("missense",d$consequence) & !is.na(d$MPC) & d$MPC>=1 & !is.na(d$CADD) & d$CADD>25 & d$gene %in% gene.sets[[myset]],]
                  }
                  if(myclass %in% "missense.MPC.1.to.2.CADDgt25"){
                    e=d[1:nrow(d) %in% grep("missense",d$consequence) & !is.na(d$MPC) & d$MPC>=1 & d$MPC<2 & !is.na(d$CADD) & d$CADD>25 & d$gene %in% gene.sets[[myset]],]
                  }
                  if(myclass %in% "missense.MPC.gt.2.CADDgt25"){
                    e=d[1:nrow(d) %in% grep("missense",d$consequence) & !is.na(d$MPC) & d$MPC>=2 & !is.na(d$CADD) & d$CADD>25 & d$gene %in% gene.sets[[myset]],]
                  }
                  if(myclass %in% "all.LoF.SNV"){
                    e=d[1:nrow(d) %in% grep("stop_gained|splice_donor|splice_acceptor",d$consequence)  & d$gene %in% gene.sets[[myset]],]
                  }
                  if(myclass %in% "LoF.and.missense.SNV"){
                    e=d[1:nrow(d) %in% grep("stop_gained|splice_donor|splice_acceptor|missense_variant",d$consequence)  & d$gene %in% gene.sets[[myset]],]
                  }

#subset variants within data frame to rarer variants
                  if(maf.cutoff<0.001){
                    e=e[e$DDD.EUR.parents.MAF<maf.cutoff & e$gnomAD.POPMAX.MAF<maf.cutoff,]
                  }
                  #subset variants within data frame e
                  if(remove.gnomad.hemi){
                    e=e[e$total.n.hemi.gnomAD==0,]
                  }
#####subset samples within data frame e
#restrict to trio probands
                  if(trios.only){
                    e=e[e$dad|e$is.trio.proband,]
                  }
#restrict to those with suspected X-linked inheritance
                  if(candidate.xlinked){
                    e=e[e$dad | e$is.suspected.xlinked,]
                  }
# now count the number of variants and the number of individuals with variants in dads and probands
                  n.vars=c(nrow(unique(e[e$dad,c("proband.id","varname")])),nrow(unique(e[e$proband,c("proband.id","varname")])))
                  counts=rbind(counts,n.vars)
                }
              }
            }
          }
    }
  }
}

d=d.indicating.denovos

#add some extra columns on to the variant counts
counts=as.data.frame(counts)
colnames(counts)[1:2]=c("n.vars.dads","n.vars.probands")
counts$class=classes
counts$gene.set=sets
counts$n.dads=n.dads.all
counts$n.probands=my.n.probands

### Having counted variants in father and male probands, let's now calculate rates
counts$rate.vars.per.dad=counts$n.vars.dads/counts$n.dads
counts$rate.vars.per.proband=counts$n.vars.probands/counts$n.probands
### Calculate the ratio of these
counts$ratio.rate.vars=counts$rate.vars.per.proband/counts$rate.vars.per.dad
### Calculate the excess=observed-expected
counts$excess.vars=(counts$n.vars.probands - counts$n.probands * counts$rate.vars.per.dad)
### Calculate the attributable fraction = excess/N
counts$excess.vars.over.N=counts$excess.vars/counts$n.probands
### Label the rows
rownames(counts)=paste(counts$class,counts$gene.set,sep="_")
### Calculate the positive predictive value = (observed-expected)/observed
counts$prop.vars.pathogenic=counts$excess.vars/counts$n.vars.probands
### Calculate a Poisson p-value for the enrichment of variants in probands
counts$pois.p.ut = ppois(counts$n.vars.probands-1,counts$rate.vars.per.dad * counts$n.probands,lower.tail=F)
### Calculate confidence intervals for these various metrics
counts$att.frac.lt=NA
counts$att.frac.ut=NA
counts$burden.lt=NA
counts$burden.ut=NA
counts$ppv.lt=NA
counts$ppv.ut=NA
for(i in 1:nrow(counts)){
  #att.fraction = O/n-E/n
  ### for case/control, calculate CI on difference in poisson rates moverci(obs_child,n_child,obs_dads,n_dads,distrib = ‘poi’,contrast="RD)
    ci.att.frac=moverci(counts[i,"n.vars.probands"],counts[i,"n.probands"],counts[i,"n.vars.dads"],counts[i,"n.dads"],distrib = "poi",contrast="RD")
      counts$att.frac.lt[i] = ci.att.frac[1,"Lower"]
      counts$att.frac.ut[i] = ci.att.frac[1,"Upper"]

    #burden=O/E
    ### for case/control, calculate CI on ratio of poisson rates moverci(obs_child,n_child,obs_dads,n_dads,distrib = ‘poi’,contrast="RR")
      ci.burden=moverci(counts[i,"n.vars.probands"],counts[i,"n.probands"],counts[i,"n.vars.dads"],counts[i,"n.dads"],distrib = "poi",contrast="RR")
      counts$burden.lt[i] = ci.burden[1,"Lower"]
      counts$burden.ut[i] = ci.burden[1,"Upper"]

    #ppv=(O-E)/O = 1 - E/O = 1-(E/n)/(O/n)
    ### for case/control, calculate CI on ratio of poisson rates for (E/n)/(O/n) using moverci(obs_dads,n_dads,obs_child,n_child,distrib = ‘poi’,contrast="RR"), then feed into formula
      ci.exp.over.obs=moverci(counts[i,"n.vars.dads"],counts[i,"n.dads"],counts[i,"n.vars.probands"],counts[i,"n.probands"],distrib = "poi",contrast="RR")
      counts$ppv.lt[i] = 1-ci.exp.over.obs[1,"Upper"]
      counts$ppv.ut[i] = 1-ci.exp.over.obs[1,"Lower"]
  }

write.table(counts,"burden_analysis_on_males.from_case_control_analysis.txt",quote=F,sep="\t",row.names=T)


#### Bootstrap cases and controls in order to get confidence intervals on some key parameters
### load list of trios chosen for bootstrapping in script3.burden_analysis_and_per_gene_tests_on_de_novos.R  
load("trios_selected_for_bootstrapping.RData")

### load list of all IDs of probands and fathers from the case/control analysis
load("all_probands_and_dads_with_randomised_IDs_used_in_released_files.RData")
probands=indivs.to.use.in.case.control.analysis$probands
dads=indivs.to.use.in.case.control.analysis$dads

#### select n.iterations sets of all probands and of dads --> these will be used for the boostrapping
### This bit is slow - if you want to just load the set of individuals that were already prepared, change this to if(FALSE) in order to load individuals_selected_for_bootstrapping_in_case_control_analysis.RData
if(TRUE){
  boot.all.probs=NULL
  boot.dads=NULL
  for(i in 1:n.iterations){
    if(i %/% 50 - i/50 ==0){
      cat("selecting proband and dad set #",i,"\n")
    }
    boot.all.probs=cbind(boot.all.probs,probands[sample(n.probands.all,n.probands.all,replace=T)] )
    boot.dads=cbind(boot.dads,dads[sample(n.dads.all,n.dads.all,replace=T)])
  }
  boot.indivs=list("trio.probands"=bootstrapped.trios[["M"]],"all.probands"=boot.all.probs,"dads"=boot.dads)
  save(boot.indivs,file="individuals_selected_for_bootstrapping_in_case_control_analysis.RData")
} else {
  load("individuals_selected_for_bootstrapping_in_case_control_analysis.RData")
}

#### restrict to variants with 0 hemizygotes in gnomAD
d=d[d$total.n.hemi.gnomAD==0,]

#### split data frame of variants up into list according to individual ID
d.trio.probands=d[d$is.trio.proband,]
d.trio.probands.list=split(d.trio.probands,d.trio.probands$randomised.trio.ID)

d.all.probands=d[!d$dad,]
d.all.probands.list=split(d.all.probands,d.all.probands[,"proband.id"])

d.dads=d[d$dad,]
d.dad.list=split(d.dads,d.dads[,"proband.id"])

### now bootstrap probands and dads and repeat counts of variants
cat("ready to start bootstrapping\n")
boot.counts.list=list("trios.only_FALSE_LoF.and.functional" =NULL, "trios.only_TRUE_LoF.and.functional"  =NULL,"trios.only_TRUE_LoF.and.missense.SNV"=NULL)
#we will do this both with all male probands, or just male trio probands
trios.only.options=c(FALSE,TRUE)
for(trios.only in trios.only.options){
    for(b in 1:n.iterations){
        if(b %/% 50 - b/50 ==0){
            cat(trios.only,"\t",b,"\n")
        }
### assemble dataframe of variants from bootstrapped dads and bootstrapped probands
    these.dads=boot.indivs[["dads"]][,b]
    dads.x=d.dads[unlist(sapply(sort(these.dads),function(x){which(d.dads$proband.id %in% x)})),]
    if(trios.only){
      these.probands=boot.indivs[["trio.probands"]][,b]
      # note that for this trio-only analysis, we are interested only in SNVs because this is for testing Haldane's hypothesis which relies on an assumption about the mutation rate ratio in male versus females,
      # which is only known for SNVs
      my.classes=c("LoF.and.functional","LoF.and.missense.SNV")
      n.probs=n.trio.probands.all
      probs.x=d.trio.probands[unlist(sapply(sort(these.probands),function(x){which(d.trio.probands$randomised.trio.ID %in% x)})),]
      probs.x=probs.x[,!colnames(probs.x) %in%  "person_stable_id"]
    }else {
      these.probands=boot.indivs[["all.probands"]][,b]
      # for this analysis of all trios, we are interested in the overall burden of all LoF+missense/inframe variants
      my.classes=c("LoF.and.functional")
      n.probs=n.probands.all
      probs.x=d.all.probands[unlist(sapply(sort(these.probands),function(x){which(d.all.probands$proband.id %in% x)})),]
    }
    #create a new data.frame of variants from dads and probands
    my.boostrapped.data=as.data.frame(rbind(probs.x,dads.x))
    for(myclass in my.classes){
        #restrict to variants in the right functional categories
        if(myclass %in% c("LoF.and.functional")){
            d=my.boostrapped.data[my.boostrapped.data$class %in% c("LoF","functional"),]
        }   
        if(myclass %in% "LoF.and.missense.SNV"){
            d=my.boostrapped.data[1:nrow(my.boostrapped.data) %in% grep("stop_gained|splice_donor|splice_acceptor|missense_variant",my.boostrapped.data$consequence),]
        }
        #now count variants in either all X-linked genes, known X-linked DD-associated genes, or known XLR DD-associated genes
        count.per.set=c()
        for(myset in  c("all","XL","Hemi.only")){
            #subset variants from data frame d to those in the relevant gene set
            e=d[d$gene %in% gene.sets[[myset]],]
            #count the variants in dads and probands
            n.vars=c(nrow(e[e$dad,c("proband.id","varname")]),nrow(e[e$proband,c("proband.id","varname")]))
            count.per.set=c(count.per.set,n.vars)
        }
    #save the counts for this boostrap
        boot.counts.list[[paste("trios.only",trios.only,myclass,sep="_")]]=rbind( boot.counts.list[[paste("trios.only",trios.only,myclass,sep="_")]],count.per.set)
    }
  }
}

boot.counts.list.new=list()
for(trios.only in trios.only.options){
    if(trios.only){
        my.classes=c("LoF.and.functional","LoF.and.missense.SNV")
    } else {
        my.classes=c("LoF.and.functional")
    }
    for(myclass in my.classes){
        boot.counts=boot.counts.list[[paste("trios.only",trios.only,myclass,sep="_")]]
        colnames(boot.counts)[1:4]=c("n.vars.dads.all","n.vars.probands.all","n.vars.dads.XL","n.vars.probands.XL")
        if(trios.only & myclass=="LoF.and.missense.SNV"){
          colnames(boot.counts)[5:6]=c("n.vars.dads.Hemi.only","n.vars.probands.Hemi.only")
        }
  # calculate expected number in probands in the two different gene sets
        rownames(boot.counts)=NULL
        boot.counts=as.data.frame(boot.counts)
        if(trios.only){
            n.probs=n.trio.probands.all
            if( myclass=="LoF.and.missense.SNV"){
              boot.counts$exp.vars.probands.Hemi.only=n.probs* boot.counts$n.vars.dads.Hemi.only/n.dads.all
              boot.counts$excess.Hemi.only=boot.counts$n.vars.probands.Hemi.only-boot.counts$exp.vars.probands.Hemi.only
            }
        }else {
            n.probs=n.probands.all
        }
        boot.counts$exp.vars.probands.all=n.probs* boot.counts$n.vars.dads.all/n.dads.all
        boot.counts$exp.vars.probands.XL=n.probs* boot.counts$n.vars.dads.XL/n.dads.all
  # calculate the excess in the two different gene sets
        boot.counts$excess.all=boot.counts$n.vars.probands.all-boot.counts$exp.vars.probands.all
        boot.counts$excess.XL=boot.counts$n.vars.probands.XL-boot.counts$exp.vars.probands.XL
  # calculate the fraction of the overall excess that is in known X-linked DDG2P genes
        boot.counts$fraction.excess.known=boot.counts$excess.XL/boot.counts$excess.all
  # save the results of boostrapping in this set of probands
        boot.counts.list.new[[paste("trios.only",trios.only,myclass,sep="_")]]=boot.counts
    }
}
save(boot.counts.list.new,file="results_of_bootstrapping_in_case_control_analysis.all.RData")

