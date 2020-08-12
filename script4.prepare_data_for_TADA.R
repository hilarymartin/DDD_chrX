setwd("/lustre/scratch115/projects/ddd/users/hcm/DDD/chrX_analyses/code_and_data_to_release/")
n.trio.probands.all=5138

#load variants from case/control analysis
load("annotated_variants_for_case_control_analysis.RData")   

#load de novos
load("de_novos_in_both_sexes.post_filtering.0_hemis_gnomAD.MAX_AF_0.001.RData")
x.denovos=d.sep2[d.sep2$sex=="M",]

#load mutation rates for chrX
load("mutation_rates_per_individual.males.per_longest_Ensembl_transcript.adjusted_for_sex_specific_coverage_in_V5.with_approximate_indel_rates.RData")

# create a data frame of variants from the case control analysis, excluding variants that passed the de novo fitlering (since these will be treated separated in TADA)      
# also restrict to variants with 0 hemizygotes in gnomAD, so we focus in on variants most likely to be pathogenic
# remove a few genes that didn't have mutation rates available
d.cc=d[!d$de.novo.passed & d$total.n.hemi.gnomAD==0 & d$gene %in% rownames(adj.rates.denovos.in.males),]

#### Now do the counts per gene of de novos and for the case/control analysis in TADA
# de novo counts, using stringent de novos that passed our filtering optimised to detect de novos
lof.dnm.count=aggregate(as.numeric(is.lof)~symbol,sum,data=x.denovos)
func.dnm.count=aggregate(as.numeric(is.functional)~symbol,sum,data=x.denovos)
rownames(lof.dnm.count)=lof.dnm.count[,1]
rownames(func.dnm.count)=func.dnm.count[,1]
colnames(lof.dnm.count) = c("symbol","count")
colnames(func.dnm.count) = c("symbol","count")

# case counts for TADA
lof.case.count=aggregate(is.lof~gene,sum,data=d.cc[!d.cc$dad,])
func.case.count=aggregate(is.func~gene,sum,data=d.cc[!d.cc$dad,])

# control counts (dads) for TADA
lof.control.count=aggregate(is.lof~gene,sum,data=d.cc[d.cc$dad,])
func.control.count=aggregate(is.func~gene,sum,data=d.cc[d.cc$dad,])

### set up an empty dataframe of TADA input with 0s, which we will then populate with the observe counts
# need to calculate the expected number of de novos for LoF (cls1) and functional (cls2) [mut.cls1 and mut.cls2] by multiplying the rates by the number of trio probands
tada.data=data.frame(gene.id=rownames(adj.rates.denovos.in.males),
  mut.cls1=adj.rates.denovos.in.males$lof * n.trio.probands.all,
  mut.cls2=adj.rates.denovos.in.males$functional * n.trio.probands.all,
  dn.cls1=rep(0,nrow(adj.rates.denovos.in.males)),
  dn.cls2=rep(0,nrow(adj.rates.denovos.in.males)),
  case.cls1=rep(0,nrow(adj.rates.denovos.in.males)),
  case.cls2=rep(0,nrow(adj.rates.denovos.in.males)),
  ctrl.cls1=rep(0,nrow(adj.rates.denovos.in.males)),
  ctrl.cls2=rep(0,nrow(adj.rates.denovos.in.males)))
rownames(tada.data) = tada.data$gene.id

### now populate the TADA dataframe with the observed counts
tada.data[lof.dnm.count[,1],"dn.cls1"] = lof.dnm.count[,2]
tada.data[func.dnm.count[,1],"dn.cls2"] = func.dnm.count[,2]
tada.data[lof.case.count[,1],"case.cls1"] = lof.case.count[,2]
tada.data[func.case.count[,1],"case.cls2"] = func.case.count[,2]
tada.data[lof.control.count[,1],"ctrl.cls1"] = lof.control.count[,2]
tada.data[func.control.count[,1],"ctrl.cls2"] = func.control.count[,2]

save(tada.data,file="input_for_TADA.males.case_control_and_de_novos.RData")
