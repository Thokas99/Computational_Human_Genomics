library(data.table)
library(CLONETv2)
library(TPES)
library(ggplot2)


control = fread("Control.csv",data.table=F)
control$af = control$altCount/control$totalCount
tumor = fread("Tumor.csv",data.table=F)
tumor$af = tumor$altCount/tumor$totalCount

pileup.control = control[,c(1,2,4,5,14,8)]
colnames(pileup.control) = c("chr","pos","ref","alt","af","cov")

pileup.tumor = tumor[,c(1,2,4,5,14,8)]
colnames(pileup.tumor) = c("chr","pos","ref","alt","af","cov")

#segmentation data from previous lesson
#segmentation info and copy number for each region (mean segment -> log2R)
seg.tb <- fread("SCNA.copynumber.called.seg",data.table=F)

#use CLONET now
bt <- compute_beta_table(seg.tb, pileup.tumor, pileup.control)

## Compute ploidy table with default parameters
pl.table <- compute_ploidy(bt)
#ploidy 1.84, so probably diploid
#admixture (1-purity) and confidence interval
adm.table <- compute_dna_admixture(beta_table = bt, ploidy_table = pl.table)

#allele specific copy number events
#corrected log2R
allele_specific_cna_table <- compute_allele_specific_scna_table(beta_table = bt,
                                                                ploidy_table = pl.table, 
                                                                admixture_table = adm.table)


check.plot <- check_ploidy_and_admixture(beta_table = bt, ploidy_table = pl.table,
                                         admixture_table = adm.table)
print(check.plot)
#each point is a segment, beta and logR computed based on all segment snp. Red dots are points 
#in which i expect specific event given that purity and ploidy

p<-ggplot(data=allele_specific_cna_table,aes(x=cnA,y=cnB,col=log2.corr)) +
  coord_fixed() +
  geom_point()


# TPES

#take only somatic events
snv.reads = fread("somatic.pm",data.table=F)
snv.reads = snv.reads[which(snv.reads$somatic_status=="Somatic"),]
snv.reads = snv.reads[,c("chrom","position","position","tumor_reads1","tumor_reads2")]
colnames(snv.reads) = c("chr","start","end","ref.count","alt.count")
snv.reads$sample = "Sample.1"

#RMB is reference mapping bias
TPES_purity(ID = "Sample.1", SEGfile = seg.tb,
            SNVsReadCountsFile = snv.reads,
            ploidy = pl.table,
            RMB = 0.47, maxAF = 0.6, minCov = 10, minAltReads = 10, minSNVs = 1)

TPES_report(ID = "Sample.1", SEGfile = seg.tb,
            SNVsReadCountsFile = snv.reads,
            ploidy = pl.table,
            RMB = 0.47, maxAF = 0.6, minCov = 10, minAltReads = 10, minSNVs = 1)



