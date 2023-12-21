# Sort and Index

# Sort the Tumor.bam file and output the sorted BAM file as Tumor.sorted.bam
samtools sort Tumor.bam > Tumor.sorted.bam

# Index the sorted Tumor BAM file
samtools index Tumor.sorted.bam

# Sort the Control.bam file and output the sorted BAM file as Control.sorted.bam
samtools sort Control.bam > Control.sorted.bam

# Index the sorted Control BAM file
samtools index Control.sorted.bam

# Realignment (for both)
#### Control ####

# Step 1: RealignerTargetCreator - Identify intervals for realignment in the control sample
java -jar ../../../tools/GenomeAnalysisTK.jar -T RealignerTargetCreator \
     -R ../../../annotations/human_g1k_v37.fasta \
     -I ../../01_Samtools/Control.sorted.bam \
     -o realigner.intervals \
     -L ../../../annotations/Captured_Regions.bed

# Step 2: IndelRealigner - Perform indel realignment on the control sample using the generated target intervals
java -jar ../../../tools/GenomeAnalysisTK.jar -T IndelRealigner \
     -R ../../../annotations/human_g1k_v37.fasta \
     -I ../../01_Samtools/Control.sorted.bam \
     -targetIntervals realigner.intervals \
     -o Control.sorted.realigned.bam \
     -L ../../../annotations/Captured_Regions.bed

#### Tumor ####

# Step 1: RealignerTargetCreator - Identify intervals for realignment in the tumor sample
java -jar ../../../tools/GenomeAnalysisTK.jar -T RealignerTargetCreator \
     -R ../../../annotations/human_g1k_v37.fasta \
     -I ../../01_Samtools/Tumor.sorted.bam \
     -o realigner.intervals \
     -L ../../../annotations/Captured_Regions.bed

# Step 2: IndelRealigner - Perform indel realignment on the tumor sample using the generated target intervals
java -jar ../../../tools/GenomeAnalysisTK.jar -T IndelRealigner \
     -R ../../../annotations/human_g1k_v37.fasta \
     -I ../../01_Samtools/Tumor.sorted.bam \
     -targetIntervals realigner.intervals \
     -o Tumor.sorted.realigned.bam \
     -L ../../../annotations/Captured_Regions.bed

# Recalibration (for both)
#### Tumor ####

# Recalibration
## Step 1: BaseRecalibrator - Perform base quality score recalibration
java -jar ../../../tools/GenomeAnalysisTK.jar -T BaseRecalibrator \
     -R ../../../annotations/human_g1k_v37.fasta \
     -I ../../02_Realignment/Tumor/Tumor.sorted.realigned.bam \
     -knownSites ../../../annotations/hapmap_3.3.b37.vcf \
     -o recal.table \
     -L ../../../annotations/Captured_Regions.bed

## Step 2: BaseRecalibrator - Perform base quality score recalibration using the recalibration table from step 1
java -jar ../../../tools/GenomeAnalysisTK.jar -T BaseRecalibrator \
     -R ../../../annotations/human_g1k_v37.fasta \
     -I ../../02_Realignment/Tumor/Tumor.sorted.realigned.bam \
     -knownSites ../../../annotations/hapmap_3.3.b37.vcf \
     -BQSR recal.table \
     -o after_recal.table \
     -L ../../../annotations/Captured_Regions.bed

## Step 3: PrintReads - Apply base quality score recalibration to the aligned reads
java -jar ../../../tools/GenomeAnalysisTK.jar -T PrintReads \
     -R ../../../annotations/human_g1k_v37.fasta \
     -I ../../02_Realignment/Tumor/Tumor.sorted.realigned.bam \
     -BQSR recal.table \
     -o Tumor.sorted.realigned.recalibrated.bam \
     -L ../../../annotations/Captured_Regions.bed \
     --emit_original_quals

## Step 4: AnalyzeCovariates - Generate analysis of base quality score recalibration covariates
java -jar ../../../tools/GenomeAnalysisTK.jar -T AnalyzeCovariates \
     -R ../../../annotations/human_g1k_v37.fasta \
     -before recal.table \
     -after after_recal.table \
     -csv recal.csv \
     -plots recal.pdf

## Step 5: Count lines with original base quality in the recalibrated BAM file using samtools and grep
samtools view Tumor.sorted.realigned.recalibrated.bam | grep OQ | wc -l

(SAME FOR CONTROL, JUST CHANGE THE NAME OF THE INPUT FILES)

# Deduplicates (for both)
samtools flagstat Tumor.sorted.realigned.recalibrated.bam

samtools view Tumor.sorted.realigned.recalibrated.bam | wc -l

java -jar ../../../Tools/picard.jar MarkDuplicates I=../../03_Recalibration/Tumor/Tumor.sorted.realigned.recalibrated.bam \
     O=Tumor.sorted.markdup.bam \
     REMOVE_DUPLICATES=false \
     TMP_DIR=/tmp \
     METRICS_FILE= Tumor.picard.log \
     ASSUME_SORTED=true

samtools flagstat Tumor.sorted.markdup.bam

java -jar ../../../Tools/picard.jar MarkDuplicates I=Tumor.sorted.realigned.recalibrated.bam \
     O=Tumor.sorted.dedup.bam \
     REMOVE_DUPLICATES=true \
     TMP_DIR=/tmp \
     METRICS_FILE= Tumor.picard.log \
     ASSUME_SORTED=true

samtools index Tumor.sorted.dedup.bam

samtools view Tumor.sorted.dedup.bam | wc -l

samtools flagstat Tumor.sorted.dedup.bam > Tumor.dedup.flag.txt

# Somatic variant calling
# SNP calling on Control, use Varscan
samtools mpileup -B -f ../../Annotations/human_g1k_v37.fasta ../../Control.sorted.dedup.bam > Control.sorted.pileup

# mpileup2snp generate vcf file that contains SNP of given bam file
java -jar ../../Tools/VarScan.v2.3.9.jar mpileup2snp Control.sorted.pileup --p-value 0.01 --output-vcf 1 > Control.VARSCAN.vcf

# Somatic variant calling on control and tumor
samtools mpileup -q 1 -f ../../Annotations/human_g1k_v37.fasta Tumor.sorted.dedup.bam > Tumor.sorted.pileup

java -jar ../../Tools/VarScan.v2.3.9.jar somatic Control.sorted.pileup Tumor.sorted.pileup \
     --output-snp somatic.pm \
     --output-indel somatic.indel \
     --output-vcf 1

# Comparison of control SNPs with somatic events found in the tumor
# Are they all SNVs or are some only SNPs with different AF (for structural events)?

# Somatic events annotation
## Filtering at --min-meanDP 30 to keep more significant ones
vcftools --min-meanDP 30 --remove-indels --vcf somatic.pm.vcf --out somatic.pm --recode --recode-INFO-all

java -Xmx6g -jar ../../tools/snpEff/snpEff.jar -v hg19kg somatic.pm.vcf -s somatic.pm.vcf.html > somatic.pm.ann1.vcf

java -Xmx2g -jar ../../tools/snpEff/SnpSift.jar Annotate ../../annotations/hapmap_3.3.b37.vcf somatic.pm.ann1.vcf > somatic.pm.ann2.hapmap.vcf

java -Xmx6g -jar ../../tools/snpEff/SnpSift.jar Annotate ../../annotations/clinvar_Pathogenic.vcf somatic.pm.ann2.hapmap.vcf > somatic.pm.ann3.clinvar.vcf

cat somatic.pm.ann3.clinvar.vcf | java -Xmx4g -jar ../../tools/snpEff/SnpSift.jar filter "(ANN[ANY].IMPACT = 'HIGH')" >Filter_IMPACT.txt

cat somatic.pm.ann3.clinvar.vcf | java -Xmx4g -jar ../../tools/snpEff/SnpSift.jar filter "(exists CLNSIG)" >Filter_CLNSIG.txt

# Copy number variations call
samtools mpileup -q 1 -f ../../../annotations/human_g1k_v37.fasta \
    Control.sorted.dedup.bam Tumor.sorted.dedup.bam | \
    java -jar ../../Tools/VarScan.v2.3.9.jar copynumber --output-file SCNA --mpileup 1

# output SCNA.copynumber
java -jar ../../Tools/VarScan.v2.3.9.jar copyCaller SCNA.copynumber --output-file SCNA.copynumber.called

# Now use SCNA.copynumber.called in the CBS.R script (already ready for use)
# Obtain segmentation plots (SegPlot.pdf)

# Purity Ploidy estimation

# Take called heterozygous SNPs in normal
grep -E "(^#|HET=1)" ../08_Somatic variant calling/Control.VARSCAN.vcf > Control.het.vcf

# Compute read counts for alleles on het SNPs, both in tumor and control
java -jar ../Tools/GenomeAnalysisTK.jar -T ASEReadCounter -R ../Annotations/human_g1k_v37.fasta \
    -o Control.csv -I Control.sorted.dedup.bam -sites Control.het.vcf -U ALLOW_N_CIGAR_READS \
    -minDepth 20 --minMappingQuality 20 --minBaseQuality 20

java -jar ../Tools/GenomeAnalysisTK.jar -T ASEReadCounter -R ../Annotations/human_g1k_v37.fasta \
    -o Tumor.csv -I Tumor.sorted.dedup.bam -sites Control.het.vcf -U ALLOW_N_CIGAR_READS \
    -minDepth 20 --minMappingQuality 20 --minBaseQuality 20

# Somatic.indel (VarScan somatic output, used with the CSV in CLONET.R script)
# Use CLONET_project.R script (Clonet and TPES)
