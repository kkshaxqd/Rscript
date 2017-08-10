rm(list=ls())
source("http://bioconductor.org/biocLite.R") 
biocLite("DESeq")
library(DESeq)
count = read.table("stdin", header=TRUE, row.names=1 )
cond1 = c("control", "control", "control")
cond2 = c("treatment", "treatment", "treatment")
conds = factor(c(cond1, cond2))
cdata = newCountDataSet(count, conds)
esize = estimateSizeFactors(cdata)
edisp = estimateDispersions(esize)
rdata = nbinomTest(edisp, "control", "treatment")
write.table(rdata, file="", sep="\t", row.name=FALSE, quote=FALSE)

java -Xmx15g -jar $gatk \
-T BaseRecalibrator \
-R $ref_dir/hg19.fasta \
-I $bam_dir/$sample.dedupped.realigned.bam \
-knownSites $ref_dir/dbsnp_137.hg19.vcf \
-knownSites $ref_dir/Mills_and_1000G_gold_standard.indels.hg19.vcf \
-knownSites $ref_dir/1000G_phase1.indels.hg19.vcf \
-o $other_dir/$sample.recal.grp \
-plots $other_dir/$sample.recal.grp.pdf
java -Xmx15g -jar $gatk \
-T BaseRecalibrator \
-R $ref_dir/hg19.fasta \
-I $bam_dir/$sample.dedupped.realigned.bam \
-BQSR $other_dir/$sample.recal.grp \
-o $other_dir/$sample.post_recal.grp \
-plots $other_dir/$sample.post_recal.grp.pdf \
-knownSites $ref_dir/dbsnp_137.hg19.vcf \
-knownSites $ref_dir/Mills_and_1000G_gold_standard.indels.hg19.vcf \
-knownSites $ref_dir/1000G_phase1.indels.hg19.vcf
java -Xmx15g -jar $gatk \
-T PrintReads \
-R $ref_dir/hg19.fasta \
-I $bam_dir/$sample.dedupped.realigned.bam \
-BQSR $other_dir/$sample.recal.grp \
-o $bam_dir/$sample.dedupped.realigned.recal.bam

java -Xmx8g -Djava.io.tmpdir=./tmp -jar $path_to_GenomeAnalysisTK/GenomeAnalysisTK.jar  -T BaseRecalibrator -R $path_to_ref_exome/hg19.fa -I $Tumor_bam_realigned -o $samplerecalgrp -known $path_to_ref_exome/dbsnp_138.hg19.vcf  -known $path_to_ref_exome/1000G_phase1.indels.hg19.vcf -known $path_to_ref_exome/Mills_and_1000G_gold_standard.indels.hg19.vcf 
java -Xmx8g -Djava.io.tmpdir=./tmp -jar $path_to_GenomeAnalysisTK/GenomeAnalysisTK.jar  -T BaseRecalibrator -R $path_to_ref_exome/hg19.fa -I $Tumor_bam_realigned -BQSR $samplerecalgrp -o $samplepostrecalgrp -known $path_to_ref_exome/1000G_phase1.indels.hg19.vcf -known $path_to_ref_exome/Mills_and_1000G_gold_standard.indels.hg19.vcf 
java -Xmx8g -Djava.io.tmpdir=./tmp -jar $path_to_GenomeAnalysisTK/GenomeAnalysisTK.jar  -T PrintReads -R $path_to_ref_exome/hg19.fa -I $Tumor_bam_realigned -BQSR $samplerecalgrp -o $Tumor_bam_realigned_recal



