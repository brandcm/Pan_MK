#!/bin/bash

#SBATCH --partition=short	### queue to submit to (short, fat, or long--short is typical usage, fat is for big jobs, *I think* long is for jobs that take more than 24 hours)
#SBATCH --job-name=real_mk_filter    ### job name, whatever you want it to be
#SBATCH --output=real_mk_filter.out   ### file in which to store job stdout (make sure to change this if you want to resubmit the script but save the last job's .out file because it will just overw$
#SBATCH --error=real_mk_filter.err    ### file in which to store job stderr (errors will be printed in your .err file)
#SBATCH --time=12:00:00               ### wall-clock time limit, e.g., 5 minutes
#SBATCH --mem=50G              ### memory limit, per cpu, in MB
#SBATCH --cpus-per-task=2	### number of cores for each task (for a lot of smaller jobs you can leave this as 1; for genome alignment, which is a 200MB job, I did 24 CPUs/task on the$
#SBATCH --account=ting

#load any packages your script needs to use
module load bcftools
#module load anaconda3

#put your code below

# normalize and split multiallelics and combine indel/variant rows into single records from VCF
#bcftools norm -m +any -f pantro6_autosomes.fa pantro6.chr1.gatk.called.raw.vcf.gz -O z -o norm_pantro6.chr1.gatk.called.raw.vcf.gz
#bcftools norm -m +any -f pantro6_autosomes.fa pantro6.chr2A.gatk.called.raw.vcf.gz -O z -o norm_pantro6.chr2A.gatk.called.raw.vcf.gz
#bcftools norm -m +any -f pantro6_autosomes.fa pantro6.chr2B.gatk.called.raw.vcf.gz -O z -o norm_pantro6.chr2B.gatk.called.raw.vcf.gz
#bcftools norm -m +any -f pantro6_autosomes.fa pantro6.chr3.gatk.called.raw.vcf.gz -O z -o norm_pantro6.chr3.gatk.called.raw.vcf.gz
#bcftools norm -m +any -f pantro6_autosomes.fa pantro6.chr4.gatk.called.raw.vcf.gz -O z -o norm_pantro6.chr4.gatk.called.raw.vcf.gz
#bcftools norm -m +any -f pantro6_autosomes.fa pantro6.chr5.gatk.called.raw.vcf.gz -O z -o norm_pantro6.chr5.gatk.called.raw.vcf.gz
#bcftools norm -m +any -f pantro6_autosomes.fa pantro6.chr6.gatk.called.raw.vcf.gz -O z -o norm_pantro6.chr6.gatk.called.raw.vcf.gz
#bcftools norm -m +any -f pantro6_autosomes.fa pantro6.chr7.gatk.called.raw.vcf.gz -O z -o norm_pantro6.chr7.gatk.called.raw.vcf.gz
#bcftools norm -m +any -f pantro6_autosomes.fa pantro6.chr8.gatk.called.raw.vcf.gz -O z -o norm_pantro6.chr8.gatk.called.raw.vcf.gz
#bcftools norm -m +any -f pantro6_autosomes.fa pantro6.chr9.gatk.called.raw.vcf.gz -O z -o norm_pantro6.chr9.gatk.called.raw.vcf.gz
#bcftools norm -m +any -f pantro6_autosomes.fa pantro6.chr10.gatk.called.raw.vcf.gz -O z -o norm_pantro6.chr10.gatk.called.raw.vcf.gz
#bcftools norm -m +any -f pantro6_autosomes.fa pantro6.chr11.gatk.called.raw.vcf.gz -O z -o norm_pantro6.chr11.gatk.called.raw.vcf.gz
#bcftools norm -m +any -f pantro6_autosomes.fa pantro6.chr12.gatk.called.raw.vcf.gz -O z -o norm_pantro6.chr12.gatk.called.raw.vcf.gz
#bcftools norm -m +any -f pantro6_autosomes.fa pantro6.chr13.gatk.called.raw.vcf.gz -O z -o norm_pantro6.chr13.gatk.called.raw.vcf.gz
#bcftools norm -m +any -f pantro6_autosomes.fa pantro6.chr14.gatk.called.raw.vcf.gz -O z -o norm_pantro6.chr14.gatk.called.raw.vcf.gz
#bcftools norm -m +any -f pantro6_autosomes.fa pantro6.chr15.gatk.called.raw.vcf.gz -O z -o norm_pantro6.chr15.gatk.called.raw.vcf.gz
#bcftools norm -m +any -f pantro6_autosomes.fa pantro6.chr16.gatk.called.raw.vcf.gz -O z -o norm_pantro6.chr16.gatk.called.raw.vcf.gz
#bcftools norm -m +any -f pantro6_autosomes.fa pantro6.chr17.gatk.called.raw.vcf.gz -O z -o norm_pantro6.chr17.gatk.called.raw.vcf.gz
bcftools norm -m +any -f pantro6_autosomes.fa pantro6.chr18.gatk.called.raw.vcf.gz -O z -o norm_pantro6.chr18.gatk.called.raw.vcf.gz
#bcftools norm -m +any -f pantro6_autosomes.fa pantro6.chr19.gatk.called.raw.vcf.gz -O z -o norm_pantro6.chr19.gatk.called.raw.vcf.gz
#bcftools norm -m +any -f pantro6_autosomes.fa pantro6.chr20.gatk.called.raw.vcf.gz -O z -o norm_pantro6.chr20.gatk.called.raw.vcf.gz
#bcftools norm -m +any -f pantro6_autosomes.fa pantro6.chr21.gatk.called.raw.vcf.gz -O z -o norm_pantro6.chr21.gatk.called.raw.vcf.gz
#bcftools norm -m +any -f pantro6_autosomes.fa pantro6.chr22.gatk.called.raw.vcf.gz -O z -o norm_pantro6.chr22.gatk.called.raw.vcf.gz

# create index for new VCFs
#tabix -p vcf norm_pantro6.chr1.gatk.called.raw.vcf.gz
#tabix -p vcf norm_pantro6.chr2A.gatk.called.raw.vcf.gz
#tabix -p vcf norm_pantro6.chr2B.gatk.called.raw.vcf.gz
#tabix -p vcf norm_pantro6.chr3.gatk.called.raw.vcf.gz
#tabix -p vcf norm_pantro6.chr4.gatk.called.raw.vcf.gz
#tabix -p vcf norm_pantro6.chr5.gatk.called.raw.vcf.gz
#tabix -p vcf norm_pantro6.chr6.gatk.called.raw.vcf.gz
#tabix -p vcf norm_pantro6.chr7.gatk.called.raw.vcf.gz
#tabix -p vcf norm_pantro6.chr8.gatk.called.raw.vcf.gz
#tabix -p vcf norm_pantro6.chr9.gatk.called.raw.vcf.gz
#tabix -p vcf norm_pantro6.chr10.gatk.called.raw.vcf.gz
#tabix -p vcf norm_pantro6.chr11.gatk.called.raw.vcf.gz
#tabix -p vcf norm_pantro6.chr12.gatk.called.raw.vcf.gz
#tabix -p vcf norm_pantro6.chr13.gatk.called.raw.vcf.gz
#tabix -p vcf norm_pantro6.chr14.gatk.called.raw.vcf.gz
#tabix -p vcf norm_pantro6.chr15.gatk.called.raw.vcf.gz
#tabix -p vcf norm_pantro6.chr16.gatk.called.raw.vcf.gz
#tabix -p vcf norm_pantro6.chr17.gatk.called.raw.vcf.gz
#tabix -p vcf norm_pantro6.chr18.gatk.called.raw.vcf.gz
#tabix -p vcf norm_pantro6.chr19.gatk.called.raw.vcf.gz
#tabix -p vcf norm_pantro6.chr20.gatk.called.raw.vcf.gz
#tabix -p vcf norm_pantro6.chr21.gatk.called.raw.vcf.gz
#tabix -p vcf norm_pantro6.chr22.gatk.called.raw.vcf.gz

# filter by region and position/genotype quality
#bcftools filter -R CDS_autosomes.bed norm_pantro6.chr1.gatk.called.raw.vcf.gz | bcftools view -m2 -M2 -v snps | bcftools filter -g 5 -i 'QUAL>=30 & QD>2 & MQ>=30' | bcftools filter -S . -i 'FMT/DP>=10 & FMT/GQ>=30' | bcftools filter -i 'AN>=112 & AC>0 & AN!=AC' | bcftools filter -e 'FMT/DP > 92.4' -O z -o mk_filtered.pantro6.chr1.gatk.called.raw.vcf.gz
#bcftools filter -R CDS_autosomes.bed norm_pantro6.chr2A.gatk.called.raw.vcf.gz | bcftools view -m2 -M2 -v snps | bcftools filter -g 5 -i 'QUAL>=30 & QD>2 & MQ>=30' | bcftools filter -S . -i 'FMT/DP>=10 & FMT/GQ>=30' | bcftools filter -i 'AN>=112 & AC>0 & AN!=AC' | bcftools filter -e 'FMT/DP > 94.6' -O z -o mk_filtered.pantro6.chr2A.gatk.called.raw.vcf.gz
#bcftools filter -R CDS_autosomes.bed norm_pantro6.chr2B.gatk.called.raw.vcf.gz | bcftools view -m2 -M2 -v snps | bcftools filter -g 5 -i 'QUAL>=30 & QD>2 & MQ>=30' | bcftools filter -S . -i 'FMT/DP>=10 & FMT/GQ>=30' | bcftools filter -i 'AN>=112 & AC>0 & AN!=AC' | bcftools filter -e 'FMT/DP > 99.3' -O z -o mk_filtered.pantro6.chr2B.gatk.called.raw.vcf.gz
#bcftools filter -R CDS_autosomes.bed norm_pantro6.chr3.gatk.called.raw.vcf.gz | bcftools view -m2 -M2 -v snps | bcftools filter -g 5 -i 'QUAL>=30 & QD>2 & MQ>=30' | bcftools filter -S . -i 'FMT/DP>=10 & FMT/GQ>=30' | bcftools filter -i 'AN>=112 & AC>0 & AN!=AC' | bcftools filter -e 'FMT/DP > 98.3' -O z -o mk_filtered.pantro6.chr3.gatk.called.raw.vcf.gz
#bcftools filter -R CDS_autosomes.bed norm_pantro6.chr4.gatk.called.raw.vcf.gz | bcftools view -m2 -M2 -v snps | bcftools filter -g 5 -i 'QUAL>=30 & QD>2 & MQ>=30' | bcftools filter -S . -i 'FMT/DP>=10 & FMT/GQ>=30' | bcftools filter -i 'AN>=112 & AC>0 & AN!=AC' | bcftools filter -e 'FMT/DP > 99.9' -O z -o mk_filtered.pantro6.chr4.gatk.called.raw.vcf.gz
#bcftools filter -R CDS_autosomes.bed norm_pantro6.chr5.gatk.called.raw.vcf.gz | bcftools view -m2 -M2 -v snps | bcftools filter -g 5 -i 'QUAL>=30 & QD>2 & MQ>=30' | bcftools filter -S . -i 'FMT/DP>=10 & FMT/GQ>=30' | bcftools filter -i 'AN>=112 & AC>0 & AN!=AC' | bcftools filter -e 'FMT/DP > 97.8' -O z -o mk_filtered.pantro6.chr5.gatk.called.raw.vcf.gz
#bcftools filter -R CDS_autosomes.bed norm_pantro6.chr6.gatk.called.raw.vcf.gz | bcftools view -m2 -M2 -v snps | bcftools filter -g 5 -i 'QUAL>=30 & QD>2 & MQ>=30' | bcftools filter -S . -i 'FMT/DP>=10 & FMT/GQ>=30' | bcftools filter -i 'AN>=112 & AC>0 & AN!=AC' | bcftools filter -e 'FMT/DP > 98.8' -O z -o mk_filtered.pantro6.chr6.gatk.called.raw.vcf.gz
#bcftools filter -R CDS_autosomes.bed norm_pantro6.chr7.gatk.called.raw.vcf.gz | bcftools view -m2 -M2 -v snps | bcftools filter -g 5 -i 'QUAL>=30 & QD>2 & MQ>=30' | bcftools filter -S . -i 'FMT/DP>=10 & FMT/GQ>=30' | bcftools filter -i 'AN>=112 & AC>0 & AN!=AC' | bcftools filter -e 'FMT/DP > 92.6' -O z -o mk_filtered.pantro6.chr7.gatk.called.raw.vcf.gz
#bcftools filter -R CDS_autosomes.bed norm_pantro6.chr8.gatk.called.raw.vcf.gz | bcftools view -m2 -M2 -v snps | bcftools filter -g 5 -i 'QUAL>=30 & QD>2 & MQ>=30' | bcftools filter -S . -i 'FMT/DP>=10 & FMT/GQ>=30' | bcftools filter -i 'AN>=112 & AC>0 & AN!=AC' | bcftools filter -e 'FMT/DP > 95.3' -O z -o mk_filtered.pantro6.chr8.gatk.called.raw.vcf.gz
#bcftools filter -R CDS_autosomes.bed norm_pantro6.chr9.gatk.called.raw.vcf.gz | bcftools view -m2 -M2 -v snps | bcftools filter -g 5 -i 'QUAL>=30 & QD>2 & MQ>=30' | bcftools filter -S . -i 'FMT/DP>=10 & FMT/GQ>=30' | bcftools filter -i 'AN>=112 & AC>0 & AN!=AC' | bcftools filter -e 'FMT/DP > 94.1' -O z -o mk_filtered.pantro6.chr9.gatk.called.raw.vcf.gz
#bcftools filter -R CDS_autosomes.bed norm_pantro6.chr10.gatk.called.raw.vcf.gz | bcftools view -m2 -M2 -v snps | bcftools filter -g 5 -i 'QUAL>=30 & QD>2 & MQ>=30' | bcftools filter -S . -i 'FMT/DP>=10 & FMT/GQ>=30' | bcftools filter -i 'AN>=112 & AC>0 & AN!=AC' | bcftools filter -e 'FMT/DP > 92.1' -O z -o mk_filtered.pantro6.chr10.gatk.called.raw.vcf.gz
#bcftools filter -R CDS_autosomes.bed norm_pantro6.chr11.gatk.called.raw.vcf.gz | bcftools view -m2 -M2 -v snps | bcftools filter -g 5 -i 'QUAL>=30 & QD>2 & MQ>=30' | bcftools filter -S . -i 'FMT/DP>=10 & FMT/GQ>=30' | bcftools filter -i 'AN>=112 & AC>0 & AN!=AC' | bcftools filter -e 'FMT/DP > 93.3' -O z -o mk_filtered.pantro6.chr11.gatk.called.raw.vcf.gz
#bcftools filter -R CDS_autosomes.bed norm_pantro6.chr12.gatk.called.raw.vcf.gz | bcftools view -m2 -M2 -v snps | bcftools filter -g 5 -i 'QUAL>=30 & QD>2 & MQ>=30' | bcftools filter -S . -i 'FMT/DP>=10 & FMT/GQ>=30' | bcftools filter -i 'AN>=112 & AC>0 & AN!=AC' | bcftools filter -e 'FMT/DP > 94.7' -O z -o mk_filtered.pantro6.chr12.gatk.called.raw.vcf.gz
#bcftools filter -R CDS_autosomes.bed norm_pantro6.chr13.gatk.called.raw.vcf.gz | bcftools view -m2 -M2 -v snps | bcftools filter -g 5 -i 'QUAL>=30 & QD>2 & MQ>=30' | bcftools filter -S . -i 'FMT/DP>=10 & FMT/GQ>=30' | bcftools filter -i 'AN>=112 & AC>0 & AN!=AC' | bcftools filter -e 'FMT/DP > 100.1' -O z -o mk_filtered.pantro6.chr13.gatk.called.raw.vcf.gz
#bcftools filter -R CDS_autosomes.bed norm_pantro6.chr14.gatk.called.raw.vcf.gz | bcftools view -m2 -M2 -v snps | bcftools filter -g 5 -i 'QUAL>=30 & QD>2 & MQ>=30' | bcftools filter -S . -i 'FMT/DP>=10 & FMT/GQ>=30' | bcftools filter -i 'AN>=112 & AC>0 & AN!=AC' | bcftools filter -e 'FMT/DP > 94.9' -O z -o mk_filtered.pantro6.chr14.gatk.called.raw.vcf.gz
#bcftools filter -R CDS_autosomes.bed norm_pantro6.chr15.gatk.called.raw.vcf.gz | bcftools view -m2 -M2 -v snps | bcftools filter -g 5 -i 'QUAL>=30 & QD>2 & MQ>=30' | bcftools filter -S . -i 'FMT/DP>=10 & FMT/GQ>=30' | bcftools filter -i 'AN>=112 & AC>0 & AN!=AC' | bcftools filter -e 'FMT/DP > 92.1' -O z -o mk_filtered.pantro6.chr15.gatk.called.raw.vcf.gz
#bcftools filter -R CDS_autosomes.bed norm_pantro6.chr16.gatk.called.raw.vcf.gz | bcftools view -m2 -M2 -v snps | bcftools filter -g 5 -i 'QUAL>=30 & QD>2 & MQ>=30' | bcftools filter -S . -i 'FMT/DP>=10 & FMT/GQ>=30' | bcftools filter -i 'AN>=112 & AC>0 & AN!=AC' | bcftools filter -e 'FMT/DP > 90.7' -O z -o mk_filtered.pantro6.chr16.gatk.called.raw.vcf.gz
#bcftools filter -R CDS_autosomes.bed norm_pantro6.chr17.gatk.called.raw.vcf.gz | bcftools view -m2 -M2 -v snps | bcftools filter -g 5 -i 'QUAL>=30 & QD>2 & MQ>=30' | bcftools filter -S . -i 'FMT/DP>=10 & FMT/GQ>=30' | bcftools filter -i 'AN>=112 & AC>0 & AN!=AC' | bcftools filter -e 'FMT/DP > 94.5' -O z -o mk_filtered.pantro6.chr17.gatk.called.raw.vcf.gz
#bcftools filter -R CDS_autosomes.bed norm_pantro6.chr18.gatk.called.raw.vcf.gz | bcftools view -m2 -M2 -v snps | bcftools filter -g 5 -i 'QUAL>=30 & QD>2 & MQ>=30' | bcftools filter -S . -i 'FMT/DP>=10 & FMT/GQ>=30' | bcftools filter -i 'AN>=112 & AC>0 & AN!=AC' | bcftools filter -e 'FMT/DP > 97.0' -O z -o mk_filtered.pantro6.chr18.gatk.called.raw.vcf.gz
#bcftools filter -R CDS_autosomes.bed norm_pantro6.chr19.gatk.called.raw.vcf.gz | bcftools view -m2 -M2 -v snps | bcftools filter -g 5 -i 'QUAL>=30 & QD>2 & MQ>=30' | bcftools filter -S . -i 'FMT/DP>=10 & FMT/GQ>=30' | bcftools filter -i 'AN>=112 & AC>0 & AN!=AC' | bcftools filter -e 'FMT/DP > 96.0' -O z -o mk_filtered.pantro6.chr19.gatk.called.raw.vcf.gz
#bcftools filter -R CDS_autosomes.bed norm_pantro6.chr20.gatk.called.raw.vcf.gz | bcftools view -m2 -M2 -v snps | bcftools filter -g 5 -i 'QUAL>=30 & QD>2 & MQ>=30' | bcftools filter -S . -i 'FMT/DP>=10 & FMT/GQ>=30' | bcftools filter -i 'AN>=112 & AC>0 & AN!=AC' | bcftools filter -e 'FMT/DP > 99.6' -O z -o mk_filtered.pantro6.chr20.gatk.called.raw.vcf.gz
#bcftools filter -R CDS_autosomes.bed norm_pantro6.chr21.gatk.called.raw.vcf.gz | bcftools view -m2 -M2 -v snps | bcftools filter -g 5 -i 'QUAL>=30 & QD>2 & MQ>=30' | bcftools filter -S . -i 'FMT/DP>=10 & FMT/GQ>=30' | bcftools filter -i 'AN>=112 & AC>0 & AN!=AC' | bcftools filter -e 'FMT/DP > 94.1' -O z -o mk_filtered.pantro6.chr21.gatk.called.raw.vcf.gz
#bcftools filter -R CDS_autosomes.bed norm_pantro6.chr22.gatk.called.raw.vcf.gz | bcftools view -m2 -M2 -v snps | bcftools filter -g 5 -i 'QUAL>=30 & QD>2 & MQ>=30' | bcftools filter -S . -i 'FMT/DP>=10 & FMT/GQ>=30' | bcftools filter -i 'AN>=112 & AC>0 & AN!=AC' | bcftools filter -e 'FMT/DP > 98.7' -O z -o mk_filtered.pantro6.chr22.gatk.called.raw.vcf.gz

# remove temp files creates in stages 1 & 2 above
#rm norm*.gz
#rm norm*.gz.tbi