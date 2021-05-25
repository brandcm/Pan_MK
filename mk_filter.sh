#!/bin/bash

#SBATCH --partition=short
#SBATCH --job-name=mk_filter
#SBATCH --output=mk_filter.out
#SBATCH --error=mk_filter.err
#SBATCH --time=12:00:00
#SBATCH --mem=50G
#SBATCH --cpus-per-task=2
#SBATCH --account=ting

# load first module
# bcftools and anaconda do not play nicely together here so I load and unload modules to work around this
module load bcftools

# normalize and split multiallelics and combine indel/variant rows into single records from VCF
# chromosome 18 throws an error so note that there are fewer flags. This is okay because no variants have to be realigned for this chromosome.
bcftools norm -m +any -f pantro6_autosomes.fa pantro6.chr1.gatk.called.raw.vcf.gz -O z -o norm_pantro6.chr1.gatk.called.raw.vcf.gz
bcftools norm -m +any -f pantro6_autosomes.fa pantro6.chr2A.gatk.called.raw.vcf.gz -O z -o norm_pantro6.chr2A.gatk.called.raw.vcf.gz
bcftools norm -m +any -f pantro6_autosomes.fa pantro6.chr2B.gatk.called.raw.vcf.gz -O z -o norm_pantro6.chr2B.gatk.called.raw.vcf.gz
bcftools norm -m +any -f pantro6_autosomes.fa pantro6.chr3.gatk.called.raw.vcf.gz -O z -o norm_pantro6.chr3.gatk.called.raw.vcf.gz
bcftools norm -m +any -f pantro6_autosomes.fa pantro6.chr4.gatk.called.raw.vcf.gz -O z -o norm_pantro6.chr4.gatk.called.raw.vcf.gz
bcftools norm -m +any -f pantro6_autosomes.fa pantro6.chr5.gatk.called.raw.vcf.gz -O z -o norm_pantro6.chr5.gatk.called.raw.vcf.gz
bcftools norm -m +any -f pantro6_autosomes.fa pantro6.chr6.gatk.called.raw.vcf.gz -O z -o norm_pantro6.chr6.gatk.called.raw.vcf.gz
bcftools norm -m +any -f pantro6_autosomes.fa pantro6.chr7.gatk.called.raw.vcf.gz -O z -o norm_pantro6.chr7.gatk.called.raw.vcf.gz
bcftools norm -m +any -f pantro6_autosomes.fa pantro6.chr8.gatk.called.raw.vcf.gz -O z -o norm_pantro6.chr8.gatk.called.raw.vcf.gz
bcftools norm -m +any -f pantro6_autosomes.fa pantro6.chr9.gatk.called.raw.vcf.gz -O z -o norm_pantro6.chr9.gatk.called.raw.vcf.gz
bcftools norm -m +any -f pantro6_autosomes.fa pantro6.chr10.gatk.called.raw.vcf.gz -O z -o norm_pantro6.chr10.gatk.called.raw.vcf.gz
bcftools norm -m +any -f pantro6_autosomes.fa pantro6.chr11.gatk.called.raw.vcf.gz -O z -o norm_pantro6.chr11.gatk.called.raw.vcf.gz
bcftools norm -m +any -f pantro6_autosomes.fa pantro6.chr12.gatk.called.raw.vcf.gz -O z -o norm_pantro6.chr12.gatk.called.raw.vcf.gz
bcftools norm -m +any -f pantro6_autosomes.fa pantro6.chr13.gatk.called.raw.vcf.gz -O z -o norm_pantro6.chr13.gatk.called.raw.vcf.gz
bcftools norm -m +any -f pantro6_autosomes.fa pantro6.chr14.gatk.called.raw.vcf.gz -O z -o norm_pantro6.chr14.gatk.called.raw.vcf.gz
bcftools norm -m +any -f pantro6_autosomes.fa pantro6.chr15.gatk.called.raw.vcf.gz -O z -o norm_pantro6.chr15.gatk.called.raw.vcf.gz
bcftools norm -m +any -f pantro6_autosomes.fa pantro6.chr16.gatk.called.raw.vcf.gz -O z -o norm_pantro6.chr16.gatk.called.raw.vcf.gz
bcftools norm -m +any -f pantro6_autosomes.fa pantro6.chr17.gatk.called.raw.vcf.gz -O z -o norm_pantro6.chr17.gatk.called.raw.vcf.gz
bcftools norm -f pantro6_autosomes.fa pantro6.chr18.gatk.called.raw.vcf.gz -O z -o norm_pantro6.chr18.gatk.called.raw.vcf.gz
bcftools norm -m +any -f pantro6_autosomes.fa pantro6.chr19.gatk.called.raw.vcf.gz -O z -o norm_pantro6.chr19.gatk.called.raw.vcf.gz
bcftools norm -m +any -f pantro6_autosomes.fa pantro6.chr20.gatk.called.raw.vcf.gz -O z -o norm_pantro6.chr20.gatk.called.raw.vcf.gz
bcftools norm -m +any -f pantro6_autosomes.fa pantro6.chr21.gatk.called.raw.vcf.gz -O z -o norm_pantro6.chr21.gatk.called.raw.vcf.gz
bcftools norm -m +any -f pantro6_autosomes.fa pantro6.chr22.gatk.called.raw.vcf.gz -O z -o norm_pantro6.chr22.gatk.called.raw.vcf.gz

module unload bcftools

# create index for new VCFs

module load anaconda3
for f in norm*.raw.vcf.gz; do tabix -p vcf $f; done
module unload anaconda3

# filter by region and position/genotype quality, note the different max depth filters per chromosome

module load bcftools

bcftools filter -R CDS_autosomes.bed norm_pantro6.chr1.gatk.called.raw.vcf.gz | bcftools view -m2 -M2 -v snps | bcftools filter -g 5 -i 'QUAL>=30 & QD>2 & MQ>=30' | bcftools filter -S . -i 'FMT/DP>=10 & FMT/GQ>=30' | bcftools filter -i 'AN>=112 & AC>0 & AN!=AC' | bcftools filter -e 'FMT/DP > 92.4' -O z -o mk_filtered.pantro6.chr1.gatk.called.raw.vcf.gz
bcftools filter -R CDS_autosomes.bed norm_pantro6.chr2A.gatk.called.raw.vcf.gz | bcftools view -m2 -M2 -v snps | bcftools filter -g 5 -i 'QUAL>=30 & QD>2 & MQ>=30' | bcftools filter -S . -i 'FMT/DP>=10 & FMT/GQ>=30' | bcftools filter -i 'AN>=112 & AC>0 & AN!=AC' | bcftools filter -e 'FMT/DP > 94.6' -O z -o mk_filtered.pantro6.chr2A.gatk.called.raw.vcf.gz
bcftools filter -R CDS_autosomes.bed norm_pantro6.chr2B.gatk.called.raw.vcf.gz | bcftools view -m2 -M2 -v snps | bcftools filter -g 5 -i 'QUAL>=30 & QD>2 & MQ>=30' | bcftools filter -S . -i 'FMT/DP>=10 & FMT/GQ>=30' | bcftools filter -i 'AN>=112 & AC>0 & AN!=AC' | bcftools filter -e 'FMT/DP > 99.3' -O z -o mk_filtered.pantro6.chr2B.gatk.called.raw.vcf.gz
bcftools filter -R CDS_autosomes.bed norm_pantro6.chr3.gatk.called.raw.vcf.gz | bcftools view -m2 -M2 -v snps | bcftools filter -g 5 -i 'QUAL>=30 & QD>2 & MQ>=30' | bcftools filter -S . -i 'FMT/DP>=10 & FMT/GQ>=30' | bcftools filter -i 'AN>=112 & AC>0 & AN!=AC' | bcftools filter -e 'FMT/DP > 98.3' -O z -o mk_filtered.pantro6.chr3.gatk.called.raw.vcf.gz
bcftools filter -R CDS_autosomes.bed norm_pantro6.chr4.gatk.called.raw.vcf.gz | bcftools view -m2 -M2 -v snps | bcftools filter -g 5 -i 'QUAL>=30 & QD>2 & MQ>=30' | bcftools filter -S . -i 'FMT/DP>=10 & FMT/GQ>=30' | bcftools filter -i 'AN>=112 & AC>0 & AN!=AC' | bcftools filter -e 'FMT/DP > 99.9' -O z -o mk_filtered.pantro6.chr4.gatk.called.raw.vcf.gz
bcftools filter -R CDS_autosomes.bed norm_pantro6.chr5.gatk.called.raw.vcf.gz | bcftools view -m2 -M2 -v snps | bcftools filter -g 5 -i 'QUAL>=30 & QD>2 & MQ>=30' | bcftools filter -S . -i 'FMT/DP>=10 & FMT/GQ>=30' | bcftools filter -i 'AN>=112 & AC>0 & AN!=AC' | bcftools filter -e 'FMT/DP > 97.8' -O z -o mk_filtered.pantro6.chr5.gatk.called.raw.vcf.gz
bcftools filter -R CDS_autosomes.bed norm_pantro6.chr6.gatk.called.raw.vcf.gz | bcftools view -m2 -M2 -v snps | bcftools filter -g 5 -i 'QUAL>=30 & QD>2 & MQ>=30' | bcftools filter -S . -i 'FMT/DP>=10 & FMT/GQ>=30' | bcftools filter -i 'AN>=112 & AC>0 & AN!=AC' | bcftools filter -e 'FMT/DP > 98.8' -O z -o mk_filtered.pantro6.chr6.gatk.called.raw.vcf.gz
bcftools filter -R CDS_autosomes.bed norm_pantro6.chr7.gatk.called.raw.vcf.gz | bcftools view -m2 -M2 -v snps | bcftools filter -g 5 -i 'QUAL>=30 & QD>2 & MQ>=30' | bcftools filter -S . -i 'FMT/DP>=10 & FMT/GQ>=30' | bcftools filter -i 'AN>=112 & AC>0 & AN!=AC' | bcftools filter -e 'FMT/DP > 92.6' -O z -o mk_filtered.pantro6.chr7.gatk.called.raw.vcf.gz
bcftools filter -R CDS_autosomes.bed norm_pantro6.chr8.gatk.called.raw.vcf.gz | bcftools view -m2 -M2 -v snps | bcftools filter -g 5 -i 'QUAL>=30 & QD>2 & MQ>=30' | bcftools filter -S . -i 'FMT/DP>=10 & FMT/GQ>=30' | bcftools filter -i 'AN>=112 & AC>0 & AN!=AC' | bcftools filter -e 'FMT/DP > 95.3' -O z -o mk_filtered.pantro6.chr8.gatk.called.raw.vcf.gz
bcftools filter -R CDS_autosomes.bed norm_pantro6.chr9.gatk.called.raw.vcf.gz | bcftools view -m2 -M2 -v snps | bcftools filter -g 5 -i 'QUAL>=30 & QD>2 & MQ>=30' | bcftools filter -S . -i 'FMT/DP>=10 & FMT/GQ>=30' | bcftools filter -i 'AN>=112 & AC>0 & AN!=AC' | bcftools filter -e 'FMT/DP > 94.1' -O z -o mk_filtered.pantro6.chr9.gatk.called.raw.vcf.gz
bcftools filter -R CDS_autosomes.bed norm_pantro6.chr10.gatk.called.raw.vcf.gz | bcftools view -m2 -M2 -v snps | bcftools filter -g 5 -i 'QUAL>=30 & QD>2 & MQ>=30' | bcftools filter -S . -i 'FMT/DP>=10 & FMT/GQ>=30' | bcftools filter -i 'AN>=112 & AC>0 & AN!=AC' | bcftools filter -e 'FMT/DP > 92.1' -O z -o mk_filtered.pantro6.chr10.gatk.called.raw.vcf.gz
bcftools filter -R CDS_autosomes.bed norm_pantro6.chr11.gatk.called.raw.vcf.gz | bcftools view -m2 -M2 -v snps | bcftools filter -g 5 -i 'QUAL>=30 & QD>2 & MQ>=30' | bcftools filter -S . -i 'FMT/DP>=10 & FMT/GQ>=30' | bcftools filter -i 'AN>=112 & AC>0 & AN!=AC' | bcftools filter -e 'FMT/DP > 93.3' -O z -o mk_filtered.pantro6.chr11.gatk.called.raw.vcf.gz
bcftools filter -R CDS_autosomes.bed norm_pantro6.chr12.gatk.called.raw.vcf.gz | bcftools view -m2 -M2 -v snps | bcftools filter -g 5 -i 'QUAL>=30 & QD>2 & MQ>=30' | bcftools filter -S . -i 'FMT/DP>=10 & FMT/GQ>=30' | bcftools filter -i 'AN>=112 & AC>0 & AN!=AC' | bcftools filter -e 'FMT/DP > 94.7' -O z -o mk_filtered.pantro6.chr12.gatk.called.raw.vcf.gz
bcftools filter -R CDS_autosomes.bed norm_pantro6.chr13.gatk.called.raw.vcf.gz | bcftools view -m2 -M2 -v snps | bcftools filter -g 5 -i 'QUAL>=30 & QD>2 & MQ>=30' | bcftools filter -S . -i 'FMT/DP>=10 & FMT/GQ>=30' | bcftools filter -i 'AN>=112 & AC>0 & AN!=AC' | bcftools filter -e 'FMT/DP > 100.1' -O z -o mk_filtered.pantro6.chr13.gatk.called.raw.vcf.gz
bcftools filter -R CDS_autosomes.bed norm_pantro6.chr14.gatk.called.raw.vcf.gz | bcftools view -m2 -M2 -v snps | bcftools filter -g 5 -i 'QUAL>=30 & QD>2 & MQ>=30' | bcftools filter -S . -i 'FMT/DP>=10 & FMT/GQ>=30' | bcftools filter -i 'AN>=112 & AC>0 & AN!=AC' | bcftools filter -e 'FMT/DP > 94.9' -O z -o mk_filtered.pantro6.chr14.gatk.called.raw.vcf.gz
bcftools filter -R CDS_autosomes.bed norm_pantro6.chr15.gatk.called.raw.vcf.gz | bcftools view -m2 -M2 -v snps | bcftools filter -g 5 -i 'QUAL>=30 & QD>2 & MQ>=30' | bcftools filter -S . -i 'FMT/DP>=10 & FMT/GQ>=30' | bcftools filter -i 'AN>=112 & AC>0 & AN!=AC' | bcftools filter -e 'FMT/DP > 92.1' -O z -o mk_filtered.pantro6.chr15.gatk.called.raw.vcf.gz
bcftools filter -R CDS_autosomes.bed norm_pantro6.chr16.gatk.called.raw.vcf.gz | bcftools view -m2 -M2 -v snps | bcftools filter -g 5 -i 'QUAL>=30 & QD>2 & MQ>=30' | bcftools filter -S . -i 'FMT/DP>=10 & FMT/GQ>=30' | bcftools filter -i 'AN>=112 & AC>0 & AN!=AC' | bcftools filter -e 'FMT/DP > 90.7' -O z -o mk_filtered.pantro6.chr16.gatk.called.raw.vcf.gz
bcftools filter -R CDS_autosomes.bed norm_pantro6.chr17.gatk.called.raw.vcf.gz | bcftools view -m2 -M2 -v snps | bcftools filter -g 5 -i 'QUAL>=30 & QD>2 & MQ>=30' | bcftools filter -S . -i 'FMT/DP>=10 & FMT/GQ>=30' | bcftools filter -i 'AN>=112 & AC>0 & AN!=AC' | bcftools filter -e 'FMT/DP > 94.5' -O z -o mk_filtered.pantro6.chr17.gatk.called.raw.vcf.gz
bcftools filter -R CDS_autosomes.bed norm_pantro6.chr18.gatk.called.raw.vcf.gz | bcftools view -m2 -M2 -v snps | bcftools filter -g 5 -i 'QUAL>=30 & QD>2 & MQ>=30' | bcftools filter -S . -i 'FMT/DP>=10 & FMT/GQ>=30' | bcftools filter -i 'AN>=112 & AC>0 & AN!=AC' | bcftools filter -e 'FMT/DP > 97.0' -O z -o mk_filtered.pantro6.chr18.gatk.called.raw.vcf.gz
bcftools filter -R CDS_autosomes.bed norm_pantro6.chr19.gatk.called.raw.vcf.gz | bcftools view -m2 -M2 -v snps | bcftools filter -g 5 -i 'QUAL>=30 & QD>2 & MQ>=30' | bcftools filter -S . -i 'FMT/DP>=10 & FMT/GQ>=30' | bcftools filter -i 'AN>=112 & AC>0 & AN!=AC' | bcftools filter -e 'FMT/DP > 96.0' -O z -o mk_filtered.pantro6.chr19.gatk.called.raw.vcf.gz
bcftools filter -R CDS_autosomes.bed norm_pantro6.chr20.gatk.called.raw.vcf.gz | bcftools view -m2 -M2 -v snps | bcftools filter -g 5 -i 'QUAL>=30 & QD>2 & MQ>=30' | bcftools filter -S . -i 'FMT/DP>=10 & FMT/GQ>=30' | bcftools filter -i 'AN>=112 & AC>0 & AN!=AC' | bcftools filter -e 'FMT/DP > 99.6' -O z -o mk_filtered.pantro6.chr20.gatk.called.raw.vcf.gz
bcftools filter -R CDS_autosomes.bed norm_pantro6.chr21.gatk.called.raw.vcf.gz | bcftools view -m2 -M2 -v snps | bcftools filter -g 5 -i 'QUAL>=30 & QD>2 & MQ>=30' | bcftools filter -S . -i 'FMT/DP>=10 & FMT/GQ>=30' | bcftools filter -i 'AN>=112 & AC>0 & AN!=AC' | bcftools filter -e 'FMT/DP > 94.1' -O z -o mk_filtered.pantro6.chr21.gatk.called.raw.vcf.gz
bcftools filter -R CDS_autosomes.bed norm_pantro6.chr22.gatk.called.raw.vcf.gz | bcftools view -m2 -M2 -v snps | bcftools filter -g 5 -i 'QUAL>=30 & QD>2 & MQ>=30' | bcftools filter -S . -i 'FMT/DP>=10 & FMT/GQ>=30' | bcftools filter -i 'AN>=112 & AC>0 & AN!=AC' | bcftools filter -e 'FMT/DP > 98.7' -O z -o mk_filtered.pantro6.chr22.gatk.called.raw.vcf.gz

# remove temp files creates in stages 1 & 2 above
rm norm*.gz
rm norm*.gz.tbi