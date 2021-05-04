# Build a database (or use an existing one)
java -jar snpEff/snpEff.jar build -gff3 -v chimp

# Annotate VCFs
for f in *.vcf.gz; do java -Xmx4g -jar snpEff/snpEff.jar -v chimp snpEff/data/chimp/$f > $f.annotated

# Rename VCFs
for file in *.annotated; do mv "$file" "${file%vcf.gz.annotated}annotated.vcf"; done

# Generate allele frequency files 
for f in *.annotated.vcf; do vcftools/bin/vcftools --vcf $f --keep hss_ids.txt --freq --out hss.$f; done
for f in *.annotated.vcf; do vcftools/bin/vcftools --vcf $f --keep ppn_ids.txt --freq --out ppn.$f; done
for f in *.annotated.vcf; do vcftools/bin/vcftools --vcf $f --keep pte_ids.txt --freq --out pte.$f; done
for f in *.annotated.vcf; do vcftools/bin/vcftools --vcf $f --keep pts_ids.txt --freq --out pts.$f; done
for f in *.annotated.vcf; do vcftools/bin/vcftools --vcf $f --keep ptt_ids.txt --freq --out ptt.$f; done
for f in *.annotated.vcf; do vcftools/bin/vcftools --vcf $f --keep ptv_ids.txt --freq --out ptv.$f; done

# Run the analysis. Note that this will run a MK test per gene, a MK test with a MAF filter (0.1) per gene, and a MK test per exon.
while read f; do Rscript command_line_mk_script.R ann_vcfs/$f.annotated.vcf ppn.$f.annotated.vcf.frq hss.$f.annotated.vcf.frq $f ppn_$f.csv ppn_mk_$f.csv ppn_mk_maf_$f.csv ppn_exon_$f.csv; done < chr_names.txt
while read f; do Rscript command_line_mk_script.R ann_vcfs/$f.annotated.vcf pte.$f.annotated.vcf.frq hss.$f.annotated.vcf.frq $f pte_$f.csv pte_mk_$f.csv pte_mk_maf_$f.csv pte_exon_$f.csv; done < chr_names.txt
while read f; do Rscript command_line_mk_script.R ann_vcfs/$f.annotated.vcf pts.$f.annotated.vcf.frq hss.$f.annotated.vcf.frq $f pts_$f.csv pts_mk_$f.csv pts_mk_maf_$f.csv pts_exon_$f.csv; done < chr_names.txt
while read f; do Rscript command_line_mk_script.R ann_vcfs/$f.annotated.vcf ptt.$f.annotated.vcf.frq hss.$f.annotated.vcf.frq $f ptt_$f.csv ptt_mk_$f.csv ptt_mk_maf_$f.csv ptt_exon_$f.csv; done < chr_names.txt
while read f; do Rscript command_line_mk_script.R ann_vcfs/$f.annotated.vcf ptv.$f.annotated.vcf.frq hss.$f.annotated.vcf.frq $f ptv_$f.csv ptv_mk_$f.csv ptv_mk_maf_$f.csv ptv_exon_$f.csv; done < chr_names.txt

# Combine files
awk 'FNR==1 && NR!=1{next;}{print}' ppn_mk_chr*.csv > ppn_mk.csv
awk 'FNR==1 && NR!=1{next;}{print}' pte_mk_chr*.csv > pte_mk.csv
awk 'FNR==1 && NR!=1{next;}{print}' pts_mk_chr*.csv > pts_mk.csv
awk 'FNR==1 && NR!=1{next;}{print}' ptt_mk_chr*.csv > ptt_mk.csv
awk 'FNR==1 && NR!=1{next;}{print}' ptv_mk_chr*.csv > ptv_mk.csv

awk 'FNR==1 && NR!=1{next;}{print}' ppn_mk_maf_chr*.csv > ppn_mk_maf.csv
awk 'FNR==1 && NR!=1{next;}{print}' pte_mk_maf_chr*.csv > pte_mk_maf.csv
awk 'FNR==1 && NR!=1{next;}{print}' pts_mk_maf_chr*.csv > pts_mk_maf.csv
awk 'FNR==1 && NR!=1{next;}{print}' ptt_mk_maf_chr*.csv > ptt_mk_maf.csv
awk 'FNR==1 && NR!=1{next;}{print}' ptv_mk_maf_chr*.csv > ptv_mk_maf.csv

awk 'FNR==1 && NR!=1{next;}{print}' ppn_exon_chr*.csv > ppn_mk_exon.csv
awk 'FNR==1 && NR!=1{next;}{print}' pte_exon_chr*.csv > pte_mk_exon.csv
awk 'FNR==1 && NR!=1{next;}{print}' pts_exon_chr*.csv > pts_mk_exon.csv
awk 'FNR==1 && NR!=1{next;}{print}' ptt_exon_chr*.csv > ptt_mk_exon.csv
awk 'FNR==1 && NR!=1{next;}{print}' ptv_exon_chr*.csv > ptv_mk_exon.csv