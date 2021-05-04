#!/usr/bin/env Rscript
# Line above is for running via the command line or in a SLURM script

# This script runs per contig per lineage so we'll use arguments to make running this a little easier
args <- commandArgs(TRUE)

# List of command arguments 
# args[1] - path to VCF file
# args[2] - path to target population VCFtools allele frequency file (.frq)
# args[3] - path to outgroup population VCFtools allele frequency file (.frq)
# args[4] - target contig/chromosome string (e.g., "chr1", *Must* match what is listed in VCF file)
# args[5] - filename for CSV of variants analyzed
# args[6] - filename for MK results in a CSV
# args[7] - filename for MK with MAF results in a CSV
# args[8] - filename for MK by exon results in a CSV

# Create dataframe from GFF to get gene names (Official Gene Symbol) and the number of exons. This step assumes there is a GFF that only contains only CDS rows and one transcript per gene. See "single_transcript_per_gene_filtered_gff.R" to build a GFF that retains the longest transcript. 

# load libraries
library(tidyverse)
library(stringr)

gff_file <- read.delim("pantro6_one_transcript_per_gene.gff", header=FALSE, sep="\t")
colnames(gff_file) <- c("A","B","C","D","E","F","G","H","I")

gff_file$D <- as.numeric(as.character(gff_file$D))
gff_file$E <- as.numeric(as.character(gff_file$E))
gff_file <- gff_file[order(gff_file$A, gff_file$D),] #order things nicely (i.e., by position as transcripts in GFFs are rarely ordered this way and we need SNVs to be in order downstream) 

gff_file$ID <- str_extract(gff_file$I, "ID=(.*?);")
gff_file$Gene <- str_extract(gff_file$I, "gene=(.*?);")

gff_file$Gene <- gsub("gene=", "", gff_file$Gene)
gff_file$Gene <- gsub(";", "", gff_file$Gene)

gff_file$exon_no <- ave(gff_file$Gene, gff_file$Gene, FUN=seq_along)
gff_file$exon <- paste(gff_file$Gene,"_exon_", gff_file$exon_no)

# add column with calculated exon length for downstream analyses

gff_file$exon_length <- (gff_file$E - gff_file$D) + 1


genes <- gff_file[c('ID','Gene')]
genes$ID <- gsub("ID=", "", genes$ID)
genes$ID <- gsub(";", "", genes$ID)
genes$Gene <- gsub("gene=", "", genes$Gene)
genes$Gene <- gsub(";", "", genes$Gene)

# Remove duplicate rows from dictionary

#genes$ID<-factor(genes$ID)
#genes$Gene<-factor(genes$Gene)
#genes <- unique(genes)

# Read in VCF file and rename column 1 to "contig"
vcf_path=args[1]
vcf_file<-readLines(vcf_path)
vcf_SNVs<-data.frame(vcf_file[grep(pattern="#CHROM",vcf_file):length(vcf_file)])
vcf_SNVs <- data.frame(do.call('rbind', strsplit(as.character(vcf_SNVs[,1]),'\t',fixed=TRUE)))
colnames(vcf_SNVs) <- as.character(unlist(vcf_SNVs[1,]))
vcf_SNVs<-vcf_SNVs[-1,]
colnames(vcf_SNVs)[1]="chrom"

# Get rid of columns we don't need
vcf_SNVs <- within(vcf_SNVs, rm(ID, QUAL, FILTER, FORMAT)) 
vcf_SNVs <- within(vcf_SNVs, rm('Akwaya-Jean', 'Alfred', 'Alice', 'Andromeda', 'Annie', 'Athanga', 'Banyo', 'Basho', 'Berta', 'Bihati', 'Blanquita', 'Bono', 'Bosco', 'Brigitta', 'Bwamble', 'Catherine', 'Chipita', 'Cindy_schwein', 'Cindy_troglodytes', 'Cindy_verus', 'Clara', 'Cleo', 'Clint', 'Coco_chimp', 'Damian', 'Desmond', 'Doris', 'Dzeeta', 'Frederike', 'Gamin', 'Harriet', 'Hermien', 'HG00513', 'Hortense', 'Ikuru', 'Jimmie', 'Julie_A959', 'Julie_LWC21', 'Kidongo', 'Koby', 'Kombote', 'Kopongo', 'Kosana', 'Koto', 'Kumbuka', 'Lara', 'LB502', 'Linda', 'Luky', 'Marlin', 'Maya', 'Mgbadolite', 'Mike', 'Mirinda', 'Nakuu', 'Natalie', 'Negrita', 'Noemie', 'Padda', 'Paquita', 'Salonga', 'SeppToni', 'Taweh', 'Tibe', 'Tobi', 'Tongo', 'Trixie', 'Ula', 'Vaillant', 'Vincent', 'Washu', 'Yogui'))

# Read in allele frequencies for target Pan lineage
# Allele frequencies here are assumed to be in a .frq format file output from VCFtools. Allele frequency data in other formats will require revision of the code below. 
pan_freq_path=args[2]
pan_freq<-read.delim(pan_freq_path, header=FALSE, skip=1)
colnames(pan_freq)<-c("chrom","pos","n_alleles","n_chromosomes","ref","alt")
pan_freq <- pan_freq[grep('^[chr]', pan_freq$chrom),]
pan_freq <- pan_freq[c(1:6)]
pan_freq<-separate(pan_freq,ref,c("ref_allele","ref_freq"), ":")
pan_freq<-separate(pan_freq,alt,c("alt_allele","alt_freq"), ":")
pan_freq$ref_freq<-as.numeric(pan_freq$ref_freq)
pan_freq$alt_freq<-as.numeric(pan_freq$alt_freq)

# Read in a human genome as the outgroup
hss_path=args[3]
hss_freq<-read.delim(hss_path, header=FALSE, skip=1)
colnames(hss_freq)<-c("chrom","pos","n_alleles","n_chromosomes","ref","alt")
hss_freq <- hss_freq[grep('^[chr]', hss_freq$chrom),]
hss_freq <- hss_freq[c(1:6)]
hss_freq<-separate(hss_freq,ref,c("ref_allele","ref_freq"), ":")
hss_freq<-separate(hss_freq,alt,c("alt_allele","alt_freq"), ":")
hss_freq$ref_freq<-as.numeric(hss_freq$ref_freq)
hss_freq$alt_freq<-as.numeric(hss_freq$alt_freq)

# Attach allele frequencies to SNV file
vcf_SNVs$pan_ref_freq<-pan_freq$ref_freq
vcf_SNVs$pan_alt_freq<-pan_freq$alt_freq
vcf_SNVs$hss_ref_freq<-hss_freq$ref_freq
vcf_SNVs$hss_alt_freq<-hss_freq$alt_freq

# Only keep SNVs that are either synonymous or nonsynonymous
vcf_SNVs<-vcf_SNVs[grep(pattern="synonymous_variant|missense_variant",vcf_SNVs$INFO),]

suppressWarnings(vcf_SNVs$annotation<-separate(vcf_SNVs,INFO,c("A","B","C","D","E","F","G","H","I","J","K","L","M","N","O","P","Q","R"), "\\|")[,"B"]) # Random number of columns here. Warnings are suppressed because the number of split columns will vary depending on the number of annotation results for each SNV. We only need the second column as we are focusing on the first annotation result so all but the first few columns matter (although we will borrow information from the fourth, see below).

# Extract the gene name to a new column. We will run an MK test per gene so we need to identify which SNV belongs to which gene

suppressWarnings(vcf_SNVs$gene<-separate(vcf_SNVs,INFO,c("A","B","C","D","E","F","G","H","I","J","K","L","M","N","O","P","Q","R"), "\\|")[,"D"]) # Warnings suppressed here again. See above for explanation.
vcf_SNVs$gene_name<-NA
vcf_SNVs$gene_name<-genes$Gene[match(vcf_SNVs$gene, genes$ID)]

# Slap the exon number on that data frame as well. Same as above. We will run an MK test per exon (per gene) so we add an identifier. 

vcf_SNVs$POS <- as.numeric(as.character(vcf_SNVs$POS))
vcf_SNVs$exon <- NA

subset_gff_file <- subset(gff_file, gff_file$A == args[4])

vcf_SNVs$exon <- sapply(vcf_SNVs$POS, function(v.element) subset_gff_file[v.element >= subset_gff_file$D & v.element <= subset_gff_file$E,"exon"])

vcf_SNVs$chrom<-factor(vcf_SNVs$chrom) # reset factor level
vcf_SNVs$gene_name<-factor(vcf_SNVs$gene_name) # reset factor level
vcf_SNVs$exon<-as.factor(as.character(vcf_SNVs$exon)) # reset factor level

# I used bcftools -R to subset the CDS loci from a larger VCF file when filtering for this analysis. However, overlap (and possibly improper sorting) of the BED file used to subset loci can result in unsorted and duplicate SNVs; see bcftools documentation - http://samtools.github.io/bcftools/bcftools.html. This could have a large impact on an MK test by inflating the number of SNVs. Here we will cull any duplicates that have made it this far. As this script is run per contig, we only need to check using the "POS" field.  

print("The number of SNVs prior to filtering for duplicates is:", quote = FALSE)
print(nrow(vcf_SNVs))

vcf_SNVs<-vcf_SNVs[!duplicated(vcf_SNVs$POS), ]

# Count the number of annotations per variant. Given that SNVs may have multiple annotations, it will be useful to assess the distribution of annotation counts per contig or summed across the entire genome via a histogram to consider how this variation impacts downstream analyses

library(stringr)
vcf_SNVs$ann_count<-str_count(vcf_SNVs$INFO, ',')
vcf_SNVs$ann_count<-vcf_SNVs$ann_count + 1

# Save the vcf, from which we can generate a histogram of annotation counts per SNV later

write.csv(vcf_SNVs, args[5])

# Count the number of non-duplicated SNVs in our non-MAF MK analyses
print("The number of non-duplicated SNVs:", quote = FALSE)
print(nrow(vcf_SNVs))

#===============================================================================

# Load necessary information for downstream data clean up here. Specifically, we need a list of genes to throw out because they have poor coverage in our total sample and the lengths of genes. We already have exon lengths per the code above. Note that these are not command line arguments because these will not change across contigs/populations and there are already enough args. 

exclude <- read.delim("genes_to_exclude.txt", header = TRUE, sep = "\t")
gene_length <- read.delim("gene_length.txt", header = TRUE, sep = "\t")

#===============================================================================

# Now do the actual analysis

# Create columns with 0s
vcf_SNVs$pN<-0
vcf_SNVs$pS<-0
vcf_SNVs$dN<-0
vcf_SNVs$dS<-0

# Exclude rows of missing data 
vcf_SNVs <- na.exclude(vcf_SNVs)

# Count the number of SNVs in our non-MAF MK analyses
print("The number of SNVs included in our non-MAF MK analyses:", quote = FALSE)
print(nrow(vcf_SNVs))

# Check the category per SNV
for(i in 1:nrow(vcf_SNVs)) {
  if((vcf_SNVs$pan_ref_freq[i]==0) && (vcf_SNVs$hss_alt_freq[i]==0) && "synonymous_variant" %in% vcf_SNVs$annotation[i]) { vcf_SNVs$dS[i]<-1 }
  if((vcf_SNVs$pan_alt_freq[i]==0) && (vcf_SNVs$hss_ref_freq[i]==0) && "synonymous_variant" %in% vcf_SNVs$annotation[i]) { vcf_SNVs$dS[i]<-1 }
  if((vcf_SNVs$pan_ref_freq[i]==0) && (vcf_SNVs$hss_alt_freq[i]==0) && "missense_variant" %in% vcf_SNVs$annotation[i]) { vcf_SNVs$dN[i]<-1 }
  if((vcf_SNVs$pan_alt_freq[i]==0) && (vcf_SNVs$hss_ref_freq[i]==0) && "missense_variant" %in% vcf_SNVs$annotation[i]) { vcf_SNVs$dN[i]<-1 }
  if((vcf_SNVs$pan_ref_freq[i]!=0 && vcf_SNVs$pan_alt_freq[i]!=0) && "synonymous_variant" %in% vcf_SNVs$annotation[i]) { vcf_SNVs$pS[i]<-1 }
  if((vcf_SNVs$pan_ref_freq[i]!=0 && vcf_SNVs$pan_alt_freq[i]!=0) && "missense_variant" %in% vcf_SNVs$annotation[i]) { vcf_SNVs$pN[i]<-1 }
}

# Make empty data frame
MKT<-data.frame(chrom=factor(), gene=factor(), gene_length=numeric(), pN=numeric(), pS=numeric(), dN=numeric(), dS=numeric())

# Reset factor level (no idea why but the script works when this line is here)
vcf_SNVs$gene_name<-factor(vcf_SNVs$gene_name)

# Sum pNs, pSs, dNs, and dSs for each gene across SNVs
for(i in 1:length(levels(vcf_SNVs$gene_name))) {
  temp<-vcf_SNVs[which(vcf_SNVs$gene_name==levels(vcf_SNVs$gene_name)[i]),]
  MKT<-rbind(MKT,data.frame(
    chrom=as.character(temp$chrom[1]),
    gene=as.character(temp$gene_name[1]),
    pN=sum(temp[,"pN"]),
    pS=sum(temp[,"pS"]),
    dN=sum(temp[,"dN"]),
    dS=sum(temp[,"dS"])))
}

# Disable scientific notation in p-value column
options(scipen = 999)

# Calculate some stuff
MKT$pN_pS <- MKT$pN+MKT$pS
MKT$pN_dN <- MKT$pN+MKT$dN
MKT$dS_dN <- MKT$dS+MKT$dN
MKT$dS_pS <- MKT$dS+MKT$pS
MKT$gene_length <- gene_length$length[match(MKT$gene, gene_length$gene)]
MKT$SNVs_per_bp <- (MKT$pN+MKT$pS+MKT$dN+MKT$dS)/MKT$gene_length

# Perform MK test
MKT$fisher.test.P<-99  # create new column for p-values

for(i in 1:nrow(MKT)){
  MKT$fisher.test.P[i]<-fisher.test(matrix(as.numeric(MKT[i,c(6,4,5,3)]), ncol=2))$p.value
  if((MKT$pN[i] == 0 && MKT$dN[i] == 0) || (MKT$pS[i] == 0 && MKT$dS[i] == 0) || (MKT$pS[i] == 0 && MKT$pN[i] == 0) || (MKT$dS[i] == 0 && MKT$dN[i] == 0)) { MKT$fisher.test.P[i]<-NA} # this lines assigns an NA to all p-values that are meaningless, because the contingency table was incomplete
}

# Division from last chunk creates 'Inf's, replace with NAs
is.na(MKT) <- sapply(MKT, is.infinite)

# Calculate neutrality index (Rand and Kann 1996), replace 'Inf's again
MKT$neut_index=(MKT$pN*MKT$dS)/(MKT$pS*MKT$dN)
is.na(MKT$neut_index) <- sapply(MKT$neut_index, is.infinite)

# Calculate direction of selection (Stoletizki and Eyre-Walker 2011)
MKT$dos=((MKT$dN/(MKT$dN+MKT$dS))-(MKT$pN/(MKT$pN+MKT$pS)))

# Remove stuff
MKT <- MKT[!is.na(MKT$fisher.test.P), ]
MKT <- filter(MKT, !(gene %in% exclude$gene))

# Save MK table
write.csv(MKT, args[6])

### MAF ===================================================================

# Redo above analysis using a MAF of 0.1 for Pan
MKT_MAF<-data.frame(chrom=factor(), gene=factor(), pN=numeric(), pS=numeric(), dN=numeric(), dS=numeric())

vcf_SNVs_MAF <- subset(vcf_SNVs, pan_ref_freq == 1 | pan_alt_freq == 1 && pan_ref_freq >= 0.1 | pan_alt_freq >= 0.1)

# Count the number of SNVs in our MAF MK analyses
print("The number of SNVs included in our MAF MK analyses:", quote = FALSE)
print(nrow(vcf_SNVs_MAF))

for(i in 1:length(levels(vcf_SNVs_MAF$gene_name))) {
  temp<-vcf_SNVs_MAF[which(vcf_SNVs_MAF$gene_name==levels(vcf_SNVs_MAF$gene_name)[i]),]
  MKT_MAF<-rbind(MKT_MAF,data.frame(
    chrom=as.character(temp$chrom[1]),
    gene=as.character(temp$gene_name[1]),
    pN=sum(temp[,"pN"]),
    pS=sum(temp[,"pS"]),
    dN=sum(temp[,"dN"]),
    dS=sum(temp[,"dS"])))
}

# Calculate some stuff
MKT_MAF$pN_pS <- MKT_MAF$pN+MKT_MAF$pS
MKT_MAF$pN_dN <- MKT_MAF$pN+MKT_MAF$dN
MKT_MAF$dS_dN <- MKT_MAF$dS+MKT_MAF$dN
MKT_MAF$dS_pS <- MKT_MAF$dS+MKT_MAF$pS
MKT_MAF$gene_length <- gene_length$length[match(MKT_MAF$gene, gene_length$gene)]
MKT_MAF$SNVs_per_bp <- (MKT_MAF$pN+MKT_MAF$pS+MKT_MAF$dN+MKT_MAF$dS)/MKT_MAF$gene_length

MKT_MAF$fisher.test.P<-99

for(i in 1:nrow(MKT_MAF)){
  MKT_MAF$fisher.test.P[i]<-fisher.test(matrix(as.numeric(MKT_MAF[i,c(6,4,5,3)]), ncol=2))$p.value 
  if((MKT_MAF$pN[i] == 0 && MKT_MAF$dN[i] == 0) || (MKT_MAF$pS[i] == 0 && MKT_MAF$dS[i] == 0) || (MKT_MAF$pS[i] == 0 && MKT_MAF$pN[i] == 0) || (MKT_MAF$dS[i] == 0 && MKT_MAF$dN[i] == 0)) { MKT_MAF$fisher.test.P[i]<-NA} 
} 

is.na(MKT_MAF) <- sapply(MKT_MAF, is.infinite)

MKT_MAF$neut_index=(MKT_MAF$pN*MKT_MAF$dS)/(MKT_MAF$pS*MKT_MAF$dN)
is.na(MKT_MAF$neut_index) <- sapply(MKT_MAF$neut_index, is.infinite)

MKT_MAF$dos=((MKT_MAF$dN/(MKT_MAF$dN+MKT_MAF$dS))-(MKT_MAF$pN/(MKT_MAF$pN+MKT_MAF$pS)))

MKT_MAF <- MKT_MAF[!is.na(MKT_MAF$fisher.test.P), ]
MKT_MAF <- filter(MKT_MAF, !(gene %in% exclude$gene))

write.csv(MKT_MAF, args[7])

### Exons ===================================================================

MKT_exon<-data.frame(chrom=factor(), gene=factor(), exon=factor(), pN=numeric(), pS=numeric(), dN=numeric(), dS=numeric())

for(i in 1:length(levels(vcf_SNVs$exon))) {
  temp<-vcf_SNVs[which(vcf_SNVs$exon==levels(vcf_SNVs$exon)[i]),]
  MKT_exon<-rbind(MKT_exon,data.frame(
    chrom=as.character(temp$chrom[1]),
    gene=as.character(temp$gene_name[1]),
    exon=as.character(temp$exon[1]),
    pN=sum(temp[,"pN"]),
    pS=sum(temp[,"pS"]),
    dN=sum(temp[,"dN"]),
    dS=sum(temp[,"dS"])))
}

MKT_exon$pN_pS <- MKT_exon$pN+MKT_exon$pS
MKT_exon$pN_dN <- MKT_exon$pN+MKT_exon$dN
MKT_exon$dS_dN <- MKT_exon$dS+MKT_exon$dN
MKT_exon$dS_pS <- MKT_exon$dS+MKT_exon$pS
MKT_exon$exon_length <- gff_file$exon_length[match(MKT_exon$exon, gff_file$exon)]
MKT_exon$SNVs_per_bp <- (MKT_exon$pN+MKT_exon$pS+MKT_exon$dN+MKT_exon$dS)/MKT_exon$exon_length

MKT_exon$fisher.test.P<-99

for(i in 1:nrow(MKT_exon)){
  MKT_exon$fisher.test.P[i]<-fisher.test(matrix(as.numeric(MKT_exon[i,c(7,5,6,4)]), ncol=2))$p.value 
  if((MKT_exon$pN[i] == 0 && MKT_exon$dN[i] == 0) || (MKT_exon$pS[i] == 0 && MKT_exon$dS[i] == 0) || (MKT_exon$pS[i] == 0 && MKT_exon$pN[i] == 0) || (MKT_exon$dS[i] == 0 && MKT_exon$dN[i] == 0)) { MKT_exon$fisher.test.P[i]<-NA} 
}

is.na(MKT_exon) <- sapply(MKT_exon, is.infinite)

MKT_exon$neut_index=(MKT_exon$pN*MKT_exon$dS)/(MKT_exon$pS*MKT_exon$dN)
is.na(MKT_exon$neut_index) <- sapply(MKT_exon$neut_index, is.infinite)

MKT_exon$dos=((MKT_exon$dN/(MKT_exon$dN+MKT_exon$dS))-(MKT_exon$pN/(MKT_exon$pN+MKT_exon$pS)))

MKT_exon <- MKT_exon[!is.na(MKT_exon$fisher.test.P), ]
MKT_exon <- filter(MKT_exon, !(gene %in% exclude$gene))

write.csv(MKT_exon, args[8])