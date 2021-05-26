# Load libraries
library(dplyr)
library(stringr)
library(tidyr)

# Create dataframe from annotated gene list
gff_path="pantro6_autosomes.gff"
gff_file<-readLines(gff_path)
gff_file<-data.frame(gff_file[grep(pattern="chr",gff_file):length(gff_file)])
gff_file <- data.frame(do.call('rbind', strsplit(as.character(gff_file[,1]),'\t',fixed=TRUE)))
colnames(gff_file) <- c("A","B","C","D","E","F","G","H","J")
gff_file <- gff_file[grep('^[chr]', gff_file$A),]
gff_file <- subset(gff_file, gff_file$C == 'CDS')

gff_file$Parent <- str_extract(gff_file$J, "Parent=(.*?);")
gff_file$Gene <- str_extract(gff_file$J, "gene=(.*?);")
gff_file$Parent <- gsub("Parent=rna-", "", gff_file$Parent)
gff_file$Parent <- gsub(";", "", gff_file$Parent)
gff_file$Gene <- gsub("gene=", "", gff_file$Gene)
gff_file$Gene <- gsub(";", "", gff_file$Gene)

gff_file$Gene <- as.factor(as.character(gff_file$Gene))
gff_file$Parent <- as.factor(as.character(gff_file$Parent))
gff_file$D<- as.numeric(as.character(gff_file$D))
gff_file$E<- as.numeric(as.character(gff_file$E))
gff_file$CDS_range <- gff_file$E-gff_file$D

# Calculate transcript lengths and subset the longest
gff_file$transcript_sum <- with(gff_file, ave(CDS_range, Parent, FUN=sum))
gff_file$longest_transcript <- with(gff_file, ave(transcript_sum, Gene, FUN=max))
filtered_gff <- subset(gff_file, transcript_sum == longest_transcript)

# Ensure no duplicates made it through
super_filtered_gff <- filtered_gff %>% distinct(filtered_gff$A, filtered_gff$D, filtered_gff$E, .keep_all=TRUE)

# Make new gff file with one (the longest) transcript per gene
new_gff <- data.frame(super_filtered_gff[c(1:9)])
write.table(copy,'new_gff.gff',sep="\t",row.names=FALSE,col.names=FALSE,quote=FALSE)

# Make BED file for autosomal CDS
# We need this subset the VCFs for CDS regions for any MK tests
cds_autosomes <- data.frame(super_filtered_gff[c(1, 4:5)])
cds_autosomes$D <- copy2$D-1
cds_autosomes$E <- copy2$E
write.table(cds_autosomes,'CDS_autosomes.bed',sep="\t",row.names=FALSE,col.names=FALSE,quote=FALSE)