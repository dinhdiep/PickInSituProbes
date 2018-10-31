#!/usr/bin/Rscript

library(dplyr)
library(tidyr)
library(stringr)

data_dir <- "Data"

rsem_GRCh38 <- read.table(gzfile(paste0(data_dir, "/rsem_GRCh38.p2.gtf.gz")), sep = "\t", header = F, stringsAsFactors = FALSE)
colnames(rsem_GRCh38) <- make.names(c("Ref", "Source", "Region", "Start", "End", "C5", "Strand", "C7", "GeneInfo" ))
dim(rsem_GRCh38)

ncbi_GRCh38_cds <- read.table(gzfile(paste0(data_dir, "/ncbi_refseq_GRCh38_cds.gtf.gz")), sep = "\t", header = F, stringsAsFactors = FALSE)
colnames(ncbi_GRCh38_cds) <- make.names(c("Ref", "Source", "Region", "Start", "End", "C5", "Strand", "C7", "GeneInfo"))

# get only the exons
rsem_GRCh38_exon <- rsem_GRCh38 %>% filter(Region == "exon")
dim(rsem_GRCh38_exon)

gene_symbol <- sapply(rsem_GRCh38_exon$GeneInfo, function(x) 
                      strsplit(strsplit(x, "; ")[[1]][2], " ")[[1]][2] )
transcript_id <- sapply(rsem_GRCh38_exon$GeneInfo, function(x) 
                      strsplit(strsplit(x, "; ")[[1]][3], " ")[[1]][2] )
transcript_id_cds <- sapply(ncbi_GRCh38_cds$GeneInfo, function(x) 
                      strsplit(strsplit(x, "; ")[[1]][2], " ")[[1]][2] )

rsem_GRCh38_exon$GeneSymbol <- gene_symbol
rsem_GRCh38_exon$TranscriptID <- transcript_id
modexpressers <- read.table("genes.txt")$V1

print("Total number of genes in list:")
length(modexpressers)

unfound <- modexpressers[! modexpressers %in% gene_symbol]
print("Missing gene from table:")
unfound

selected_exon <- rsem_GRCh38_exon %>% filter(GeneSymbol %in% modexpressers)

ncbi_GRCh38_cds$TranscriptID <- transcript_id_cds
selected_cds <- ncbi_GRCh38_cds %>% filter(transcript_id_cds %in% selected_exon$TranscriptID) 
selected_cds$GeneSymbol <- selected_exon$GeneSymbol[match(selected_cds$TranscriptID, selected_exon$TranscriptID)]
ref_conversion_table <- read.table(paste0(data_dir, "/ref_conversion_table"), sep=",", header=F)
selected_cds$Ref_2 <- ref_conversion_table$V2[match(selected_cds$Ref, ref_conversion_table$V1)]
selected_cds <- selected_cds %>% count(GeneSymbol, Ref_2, Start, End, Strand, sort = TRUE)
write.table(selected_cds, file = "cds.txt", sep = "\t", quote = FALSE, row = FALSE)

unfound <- modexpressers[! modexpressers %in% selected_cds$GeneSymbol]
print("Genes without known CDS:")
unfound

selected_exon <- rsem_GRCh38_exon %>% filter(GeneSymbol %in% unfound)
selected_exon$ExonSize <- sapply(1:nrow(selected_exon), function(x)
                                  selected_exon$End[x] - selected_exon$Start[x])

transcript_lengths <- selected_exon %>% group_by(TranscriptID, GeneSymbol) %>% summarise_at(vars(ExonSize), funs(sum(., na.rm = TRUE))) 
min_transcripts <- transcript_lengths %>% group_by(GeneSymbol) %>% filter(ExonSize == min(ExonSize))
min_transcripts_exon <- selected_exon %>% filter(TranscriptID %in% min_transcripts$TranscriptID) %>% count(GeneSymbol, Ref, Start, End, Strand, sort = TRUE)
write.table(min_transcripts_exon, file = "shortest_isoforms.txt", sep = "\t", quote = FALSE, row = FALSE)


selected_exon <- rsem_GRCh38_exon %>% filter(GeneSymbol %in% modexpressers)
exon_counts <- selected_exon %>% count(GeneSymbol, Ref, Start, End, Strand, sort = TRUE)
max_exon_counts <- exon_counts %>% group_by(GeneSymbol) %>% filter(n == max(n))
write.table(max_exon_counts, file = "constitutive.txt", sep = "\t", quote = FALSE, row = FALSE)
