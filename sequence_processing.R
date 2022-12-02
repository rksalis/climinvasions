#18S metabarcoding sequence processing 
#R.K.Salis March 2022
BiocManager::install(version = '3.14')
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("dada2", version = "3.14")
library(dada2); packageVersion("dada2")
library("ape")
library(seqinr)
library(phyloseq)

path <- "" # directory containing the fastq files after unzipping.
list.files(path)

# Forward and reverse fastq filenames have format: SAMPLENAME_R1_001.fastq and SAMPLENAME_R2_001.fastq
fnFs <- sort(list.files(path, pattern="_R1_001.fastq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_R2_001.fastq", full.names = TRUE))
# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)

plotQualityProfile(fnFs[1:2])
plotQualityProfile(fnRs[1:2])

#Filter and trim
#Assign the filenames for the filtered fastq.gz files.
# Place filtered files in filtered/ subdirectory
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(290,260),
                     maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=TRUE) 
head(out)
#Learn the Error Rates
errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)
plotErrors(errF, nominalQ=TRUE)
#Sample Inference
dadaFs <- dada(filtFs, err=errF, multithread=TRUE)
dadaRs <- dada(filtRs, err=errR, multithread=TRUE)
dadaFs[[1]]
#Merge paired reads
mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE)
# Inspect the merger data.frame from the first sample
head(mergers[[1]])
#Construct sequence table
seqtab.18S <- makeSequenceTable(mergers)
dim(seqtab.18S)
write.csv(seqtab.18S, "seqtab.18S.csv")
# Inspect distribution of sequence lengths #amplicon size = 438bp
table(nchar(getSequences(seqtab.18S)))
#Remove chimeras
seqtab.nochim.18S <- removeBimeraDenovo(seqtab.18S, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim.18S)
sum(seqtab.nochim.18S)/sum(seqtab.18S)
saveRDS(seqtab.nochim.18S, "seqtab.nochim.18S.rds")
#Track reads through the pipeline
getN <- function(x) sum(getUniques(x))
track18S <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
colnames(track18S) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track18S) <- sample.names
head(track18S)
saveRDS(track18S, "track18S.rds")
write.csv(track18S, "track18S.csv")

#Assign taxonomy
st18S <- readRDS("seqtab.nochim.18S.rds")
#using PR2 version 4.13.0 (protist reference database)
taxa.18S.pr2.all <- assignTaxonomy(st18S, "pr2_version_4.13.0_18S_dada2.fasta.gz", multithread=TRUE, taxLevels = c("Kingdom","Supergroup","Division","Class","Order","Family","Genus","Species"))
saveRDS(taxa.18S.pr2.all, "taxa.18S.pr2.all.rds")
taxa.18S.pr2.all.print <- taxa.18S.pr2.all 
rownames(taxa.18S.pr2.all.print) <- NULL
head(taxa.18S.pr2.all.print)
#create fasta file
seqnum <- paste0("ASV", seq(ncol(st18S)))
uniqueSeqs <- as.list(colnames(st18S))
write.fasta(uniqueSeqs, seqnum, "18Sfulldataset.fasta")
#using pr2
#create phyloseq object : taxa.18S.pr2
st18S <- readRDS("seqtab_final.rds")
taxa.18S.pr2.all <- readRDS("taxa.18S.pr2.all.rds")
eDNA1_18S_pr2 = phyloseq(tax_table(taxa.18S.pr2.all), otu_table(st18S, taxa_are_rows = FALSE))
eDNA1_18S_pr2
rownames(otu_table(eDNA1_18S_pr2))
rownames(tax_table(eDNA1_18S_pr2))
## Upload metadata ##
map.eDNA1_18S_pr2 = read.csv("mapping_eDNA1_18S.csv")
head(map.eDNA1_18S_pr2)
###Rename samples for metadata
row.names(map.eDNA1_18S_pr2)<- map.eDNA1_18S_pr2$Sample #sample data row names must align with dada2 rowname outputs
head(map.eDNA1_18S_pr2)
eDNA1_18S_pr2 = phyloseq(tax_table(taxa.18S.pr2.all), otu_table(st18S, taxa_are_rows = FALSE), sample_data(map.eDNA1_18S_pr2))
eDNA1_18S_pr2
###Rename sequence variants
a.vec = as.vector(1:7282)  #total ASVs
a.nam = cbind("ASV", a.vec)
a.nam = as.data.frame(a.nam)
asv.names = paste0(a.nam$V1, a.nam$a.vec)
asv.names = as.data.frame(asv.names)
head(asv.names)
# apply ASV names to sequence table
taxa_names(eDNA1_18S_pr2) = asv.names$asv.names
taxa_names(eDNA1_18S_pr2)
eDNA1_18S_pr2
#prune samples
eDNA1_18S_pr2 <- prune_samples(sample_sums(eDNA1_18S_pr2)>=1, eDNA1_18S_pr2)
# write files for full dataset, no taxonomic filtering or sample filtering
write.csv(otu_table(eDNA1_18S_pr2), "asv_table_eDNA1_18S_pr2.csv")
write.csv(tax_table(eDNA1_18S_pr2), "tax_table_eDNA1_18S_pr2.csv")
write.csv(sample_data(eDNA1_18S_pr2), "sample_data_eDNA1_18S_pr2.csv")
saveRDS(eDNA1_18S_pr2, "eDNA1_18S_pr2.rds")

##### assign forwards only - to get those missed through merging
#Filter and trim
#Assign the filenames for the filtered fastq.gz files.
# Place filtered files in filtered/ subdirectory
filtFs <- file.path(path, "filteredfwds", paste0(sample.names, "_F_filt.fastq.gz"))
names(filtFs) <- sample.names
#Filtering paraments: maxN=0 (DADA2 requires sequences contain no Ns), truncQ = 2,  rm.phix = TRUE and maxEE=2 (this means a maximum number of “expected errors” is allowed in a read)
outfwd <- filterAndTrim(fnFs, filtFs, truncLen=290, maxN=0, maxEE=2, truncQ=2, rm.phix=TRUE,
                        compress=TRUE, multithread=TRUE) 
head(outfwd)
saveRDS(outfwd, "outfwd.rds")
#Learn the Error Rates
errF <- learnErrors(filtFs, multithread=TRUE)
plotErrors(errF, nominalQ=TRUE)
# Infer sequence variants
dds <- vector("list", length(sample.names))
names(dds) <- sample.names
for(sam in sample.names) {
  cat("Processing:", sam, "\n")
  derep <- derepFastq(filtFs[[sam]])
  dds[[sam]] <- dada(derep, err=errF, multithread=TRUE)
}
## Generate sequence table ##
seqtab.fwd <- makeSequenceTable(dds)
saveRDS(seqtab.fwd, "seqtab.fwd.rds")
dim(seqtab.fwd)
#Remove chimeras
seqtab.nochim.18S.fwd <- removeBimeraDenovo(seqtab.fwd, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim.18S.fwd)
sum(seqtab.nochim.18S.fwd)/sum(seqtab.fwd)
saveRDS(seqtab.nochim.18S.fwd, "seqtab.nochim.18S.fwd.rds")
#Track reads through the pipeline
getN <- function(x) sum(getUniques(x))
track18S.fwd <- cbind(outfwd, sapply(dds, getN), rowSums(seqtab.nochim.18S.fwd))
colnames(track18S.fwd) <- c("input", "filteredF", "denoisedF", "nonchim")
rownames(track18S.fwd) <- sample.names
head(track18S.fwd)
saveRDS(track18S.fwd, "track18S.fwd.rds")
write.csv(track18S.fwd, "track18S.fwd.csv")
#Assign taxonomy
seqtab.nochim.18S.fwd <- readRDS("seqtab.nochim.18S.fwd.rds")
#using PR2 version 4.13.0 (protist reference database)
taxa.18S.pr2.fwd <- assignTaxonomy(seqtab.nochim.18S.fwd, "pr2_version_4.13.0_18S_dada2.fasta.gz", multithread=TRUE, taxLevels = c("Kingdom","Supergroup","Division","Class","Order","Family","Genus","Species"))
saveRDS(taxa.18S.pr2.fwd, "taxa.18S.pr2.fwd.rds")
taxa.print.pr2.fwd <- taxa.18S.pr2.fwd 
rownames(taxa.print.pr2.fwd) <- NULL
head(taxa.print.pr2.fwd)
#create fasta file
seqnum <- paste0("ASV", seq(ncol(seqtab.nochim.18S.fwd)))
uniqueSeqs <- as.list(colnames(seqtab.nochim.18S.fwd))
write.fasta(uniqueSeqs, seqnum, "18Sfulldataset.fwd.fasta")
#create phyloseq object : taxa.18S.pr2.fwd
seqtab.nochim.18S.fwd <- readRDS("seqtab.nochim.18S.fwd.rds")
taxa.18S.pr2.fwd <- readRDS("taxa.18S.pr2.fwd.rds")
eDNA1_18S_pr2_fwd = phyloseq(tax_table(taxa.18S.pr2.fwd), otu_table(seqtab.nochim.18S.fwd, taxa_are_rows = FALSE))
eDNA1_18S_pr2_fwd
rownames(otu_table(eDNA1_18S_pr2_fwd))
rownames(tax_table(eDNA1_18S_pr2_fwd))
## Upload metadata ##
map.eDNA1_18S_pr2_fwd = read.csv("mapping_eDNA1_18S.csv")
head(map.eDNA1_18S_pr2_fwd)
###Rename samples for metadata
row.names(map.eDNA1_18S_pr2_fwd)<- map.eDNA1_18S_pr2_fwd$Sample #sample data row names must align with dada2 rowname outputs
head(map.eDNA1_18S_pr2_fwd)
eDNA1_18S_pr2_fwd = phyloseq(tax_table(taxa.18S.pr2.fwd), otu_table(seqtab.nochim.18S.fwd, taxa_are_rows = FALSE), sample_data(map.eDNA1_18S_pr2_fwd, phy_tree(random_tree)))
eDNA1_18S_pr2_fwd
###Rename sequence variants
a.vec = as.vector(1:5182)  #total ASVs
a.nam = cbind("ASVf", a.vec)
a.nam = as.data.frame(a.nam)
asv.names = paste0(a.nam$V1, a.nam$a.vec)
asv.names = as.data.frame(asv.names)
head(asv.names)
# apply ASV names to sequence table
taxa_names(eDNA1_18S_pr2_fwd) = asv.names$asv.names
taxa_names(eDNA1_18S_pr2_fwd)
eDNA1_18S_pr2_fwd
# write files for full dataset, no taxonomic filtering or sample filtering
write.csv(otu_table(eDNA1_18S_pr2_fwd), "asv_table_eDNA1_18S_pr2.fwd.csv")
write.csv(tax_table(eDNA1_18S_pr2_fwd), "tax_table_eDNA1_18S_pr2.fwd.csv")
write.csv(sample_data(eDNA1_18S_pr2_fwd), file = "sample_data_eDNA1_18S_pr2_fwd.csv")
saveRDS(eDNA1_18S_pr2_fwd, "eDNA1_18S_pr2_fwd.rds")

#identify and input ASVs unique to fwd reads (not retained in merged dataset - an additional 391 ASVs)
#filter rows in forward tax_table_eDNA1_18Sf that are also in tax_table_eDNA1_18S
tax_table_eDNA1_18Sf <- read.csv("tax_table_eDNA1_18S_pr2.fwd.csv")
tax_table_eDNA1_18S <- read.csv("tax_table_eDNA1_18S_pr2.csv")
tax.eDNA1_18Sfwduniques1 <- semi_join(tax_table_eDNA1_18Sf,tax_table_eDNA1_18S, by = c("Kingdom", "Supergroup", "Division", "Class", "Order", "Family", "Genus", "Species"))
tax.eDNA1_18Sfwduniques <- tax.eDNA1_18Sfwduniques1[, -1]
row.names(tax.eDNA1_18Sfwduniques) <- tax.eDNA1_18Sfwduniques1$X #sample data row names must align with dada2 rowname outputs
tax.eDNA1_18Sfwduniques = as.matrix(tax.eDNA1_18Sfwduniques)
head(tax.eDNA1_18Sfwduniques)
#asv table
asv.eDNA1_18Sf1 <- read.csv("asv_table_eDNA1_18S_pr2.fwd.csv")
head(asv.eDNA1_18Sf1)
asv.eDNA1_18Sf <- asv.eDNA1_18Sf1[, -1]
row.names(asv.eDNA1_18Sf)<- asv.eDNA1_18Sf1$X #sample data row names must align with dada2 rowname outputs
asv.eDNA1_18Sf = as.matrix(asv.eDNA1_18Sf)
head(asv.eDNA1_18Sf)
## Upload metadata ##
map.eDNA1.18Sf1 = read.csv("sample_data_eDNA1_18S_pr2_fwd.csv")
head(map.eDNA1.18Sf1)
###Rename samples for metadata
map.eDNA1.18Sf <- map.eDNA1.18Sf1[, -1]
row.names(map.eDNA1.18Sf)<- map.eDNA1.18Sf1$X #sample data row names must align with dada2 rowname outputs
map.eDNA1.18Sf = as.data.frame(map.eDNA1.18Sf)
head(map.eDNA1.18Sf)

#phyloseq object fwd uniques
eDNA1_18Sfwduniques = phyloseq(tax_table(tax.eDNA1_18Sfwduniques), otu_table(asv.eDNA1_18Sf,taxa_are_rows = FALSE), 
                               sample_data(map.eDNA1.18Sf))
eDNA1_18Sfwduniques
write.csv(otu_table(eDNA1_18Sfwduniques), "asv_table_eDNA1_18S_pr2.fwduniques.csv")
write.csv(tax_table(eDNA1_18Sfwduniques), "tax_table_eDNA1_18S_pr2.fwduniques.csv")
write.csv(sample_data(eDNA1_18Sfwduniques), file = "sample_data_eDNA1_18S_pr2.fwduniques.csv")
saveRDS(eDNA1_18Sfwduniques, "eDNA1_18S_pr2.fwduniques.rds")
eDNA1_18S_pr2<-readRDS("eDNA1_18S_pr2.rds")
eDNA1_18Sfwduniques<-readRDS("eDNA1_18S_pr2.fwduniques.rds")
#add in data from forward reads to merged dataset
eDNA1_18S_pr2_plusfwd <- merge_phyloseq(eDNA1_18S_pr2, eDNA1_18Sfwduniques)
eDNA1_18S_pr2_plusfwd
plot_bar(eDNA1_18S_pr2_plusfwd)
rarecurve(t(otu_table(eDNA1_18S_pr2_plusfwd)), step=50, cex=0.5)
C1.ord <- ordinate(eDNA1_18S_pr2_plusfwd, method ="PCoA", "jaccard")
plot_ordination(eDNA1_18S_pr2_plusfwd, C1.ord, type="samples", color="treatment")

###taxonomic filtering
rank_names(eDNA1_18S_pr2_plusfwd)
table(tax_table(eDNA1_18S_pr2_plusfwd)[, "Kingdom"], exclude = NULL)
eDNA1_18S_pr2_plusfwd
table(tax_table(eDNA1_18S_pr2_plusfwd)[, "Supergroup"], exclude = NULL)
table(tax_table(eDNA1_18S_pr2_plusfwd)[, "Division"], exclude = NULL)
#remove taxa for which Division is NA
eDNA1_18S_pr20 <- subset_taxa(eDNA1_18S_pr2_plusfwd, !is.na(Division) & !Division %in% c("", "uncharacterized"))
table(tax_table(eDNA1_18S_pr20)[, "Division"], exclude = NULL)
#remove Mammals (homo sapians) 
eDNA1_18S_pr20 <- subset_taxa(eDNA1_18S_pr20, (Class!="Craniata") | is.na(Order))
eDNA1_18S_pr20
#remove grass
eDNA1_18S_pr20 <- subset_taxa(eDNA1_18S_pr20, (Genus!="Brachypodium") | is.na(Species))
eDNA1_18S_pr20

#subset field controls (NEG)-
eDNA1_18S_pr2_NEG  = subset_samples(eDNA1_18S_pr20, treatment == "NEG")
eDNA1_18S_pr2_NEG 
eDNA1_18S_pr2_NEG  <- filter_taxa(eDNA1_18S_pr2_NEG, function (x) {sum(x > 0) > 0.5}, prune=TRUE)
eDNA1_18S_pr2_NEG
write.csv(otu_table(eDNA1_18S_pr2_NEG), "asv_table_eDNA1_18S_pr2_NEG.csv")
write.csv(tax_table(eDNA1_18S_pr2_NEG), "tax_table_eDNA1_18S_pr2_NEG.csv")
write.csv(sample_data(eDNA1_18S_pr2_NEG), "sample_data_eDNA1_18S_pr2_NEG.csv")
saveRDS(eDNA1_18S_pr2_NEG, "eDNA1_18S_pr2_NEG.rds")

#remove NEGs from dataset
eDNA1_18S_pr2_noNEG  = subset_samples(eDNA1_18S_pr20, treatment != "NEG")
eDNA1_18S_pr2_noNEG 
eDNA1_18S_pr2_mesocosms <- filter_taxa(eDNA1_18S_pr2_noNEG, function (x) {sum(x > 0) > 0.5}, prune=TRUE)
eDNA1_18S_pr2_mesocosms
write.csv(otu_table(eDNA1_18S_pr2_mesocosms), "asv_table_eDNA1_18S_pr2_mesocosms.csv")
write.csv(tax_table(eDNA1_18S_pr2_mesocosms), "tax_table_eDNA1_18S_pr2_mesocosms.csv")
write.csv(sample_data(eDNA1_18S_pr2_mesocosms), "sample_data_eDNA1_18S_pr2_mesocosms.csv")
saveRDS(eDNA1_18S_pr2_mesocosms, "eDNA1_18S_pr2_mesocosms.rds")

eDNA1_18S_pr2_mesocosms <- readRDS("eDNA1_18S_pr2_mesocosms.rds")

#Summarise asvs to species level
tax.clean_eDNA1_18S_pr2_mesocosms <- data.frame(tax_table(eDNA1_18S_pr2_mesocosms))
for (i in 1:8){ tax.clean_eDNA1_18S_pr2_mesocosms[,i] <- as.character(tax.clean_eDNA1_18S_pr2_mesocosms[,i])}
tax.clean_eDNA1_18S_pr2_mesocosms[is.na(tax.clean_eDNA1_18S_pr2_mesocosms)] <- ""
for (i in 1:nrow(tax.clean_eDNA1_18S_pr2_mesocosms)){
  if (tax.clean_eDNA1_18S_pr2_mesocosms[i,2] == ""){
    kingdom <- paste("Kingdom_", tax.clean_eDNA1_18S_pr2_mesocosms[i,1], sep = "")
    tax.clean_eDNA1_18S_pr2_mesocosms[i, 2:8] <- kingdom
  } else if (tax.clean_eDNA1_18S_pr2_mesocosms[i,3] == ""){
    Supergroup <- paste("Supergroup_", tax.clean_eDNA1_18S_pr2_mesocosms[i,2], sep = "")
    tax.clean_eDNA1_18S_pr2_mesocosms[i, 3:8] <- Supergroup
  } else if (tax.clean_eDNA1_18S_pr2_mesocosms[i,4] == ""){
    Division <- paste("Division_", tax.clean_eDNA1_18S_pr2_mesocosms[i,3], sep = "")
    tax.clean_eDNA1_18S_pr2_mesocosms[i, 4:8] <- Division
  } else if (tax.clean_eDNA1_18S_pr2_mesocosms[i,5] == ""){
    Class <- paste("Class_", tax.clean_eDNA1_18S_pr2_mesocosms[i,4], sep = "")
    tax.clean_eDNA1_18S_pr2_mesocosms[i, 5:8] <- Class
  } else if (tax.clean_eDNA1_18S_pr2_mesocosms[i,6] == ""){
    Order <- paste("Order_", tax.clean_eDNA1_18S_pr2_mesocosms[i,5], sep = "")
    tax.clean_eDNA1_18S_pr2_mesocosms[i, 6:8] <- Order
  } else if (tax.clean_eDNA1_18S_pr2_mesocosms[i,7] == ""){
    Family <- paste("Family_", tax.clean_eDNA1_18S_pr2_mesocosms[i,6], sep = "")
    tax.clean_eDNA1_18S_pr2_mesocosms[i, 7:8] <- Family
  } else if (tax.clean_eDNA1_18S_pr2_mesocosms[i,8] == ""){
    tax.clean_eDNA1_18S_pr2_mesocosms$Species[i] <- paste("Genus",tax.clean_eDNA1_18S_pr2_mesocosms$Genus[i], sep = "_")
  }
}
tax_table(eDNA1_18S_pr2_mesocosms) <- as.matrix(tax.clean_eDNA1_18S_pr2_mesocosms)
tax_table(eDNA1_18S_pr2_mesocosms)
eDNA1_18S_pr2_mesocosms_sp <- tax_glom(eDNA1_18S_pr2_mesocosms, taxrank=rank_names(eDNA1_18S_pr2_mesocosms)[8])
eDNA1_18S_pr2_mesocosms_sp
write.csv(otu_table(eDNA1_18S_pr2_mesocosms_sp), file = "asv_table_eDNA1_18S_pr2_mesocosms_sp.csv")
write.csv(tax_table(eDNA1_18S_pr2_mesocosms_sp), file = "tax_table_eDNA1_18S_pr2_mesocosms_sp.csv")
write.csv(sample_data(eDNA1_18S_pr2_mesocosms_sp), file = "sample_data_eDNA1_18S_pr2_mesocosms_sp.csv")
saveRDS(eDNA1_18S_pr2_mesocosms_sp, "eDNA1_18S_pr2_mesocosms_sp.rds")

# calculate normalised read counts with DESeq2 - varianceStabilizingTransformation - fulldataset
library(DESeq2)
eDNA1_18S_pr2_mesocosms <- readRDS("eDNA1_18S_pr2_mesocosms.rds")
dds = phyloseq_to_deseq2(eDNA1_18S_pr2_mesocosms, ~ 1)
# calculate geometric means prior to estimate size factors
gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}
geoMeans = apply(counts(dds), 1, gm_mean)
dds = estimateSizeFactors(dds, geoMeans = geoMeans)
sizeFactors(dds)
vsd <- varianceStabilizingTransformation(dds, blind=TRUE)
vs_counts <- assay(vsd)
vs_counts[vs_counts<0]=0
taxa.eDNA1_18S_pr2_mesocosms1 <- read.csv("tax_table_eDNA1_18S_pr2_mesocosms.csv")
taxa.eDNA1_18S_pr2_mesocosms <- taxa.eDNA1_18S_pr2_mesocosms1[, -1]
row.names(taxa.eDNA1_18S_pr2_mesocosms) <- taxa.eDNA1_18S_pr2_mesocosms1$X #sample data row names must align with dada2 rowname outputs
taxa.eDNA1_18S_pr2_mesocosms = as.matrix(taxa.eDNA1_18S_pr2_mesocosms)
head(taxa.eDNA1_18S_pr2_mesocosms)
eDNA1_18S_pr2_mesocosms_vs = phyloseq(tax_table(taxa.eDNA1_18S_pr2_mesocosms), otu_table(vs_counts, taxa_are_rows = TRUE))
eDNA1_18S_pr2_mesocosms_vs
map.eDNA1_18S_pr2_mesocosms_vs = read.csv("sample_data_eDNA1_18S_pr2_mesocosms.csv")
head(map.eDNA1_18S_pr2_mesocosms_vs)
row.names(map.eDNA1_18S_pr2_mesocosms_vs)<- map.eDNA1_18S_pr2_mesocosms_vs$Sample #sample data row names must align with dada2 rowname outputs
head(map.eDNA1_18S_pr2_mesocosms_vs)
eDNA1_18S_pr2_mesocosms_vs = phyloseq(tax_table(taxa.eDNA1_18S_pr2_mesocosms), otu_table(vs_counts, taxa_are_rows = TRUE), sample_data(map.eDNA1_18S_pr2_mesocosms_vs))
eDNA1_18S_pr2_mesocosms_vs
table(tax_table(eDNA1_18S_pr2_mesocosms_vs)[, "Species"], exclude = NULL)
eDNA1_18S_pr2_mesocosms_vs <- filter_taxa(eDNA1_18S_pr2_mesocosms_vs, function (x) {sum(x > 0) > 0}, prune=TRUE)
eDNA1_18S_pr2_mesocosms_vs
write.csv(otu_table(eDNA1_18S_pr2_mesocosms_vs), file = "asv_table_eDNA1_18S_pr2_mesocosms_vs.csv")
write.csv(tax_table(eDNA1_18S_pr2_mesocosms_vs), file = "tax_table_eDNA1_18S_pr2_mesocosms_vs.csv")
write.csv(sample_data(eDNA1_18S_pr2_mesocosms_vs), file = "sample_data_eDNA1_18S_pr2_mesocosms_vs.csv")
saveRDS(eDNA1_18S_pr2_mesocosms_vs, "eDNA1_18S_pr2_mesocosms_vs.rds")

eDNA1_18S_pr2_mesocosms_vs <- readRDS("eDNA1_18S_pr2_mesocosms_vs.rds")

#Summarise asvs to species level
tax.clean_eDNA1_18S_pr2_mesocosms_vs <- data.frame(tax_table(eDNA1_18S_pr2_mesocosms_vs))
for (i in 1:8){ tax.clean_eDNA1_18S_pr2_mesocosms_vs[,i] <- as.character(tax.clean_eDNA1_18S_pr2_mesocosms_vs[,i])}
tax.clean_eDNA1_18S_pr2_mesocosms_vs[is.na(tax.clean_eDNA1_18S_pr2_mesocosms_vs)] <- ""
for (i in 1:nrow(tax.clean_eDNA1_18S_pr2_mesocosms_vs)){
  if (tax.clean_eDNA1_18S_pr2_mesocosms_vs[i,2] == ""){
    kingdom <- paste("Kingdom_", tax.clean_eDNA1_18S_pr2_mesocosms_vs[i,1], sep = "")
    tax.clean_eDNA1_18S_pr2_mesocosms_vs[i, 2:8] <- kingdom
  } else if (tax.clean_eDNA1_18S_pr2_mesocosms_vs[i,3] == ""){
    Supergroup <- paste("Supergroup_", tax.clean_eDNA1_18S_pr2_mesocosms_vs[i,2], sep = "")
    tax.clean_eDNA1_18S_pr2_mesocosms_vs[i, 3:8] <- Supergroup
  } else if (tax.clean_eDNA1_18S_pr2_mesocosms_vs[i,4] == ""){
    Division <- paste("Division_", tax.clean_eDNA1_18S_pr2_mesocosms_vs[i,3], sep = "")
    tax.clean_eDNA1_18S_pr2_mesocosms_vs[i, 4:8] <- Division
  } else if (tax.clean_eDNA1_18S_pr2_mesocosms_vs[i,5] == ""){
    Class <- paste("Class_", tax.clean_eDNA1_18S_pr2_mesocosms_vs[i,4], sep = "")
    tax.clean_eDNA1_18S_pr2_mesocosms_vs[i, 5:8] <- Class
  } else if (tax.clean_eDNA1_18S_pr2_mesocosms_vs[i,6] == ""){
    Order <- paste("Order_", tax.clean_eDNA1_18S_pr2_mesocosms_vs[i,5], sep = "")
    tax.clean_eDNA1_18S_pr2_mesocosms_vs[i, 6:8] <- Order
  } else if (tax.clean_eDNA1_18S_pr2_mesocosms_vs[i,7] == ""){
    Family <- paste("Family_", tax.clean_eDNA1_18S_pr2_mesocosms_vs[i,6], sep = "")
    tax.clean_eDNA1_18S_pr2_mesocosms_vs[i, 7:8] <- Family
  } else if (tax.clean_eDNA1_18S_pr2_mesocosms_vs[i,8] == ""){
    tax.clean_eDNA1_18S_pr2_mesocosms_vs$Species[i] <- paste("Genus",tax.clean_eDNA1_18S_pr2_mesocosms_vs$Genus[i], sep = "_")
  }
}
tax_table(eDNA1_18S_pr2_mesocosms_vs) <- as.matrix(tax.clean_eDNA1_18S_pr2_mesocosms_vs)
tax_table(eDNA1_18S_pr2_mesocosms_vs)
eDNA1_18S_pr2_mesocosms_vs_sp <- tax_glom(eDNA1_18S_pr2_mesocosms_vs, taxrank=rank_names(eDNA1_18S_pr2_mesocosms_vs)[8])
eDNA1_18S_pr2_mesocosms_vs_sp
write.csv(otu_table(eDNA1_18S_pr2_mesocosms_vs_sp), file = "asv_table_eDNA1_18S_pr2_mesocosms_vs_sp.csv")
write.csv(tax_table(eDNA1_18S_pr2_mesocosms_vs_sp), file = "tax_table_eDNA1_18S_pr2_mesocosms_vs_sp.csv")
write.csv(sample_data(eDNA1_18S_pr2_mesocosms_vs_sp), file = "sample_data_eDNA1_18S_pr2_mesocosms_vs_sp.csv")
saveRDS(eDNA1_18S_pr2_mesocosms_vs_sp, "eDNA1_18S_pr2_mesocosms_vs_sp.rds")


#species richness
library(ggplot2)
library(vegan)
library(car)
library(lsmeans)
library(effects)
library(cowplot)
library(Rmisc)
library(lme4)
require(reshape2)

#richness 18S species
sp_18S1 <- read.table("asv_table_eDNA1_18S_pr2_mesocosms_vs_sp.csv", sep = ",",  header = TRUE,  na.strings = c("NR", "NA")) 
sp_18S <- sp_18S1[, -1]
rownames(sp_18S) <- sp_18S1$X
sp_18St <- t(sp_18S)
sp_18S_sampledata <- read.table("sample_data_eDNA1_18S_pr2_mesocosms_vs_sp.csv", sep = ",",  header = TRUE,  na.strings = c("NR", "NA")) 
richness <- rowSums(sp_18St > 0) 
S <-specnumber(sp_18St) 
H <- diversity(sp_18St) #shannon diversity
simpson<-diversity(sp_18St,index="simpson") #simpsons 1 - D
inversimp<-diversity(sp_18St,index="invsimpson") #inverse simpsons 1/D
J<-H/log(S) #Pielou's evenness
evenness<-inversimp/S #simpson's evenness (1/D)/S
diversity <- data.frame(richness, S, H, simpson, inversimp, J, evenness)
diversity
sp_18S_richness <- data.frame(sp_18S_sampledata,diversity)
sp_18S_richness
write.csv(sp_18S_richness, file = "Richness_18S_vs_sp.csv")
