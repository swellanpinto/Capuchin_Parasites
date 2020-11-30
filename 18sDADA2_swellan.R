
###########
## Load and install packages
###########

#source("http://bioconductor.org/biocLite.R")
#biocLite(suppressUpdates = FALSE)
#biocLite("ShortRead", suppressUpdates = FALSE)
#biocLite("devtools")
#library("devtools")
#devtools::install_github("benjjneb/dada2")

library("knitr")
library("BiocStyle")
library("dada2")
opts_chunk$set(cache = FALSE,fig.path="dadafigure/")

#ERROR# read_chunk(file.path("src", "bioinformatics.R"))

.cran_packages <- c("ggplot2", "gridExtra")
.bioc_packages <- c("dada2", "phyloseq", "DECIPHER", "phangorn")
#.inst <- .cran_packages %in% installed.packages()
#if(any(!.inst)) {
#  install.packages(.cran_packages[!.inst])
#}
#.inst <- .bioc_packages %in% installed.packages()
#if(any(!.inst)) {
#  source("http://bioconductor.org/biocLite.R")
#  biocLite(.bioc_packages[!.inst], ask = F)
#}
sapply(c(.cran_packages, .bioc_packages), require, character.only = TRUE)
set.seed(100)
setwd("~/18s")

#ps <- readRDS("ps_dada2.rds")
# metadata <- sample_data(ps)

################
##Download and sort files
################

miseq_path <- ("/Users/joseph.orkin/18s/trimmed")
filt_path <- ("/Users/joseph.orkin/18s/filtered")
fns <- sort(list.files(miseq_path, full.names = TRUE))
fnFs <- fns[grepl("_R1", fns)]
fnRs <- fns[grepl("_R2", fns)]
ii <- sample(length(fnFs), 3)

#ReadsPerFilePlot <- ggplot (ReadsPerFile, aes(x=Sample, y=Reads, color=Type)) + geom_point() + geom_vline(xintercept=62) + 
#  xlab("Sample") +ylab("Reads") + ggtitle("MelinLab MiSeq April 03, 2017")
#RPFsortPlot <- ggplot (RPFsort, aes(x=order, y=Reads, color=Type)) + geom_point()  + 
#  xlab("Sample") +ylab("Reads") +ggtitle("MelinLab MiSeq April 03, 2017")

#######
## Check error distribution
#######
for(i in ii) { print(plotQualityProfile(fnFs[i]) + labs(title="Fwd", subtitle=fnFs[i]))} 
for(i in ii) { print(plotQualityProfile(fnRs[i]) + labs(title="Rev", subtitle=fnFs[i]))}

##########
## Trim & Filter Reads
#########

if(!file_test("-d", filt_path)) dir.create(filt_path)
filtFs <- file.path(filt_path, basename(fnFs))
filtRs <- file.path(filt_path, basename(fnRs))

out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(226,226), 
                     #out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(227,227),
                     maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=TRUE, verbose=TRUE)
out
#########
## Learn Error Rates
########

errF <- learnErrors(filtFs, multithread=TRUE, pool=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE, pool=TRUE)

# plot error rates
plotErrors(errF, nominalQ=TRUE)

###########
## Dereplicate the filtered fastqs and give sample names
##########

derepFs <- derepFastq(filtFs, verbose=TRUE)
derepRs <- derepFastq(filtRs, verbose=TRUE)
sam.names <- sapply(strsplit(basename(filtFs), "_"), `[`, 1) # I changed this to the item after first _ to accomodate our read names
names(derepFs) <- sam.names
names(derepRs) <- sam.names

derepFs[[50]]
derepRs[[50]]

############
## Infer the composition of the sample (i.e. remove all sampling errors) with core dada2 algorithm
############

dadaFs <- dada(derepFs, err=errF, selfConsist=TRUE, multithread = TRUE, pool=TRUE)
dadaRs <- dada(derepRs, err=errR, selfConsist=TRUE, multithread = TRUE, pool=TRUE)

#########
## Merge the unique sequences
#########
mergers<- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose =TRUE, justConcatenate = TRUE)
#mergers<- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose =TRUE)
head(mergers[[1]])



#########
## Construct the sequence table and remove chimeras
#########

########### CHECK HERE #############

# make sequence table and remove any mock samples
seqtab.all <- makeSequenceTable(mergers)
seqtab.fwd <- makeSequenceTable(dadaFs)
seqtab.rev <- makeSequenceTable(dadaRs)

seqtab.collapsed <- collapseNoMismatch(seqtab.all, minOverlap = 20, orderBy = "abundance", vec = TRUE, verbose = TRUE)
seqtab.fwd.collapsed <- collapseNoMismatch(seqtab.fwd, minOverlap = 20, orderBy = "abundance", vec = TRUE, verbose = TRUE)
seqtab.rev.collapsed <- collapseNoMismatch(seqtab.rev, minOverlap = 20, orderBy = "abundance", vec = TRUE, verbose = TRUE)


dim(seqtab.all)
table(nchar(colnames(seqtab.all)))

dim(seqtab.fwd)
table(nchar(colnames(seqtab.fwd)))

dim(seqtab.rev)
table(nchar(colnames(seqtab.rev)))

#trim misprimed seqs
#seqtab.trim <- seqtab.all[,nchar(colnames(seqtab.all)) %in% seq(369,373)]
#seqtab.trim <- seqtab.all[,nchar(colnames(seqtab.all)) %in% seq(351,355)]
#seqtab.trim <- seqtab.all[,nchar(colnames(seqtab.all)) %in% seq(387,392)]
#table(nchar(colnames(seqtab.trim)))
#dim(seqtab.trim)

#Remove Chimeras

seqtab.nochim <- removeBimeraDenovo(seqtab.all, method="consensus", multithread=TRUE, verbose=TRUE)
seqtab.nochim.fwd <- removeBimeraDenovo(seqtab.fwd, method="consensus", multithread=TRUE, verbose=TRUE)
seqtab.nochim.rev <- removeBimeraDenovo(seqtab.rev, method="consensus", multithread=TRUE, verbose=TRUE)
#seqtab.nochim <- removeBimeraDenovo(seqtab.trim, method="consensus", multithread=TRUE, verbose=TRUE)

dim(seqtab.nochim)
dim(seqtab.nochim.fwd)
dim(seqtab.nochim.rev)

#write.csv(seqtab.nochim,"OTUtable_DADA2-PR2.csv")

#Track reads through pipeline
#getN <- function(x) sum(getUniques(x))
#track <- cbind(out, sapply(dadaFs, getN), sapply(mergers, getN), rowSums(seqtab.all), rowSums(seqtab.nochim))
#colnames(track) <- c("input", "filtered", "denoised", "merged", "tabled", "nonchim")
#rownames(track) <- sample.names
#head(track)

############
## Assign taxonomy
###########

#ref_fasta <- "~/18s/silva_nr_v123_train_set-2_ANDCustomParasite.fa"
#ref_fasta <- "~/18s/silva_nr_v123_train_set-2_PLUS_WormsFromFullSilva.fa"
ref_fasta <- "~/18s/ENVGEN/PR2_derep_3000bp_DADA2.fasta"

taxtab.cat <- assignTaxonomy(seqtab.nochim, refFasta = ref_fasta, multithread=TRUE)
taxtab.fwd <- assignTaxonomy(seqtab.nochim.fwd, refFasta = ref_fasta, multithread=TRUE)
taxtab.rev <- assignTaxonomy(seqtab.nochim.rev, refFasta = ref_fasta, multithread=TRUE)
#taxtab <- assignTaxonomy(seqtab.nochim, refFasta = ref_fasta, multithread=TRUE)

#ref_fasta <- "~/Microbiome/Databases/rdp_train_set_14-2.fa"
#ref_fasta <- "~/18s/silva_nr_v123_train_set-2.fa"
#ref_fasta <- "~/18s/silva.full_current.Nematoda.fasta"
#ref_fasta <- "~/18s/silva_nr_v123_train_set-2_ANDCustomParasite.fa
#ref_fasta <- "~/18s/WormsFromSilva.fa"
#ref_fasta <- "~/18s/silva_nr_v123_train_set-2_PLUS_WormsFromFullSilva.fa"

#write.csv(taxtab, "DADA2parasites2017-10-04-369-373-Silva_plus_WormsFromSilva.csv")
write.csv(seqtab.nochim,"DADA2-Concatenated-PR2-OTU_TAB.csv")

write.csv(taxtab.cat, "DADA2parasites2017-10-11_Silva_plus_WormsFromSilvaCONCATENATED.csv")
write.csv(taxtab.fwd, "DADA2parasites2017-10-11_Silva_plus_WormsFromSilvaFORWARD.csv")
write.csv(taxtab.rev, "DADA2parasites2017-10-11_Silva_plus_CustomWormsREVERSE.csv")


#taxtab.sp <- addSpecies(taxtab, "~/18s/silva_species_assignment_v123-2.fa", verbose=TRUE) #adding sp is much slower than the rest
#taxtab.trim <- assignTaxonomy(seqtab.trim, refFasta = ref_fasta)
#taxtab.trim.sp <- addSpecies(taxtab.trim, "rdp_species_assignment_14.fa.gz", verbose=TRUE) #adding sp is much slower than the rest

###########
## Build Tree
###########

# phangorn will construct a NJ tree then a GTR+G+I ML tree

seqs <- getSequences(seqtab.nochim)
names(seqs) <- seqs # This propagates to the tip labels of the tree
alignment <- AlignSeqs(DNAStringSet(seqs), anchor=NA)
phang.align <- phyDat(as(alignment, "matrix"), type="DNA")
dm <- dist.ml(phang.align)
treeNJ <- NJ(dm) # Note, tip order != sequence order
fit = pml(treeNJ, data=phang.align)

## negative edges length changed to 0!

fitGTR <- update(fit, k=4, inv=0.2)
fitGTR <- optim.pml(fitGTR, model="GTR", optInv=TRUE, optGamma=TRUE, rearrangement = "stochastic", control = pml.control(trace = 0)) # This can take a long time
detach("package:phangorn", unload=TRUE)


##############
## Convert to phyloseq object
##############

# Input metadata and fix row names

metadata = data.frame(read.csv("~/18s/SwellanMetadata.csv", header = TRUE))
#rownames(metadata) <- metadata[,1]
rownames(metadata) <-rownames(seqtab.nochim)
dim(seqtab.nochim)
dim(metadata)

ps <- phyloseq(tax_table(taxtab), sample_data(metadata), otu_table(seqtab.nochim, taxa_are_rows = FALSE), phy_tree(fitGTR$tree))
saveRDS(ps, "ps_18s_Swellan.rds")
#ps <- readRDS("~/18s/ps_18s_Swellan.rds")

ps0 <- ps
ps0 <- prune_samples(sample_sums(ps0) > 2000, ps0)
sample_names(ps0)
sort(sample_sums(ps0))


##############
## Filtering taxa
##############

# Take a look at what's there
)

ps1 <- subset_taxa(ps0, !is.na(Phylum) & !Phylum %in% c("", "uncharacterized"))
sort (sample_sums(ps1))
table(tax_table(ps1)[, "Phylum"], exclude = NULL)
table(tax_table(ps1)[, "Genus"], exclude = NULL)
write.csv(tax_table(ps1), "WormTaxTable.csv")
write.csv(otu_table(ps1), "WormOTUTable.csv")

tree1 = phy_tree(ps1)
ape::write.tree(tree1, "~/18s/worm.tree")


dataFiltered = ps1 %>% filter_taxa(function(x) sum(x > 0) > (0.05 * length(x)), TRUE)
table(tax_table(dataFiltered)[, "Phylum"], exclude = NULL)
table(tax_table(dataFiltered)[, "Species"], exclude = NULL)
write.csv(tax_table(dataFiltered), "WormTaxTableFiltered.csv")
write.csv(otu_table(dataFiltered), "WormOTUTableFiltered.csv")


# ps0 <- ps0 %>% 
#  subset_taxa(is.na(Family) | Family != "mitochondria") %>% 
#  subset_taxa(Phylum != "Cyanobacteria/Chloroplast")

# Compute prevalence of each feature, store as data.frame
prevdf = apply(X = otu_table(ps1),
               MARGIN = ifelse(taxa_are_rows(ps0), yes = 1, no = 2),
               FUN = function(x){sum(x > 0)})



# Add taxonomy and total read counts to this data.frame
prevdf = data.frame(Prevalence = prevdf,
                    TotalAbundance = taxa_sums(ps1),
                    tax_table(ps1))

# Compute the total and average prevalences of the features in each phylum, then remove phyla with 10 or fewer members
plyr::ddply(prevdf, "Phylum", function(df1){cbind(mean(df1$Prevalence),sum(df1$Prevalence))})
# filterPhyla = c("Fusobacteria", "Deinococcus-Thermus")

# Define phyla to filter
prevdf1 = subset(prevdf, Phylum %in% get_taxa_unique(ps1, "Phylum"))

# Define prevalence threshold as 5% of total samples
prevalenceThreshold = 0.10 * nsamples(ps0)
prevalenceThreshold



tmp <- installed.packages()
installedpkgs <- as.vector(tmp[is.na(tmp[,"Priority"]), 1])
save(installedpkgs, file="installed_old.rda")
