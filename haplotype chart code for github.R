
library(ape)
library(pegas)
library(adegenet)
library(Bios2cor)

install.packages("adegenet", type="binary")
library(adegenet)
install.packages("pegas", type="binary")


#######NEW ANALYSiS####

#rarified samples
data_new_rar<-read.dna("~/Desktop/all_helminths_season_rar.fasta", format="fasta")
data_new_rar

#determine number of haplotypes
ASVhaps_new_rar <- haplotype(data_new_rar)
ASVhaps_new_rar
summary(ASVhaps_new_rar)

#plot
ASVNet_new_rar <- haploNet(ASVhaps_new_rar)
plot(ASVNet_new_rar, size=attr(ASVNet_new_rar, "freq"), fast = FALSE, show.mutation=0, scale.ratio = 0.05, cex=0.8, labels=FALSE)

ind.hap.s<-with(
  stack(setNames(attr(ASVhaps_new_rar, "index"), rownames(ASVhaps_new_rar))), 
  table(hap=ind, pop=rownames(data_new_rar)[values])
)
ind.hap.s 

ASVNet_new_rar <- haploNet(ASVhaps_new_rar)
plot(ASVNet_new_rar, size=attr(ASVNet_new_rar, "freq"), fast = FALSE, show.mutation=0, scale.ratio = 0.007, cex=0.5, labels=TRUE, pie=ind.hap.s)



