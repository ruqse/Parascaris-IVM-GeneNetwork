
# Set working directory
setwd(glue::glue("{dirname(rstudioapi::getActiveDocumentContext()$path)}/../"))


#Loading all needed packages
suppressPackageStartupMessages({
  library(DESeq2)
  library(SummarizedExperiment)
  library(tximport)
  library(biomaRt)
  library(GenomicFeatures)
  library(AnnotationDbi)
  library(pheatmap)
  library(RColorBrewer)
  library(VennDiagram)
  library(ggplot2)
  library(hyperSpec)
  library(parallel)
  library(plotly)
  library(pvclust)
  library(vsn)
  library(here)
  library(gage)
  library(pathview)
  library(gageData)
  library(gplots)
  library(gridExtra)
  library(GEOquery)
  library(tidyverse)
  
  #library(tidyverse)
  })




load("data/metadata.rda")

load("data/txiObject.rda")

# Rename samples with intrareps IDs in the meta dataframe
rownames(meta) <- meta$intrareps
# Check if passed
all(rownames(meta) == colnames(txi$counts))




# Model build 

interactive.model = as.formula(~ expbatch + seqbatch + type + concentration + type:concentration)

dds <- DESeqDataSetFromTximport(txi, 
                                colData = meta,
                                design = interactive.model)


# save(dds, file=here("data/dds.rda"))


# Pre-filtering the dataset
# Reduce the size of the object, and to increase the speed of our functions
# Filter genes with at least 3 samples with a count of 10 or higher

nrow(dds)
keep <- rowSums(counts(dds) >= 10) >= 3

dds <- dds[keep,]
nrow(dds)

#=================================================================================================================
# Differential Gnene expresion
#===============================================================================================================                    

# Relevel
colData(dds)$type <- relevel(colData(dds)$type, ref = "anterior")

ddseq <- DESeq(dds)


all_res= results(ddseq, independentFiltering = T)

all_resdf = as.data.frame(all_res)
all_resDFnoNA <- na.omit(all_resdf)
# save(all_resDFnoNA, file = here("data/all_resDFnoNA.rda"))



####################################
#          Intestine              #
###################################

# For tissue type intestine, what is the difference between concentration IVM9 and the control :
# IVM9 - Control + type intestine interaction in ivm 9


# using Median of the normalized counts

intestine.IVM9vsCD <- results(ddseq,
                                             contrast = list(c("concentration_IVM9_vs_CD", "typeintestine.concentrationIVM9")),
                                             filter= rowMedians(counts(ddseq,normalized=TRUE)), alpha=0.05)


sum(intestine.IVM9vsCD$padj < 0.05 &
      (intestine.IVM9vsCD$log2FoldChange <= -0.5 | intestine.IVM9vsCD$log2FoldChange >= 0.5), na.rm = TRUE)

intestine.IVM9vsCD_subset = subset(intestine.IVM9vsCD, intestine.IVM9vsCD$padj < 0.05)
intestine.IVM9vsCD_subset = subset(intestine.IVM9vsCD_subset,intestine.IVM9vsCD_subset$log2FoldChange <= -0.5 |intestine.IVM9vsCD_subset$log2FoldChange >= 0.5)


intestine.IVM9vsCDDF = as.data.frame(intestine.IVM9vsCD_subset)

# write.csv(intestine.IVM9vsCDDF,  file="data/DEGs/intestine.IVM9vsCtrl.csv")

# For tissue type intestine, what is the difference between concentration IVM11 and the control :
# IVM11 - Control + type intestine interaction in ivm 11


# using Median of the normalized counts

intestine.IVM11vsCD <- results(ddseq,
                                              contrast = list(c("concentration_IVM11_vs_CD", "typeintestine.concentrationIVM11")),
                                              filter= rowMedians(counts(ddseq,normalized=TRUE)))




intestine.IVM11vsCD_subset = subset(intestine.IVM11vsCD, intestine.IVM11vsCD$padj < 0.05)
intestine.IVM11vsCD_subset = subset(intestine.IVM11vsCD_subset,intestine.IVM11vsCD_subset$log2FoldChange <= -0.5 |intestine.IVM11vsCD_subset$log2FoldChange >= 0.5)

intestine.IVM11vsCDDF = as.data.frame(intestine.IVM11vsCD_subset)

write.csv(intestine.IVM11vsCDDF, file="data/DEGs/intestine.IVM11vsCtrl.csv")
# 
# 
# 



####################################
#          Anterior end            #
###################################


# Relevel
colData(dds)$type <- relevel(colData(dds)$type, ref = "intestine")



ddseq <- DESeq(dds)


# For tissue type anterior, what is the difference between concentration IVM9 and the control :
# IVM9 - Control + type anterior interaction in ivm 9


#using Median of the normalized counts
# 
anterior.IVM9vsCD <- results(ddseq,
                                        contrast = list(c("concentration_IVM9_vs_CD", "typeanterior.concentrationIVM9")),
                                        filter= rowMedians(counts(ddseq,normalized=TRUE)))


# metadata(anterior.IVM9vsCD)$filterThreshold

# mcols(ddseq ,use.names=TRUE)[1:5,1:10]

sum(anterior.IVM9vsCD$padj < 0.05 &
      (anterior.IVM9vsCD$log2FoldChange <= -0.5 | anterior.IVM9vsCD$log2FoldChange >= 0.5), na.rm = TRUE)

anterior.IVM9vsCD_subset = subset(anterior.IVM9vsCD, anterior.IVM9vsCD$padj < 0.05)
anterior.IVM9vsCD_subset = subset(anterior.IVM9vsCD_subset,anterior.IVM9vsCD_subset$log2FoldChange <= -0.5 |anterior.IVM9vsCD_subset$log2FoldChange >= 0.5)


anterior.IVM9vsCDDF = as.data.frame(anterior.IVM9vsCD_subset)

write.csv(anterior.IVM9vsCDDF, file = "data/DEGs/anterior.IVM9vsCtrl.csv")



# For tissue type anterior, what is the difference between concentration IVM11 and the control :
# IVM11 - Control + type intestine interaction in ivm 11




# using Median of the normalized counts

anterior.IVM11vsCD <- results(ddseq,
                                         contrast = list(c("concentration_IVM11_vs_CD", "typeanterior.concentrationIVM11")),
                                         filter= rowMedians(counts(ddseq,normalized=TRUE)))


sum(anterior.IVM11vsCD$padj < 0.05 &
      (anterior.IVM11vsCD$log2FoldChange <= -0.5 | anterior.IVM11vsCD$log2FoldChange >= 0.5), na.rm = TRUE)

anterior.IVM11vsCD_subset = subset(anterior.IVM11vsCD, anterior.IVM11vsCD$padj < 0.05)
anterior.IVM11vsCD_subset = subset(anterior.IVM11vsCD_subset,anterior.IVM11vsCD_subset$log2FoldChange <= -0.5 |anterior.IVM11vsCD_subset$log2FoldChange >= 0.5)

anterior.IVM11vsCDDF = as.data.frame(anterior.IVM11vsCD_subset)

write.csv(anterior.IVM11vsCDDF, file = "data/DEGs/anterior.IVM11vsCtrl.csv")



# create all list of all results
Results_DEGs_list <- list(intestine.IVM11vsCDDF, intestine.IVM9vsCDDF, anterior.IVM11vsCDDF, anterior.IVM9vsCDDF)

# label the result Df with respective names
names(Results_DEGs_list ) <- c("intestine.IVM11vsCtrl", "intestine.IVM9vsCtrl", "anterior.IVM11vsCtrl", "anterior.IVM9vsCtrl")
save(Results_DEGs_list , file = "data/Results_DEGs_list.rda")





sessionInfo()