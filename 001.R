library(ChIPseeker)
setwd("e:/K9homer/")
mergeLK9KO4LWT12<-readPeakFile("mergeLK9KO4_LWT12.bed")
LK9KO4_LWT12<-readPeakFile("LK9KO4vsLWT12_USE.bed")
LW12_LK9KO4<-readPeakFile("LWT12vsLK9KO4_USE.bed")

mergeCK9KO2LWT12<-readPeakFile("mergeCK9KO2_LWT12.bed")
CK9KO2_LWT12<-readPeakFile("CK9KO2vsLWT12_USE.bed")
LW12_CK9KO2<-readPeakFile("LWT12vsCK9KO2_USE.bed")

mergeLK9KO4LWT15<-readPeakFile("mergeLK9KO4_LWT15.bed")
LK9KO4_LWT15<-readPeakFile("LK9KO4vsLWT15_USE.bed")
LW15_LK9KO4<-readPeakFile("LWT15vsLK9KO4_USE.bed")

mergeCK9KO2LWT15<-readPeakFile("mergeCK9KO2_LWT15.bed")
CK9KO2_LWT15<-readPeakFile("CK9KO2vsLWT15_USE.bed")
LW15_CK9KO2<-readPeakFile("LWT15vsCK9KO2_USE.bed")


library(org.Mm.eg.db)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(clusterProfiler)

peakAnno_LK9KO4_LWT12 <- annotatePeak(
  LK9KO4_LWT12,
  tssRegion = c(-3000, 3000),
  TxDb =TxDb.Mmusculus.UCSC.mm10.knownGene,
  annoDb = "org.Mm.eg.db")
write.table(peakAnno_LK9KO4_LWT12@anno@elementMetadata@listData,file = "peakAnno_LK9KO4_LWT12.txt",sep = "\t")
peakAnno_LK9KO4_LWT15 <- annotatePeak(
  LK9KO4_LWT15,
  tssRegion = c(-3000, 3000),
  TxDb =TxDb.Mmusculus.UCSC.mm10.knownGene,
  annoDb = "org.Mm.eg.db")
write.table(peakAnno_LK9KO4_LWT15@anno@elementMetadata@listData,file = "peakAnno_LK9KO4_LWT15.txt",sep = "\t")
peakAnno_CK9KO2_LWT12 <- annotatePeak(
  CK9KO2_LWT12,
  tssRegion = c(-3000, 3000),
  TxDb =TxDb.Mmusculus.UCSC.mm10.knownGene,
  annoDb = "org.Mm.eg.db")
write.table(peakAnno_CK9KO2_LWT12@anno@elementMetadata@listData,file = "peakAnno_CK9KO2_LWT12.txt",sep = "\t")
peakAnno_CK9KO2_LWT15 <- annotatePeak(
  CK9KO2_LWT15,
  tssRegion = c(-3000, 3000),
  TxDb =TxDb.Mmusculus.UCSC.mm10.knownGene,
  annoDb = "org.Mm.eg.db")
write.table(peakAnno_CK9KO2_LWT15@anno@elementMetadata@listData,file = "peakAnno_CK9KO2_LWT15.txt",sep = "\t")
