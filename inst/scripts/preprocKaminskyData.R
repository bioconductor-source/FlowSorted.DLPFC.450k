
library(minfi)
pd1 = read.delim("/dcs01/lieber/ajaffe/450k/Zack_450k/FinalReport_zk_methylation450 data1_samples table_03232012.txt",
		header = TRUE, as.is=TRUE, skip = 8)
pd1$X = NULL # drop extra blank column
pd2 = read.delim("/dcs01/lieber/ajaffe/450k/Zack_450k/FinalReport_zk_methylation450 data2_samples table-04042012.txt",
		header = TRUE, as.is=TRUE, skip = 8)
		
pd = rbind(pd1,pd2)
pd = pd[,c(2,4,5)]

pheno = read.delim("/dcs01/lieber/ajaffe/450k/Zack_450k/CETS_Annotation_in.txt",
	header=TRUE, as.is=TRUE)
colnames(pd)[1] = "Microarray.ID"
pd = merge(pd, pheno, by="Microarray.ID", all=TRUE, sort= FALSE)

pd$BasePath = paste("/dcs01/lieber/ajaffe/450k/Zack_450k/IDAT/",
	pd$Sentrix.Barcode,"/",
	pd$Sentrix.Barcode,"_",pd$Sample.Section,sep="")

# only control samples
pd = pd[pd$celltype %in% c("N","G") &
		pd$diag == "Control",]
pd$SampleID = pd$ID2
pd$BrainID = pd$IDraw
pd$CellType = ifelse(pd$celltype=="N", "NeuN_pos","NeuN_neg")

pd$ID1 = pd$Microarray.ID = pd$ID2 = NULL
pd$celltype = pd$ID3 = pd$IDraw = NULL

pd = pd[,c(9:11, 1:8)]	
colnames(pd)[1:2] = c("Sample_Name", "SampleID")
pd$Sample_Name = gsub("-", "_", pd$Sample_Name)

FlowSorted.DLPFC.450k <- read.450k(pd$BasePath, verbose=TRUE)
pData(FlowSorted.DLPFC.450k) = pd
sampleNames(FlowSorted.DLPFC.450k) <- pd$Sample_Name 
save(FlowSorted.DLPFC.450k, file = "data/FlowSorted.DLPFC.450k.rda")
