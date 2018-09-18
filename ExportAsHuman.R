library('biomaRt')
library('limma')

load('countsDat.rda')

convertMouseGeneList <- function(x){
  
  human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")
  
  genesV2 = getLDS(attributes = c("mgi_symbol"), filters = "mgi_symbol", values = x , mart = mouse, attributesL = c("hgnc_symbol"), martL = human, uniqueRows=F)
  #humanx <- unique(genesV2[, 2])
  
  # Print the first 6 genes found to the screen
  #print(head(humanx))
  return(genesV2)
}

humanGenes <- convertMouseGeneList(row.names(countsDat))

countsDatHG <- countsDat[row.names(countsDat) %in% humanGenes[,1],]
row.names(countsDatHG) <- humanGenes[match(row.names(countsDatHG),humanGenes[,1]),2]

write.table(cbind(row.names(countsDatHG),countsDatHG),sep=",",file='countsDatHG.csv',row.names=F)

xcell <- read.table('xCell_countsDatHG_xCell_1034071818.txt',header=T,sep="\t",row.names=1)
xcell <- xcell[apply(xcell,1,max)>0.01,]

tumorSamp <- grep('WT',grep('LVR',colnames(xcell),value=T),value=T,invert=T)
tumorType <- gsub('_','',gsub('[0-9]','',tumorSamp))

mm <- model.matrix(~0+tumorType) 
lmf <- lmFit(xcell[,tumorSamp],mm)
ct <- makeContrasts(tumorTypeDQDLVR+tumorTypeLMTLVR-tumorTypeFGLVR-tumorTypeLMLVR,levels=mm)

topTable(eBayes(contrasts.fit(lmf,contrasts = ct)))
topTable(eBayes(contrasts.fit(lmf[c('ImmuneScore','StromaScore','MicroenvironmentScore'),],contrasts = ct)))


