library('ComplexHeatmap')

read.table('xCell_countsDatHG_xCell_1034071818.txt',
           header=T, sep="\t", row.names=1, stringsAsFactors = F) -> xcell

xcell <- xcell[apply(xcell, 1, sd) > 1e-2,]

Heatmap(xcell, clustering_distance_rows = 'pearson', clustering_distance_columns = 'pearson')

Lung.Normal <- grep('WT[1-9]LNG',colnames(xcell),value=T)
Lung.Mut <- grep('WT',grep('LNG',colnames(xcell),value=T),value=T,invert=T)

Lung.Met.Cells.P <- p.adjust(apply(xcell,1,function(x){wilcox.test(x[Lung.Mut],x[Lung.Normal])$p.value}),method='BH')

Liver.Normal <- grep('WT[1-9]LVR',colnames(xcell),value=T)
Liver.Mut <- grep('WT',grep('LVR',colnames(xcell),value=T),value=T,invert=T)

Liver.Cells.P <- p.adjust(apply(xcell,1,function(x){wilcox.test(x[Liver.Mut],x[Liver.Normal])$p.value}),method='BH')

DQDLMT <- c(grep('DQD[1-9]LVR',colnames(xcell),value=T),grep('LMT[1-9]LVR',colnames(xcell),value=T))
F191GLM <- c(grep('F191G[1-9]LVR',colnames(xcell),value=T),grep('LM[1-9]LVR',colnames(xcell),value=T))
Liver.MetType.Cells.P <- p.adjust(apply(xcell,1,function(x){wilcox.test(x[F191GLM],x[DQDLMT])$p.value}),method='BH')
