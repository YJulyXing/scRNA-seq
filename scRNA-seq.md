# scRNA-seq

## Make index
### Get and filter files

```
wget ftp://ftp.flybase.net/genomes/Drosophila_melanogaster/dmel_r6.35_FB2020_04/fasta/dmel-all-chromosome-r6.35.fasta.gz
wget ftp://ftp.flybase.net/genomes/Drosophila_melanogaster/dmel_r6.35_FB2020_04/gtf/dmel-all-r6.35.gtf.gz
```

Made gene.list which contains gene IDs of all genes having CDS. The mRNAs with unknown strand (".") were deleted removed.

```
rm -f dmel-coding-r6.35.gtf
cat gene.list | while read line
do
grep "\"${line}\"" dmel-all-r6.35.gtf >> dmel-coding-r6.35.gtf
done
```

### Add GFP

```
cat GFP.fa | grep -v "^>" | tr -d "\n" | wc -c

echo -e 'GFP\tunknown\texon\t1\t720\t.\t+\t.\tgene_id "GFP"; transcript_id "GFP"; gene_name "GFP"; gene_biotype "protein_coding";' > GFP.gtf

cp dmel-all-chromosome-r6.35.fasta dmel-all-GFP-r6.35.fasta
cat GFP.fa >> dmel-all-GFP-r6.35.fasta

cp dmel-coding-r6.35.gtf dmel-coding-GFP-r6.35.gtf
cat GFP.gtf >> dmel-coding-GFP-r6.35.gtf
```

### Build index

```
cellranger mkref \
--genome=dmel-GFP-r6.35 \
--fasta=dmel-all-GFP-r6.35.fasta \
--genes=dmel-coding-GFP-r6.35.gtf
```

## CellRanger

Run the codes for each sample:

```	
cellranger count --id=Sample_P \
--fastqs=Sample_P \
--sample=Sample_P \
--transcriptome=dmel-GFP-r6.35
```

## Seurat

```
library(Seurat)

ids=c("Sample_P","Sample_KC","Sample_KC_12","Sample_12")

for (i in 1:length(ids)) {

	name=paste0(ids[i],"/outs/filtered_feature_bc_matrix")
	
	seurat_data <- Read10X(data.dir = name)
		
	sce <- CreateSeuratObject(counts = seurat_data,
	                                    project = ids[i])

	#sce
	a=sce[["RNA"]]@counts
	a=as.matrix(a)
	a=a[rownames(a)=="GFP",]
	a[a>0]="G"
	a[a==0]="nG"

	sce <- AddMetaData(object = sce, metadata = a, col.name = "gfp_marker")

	#The mt gene IDs
	sce[["percent.mt"]] <- PercentageFeatureSet(sce, features  = c("FBgn0013703","FBgn0013690","FBgn0013710","FBgn0013693","FBgn0013684","FBgn0013695","FBgn0013683","FBgn0013702","FBgn0013679","FBgn0013698","FBgn0013686","FBgn0013708","FBgn0013688","FBgn0262952","FBgn0013696","FBgn0013700","FBgn0013680","FBgn0013709","FBgn0013674","FBgn0013699","FBgn0013675","FBgn0013697","FBgn0013691","FBgn0013673","FBgn0013672","FBgn0013676","FBgn0013694","FBgn0013681","FBgn0013689","FBgn0013704","FBgn0013701","FBgn0013705","FBgn0013692","FBgn0013707","FBgn0013685","FBgn0013678","FBgn0013706","FBgn0013687"))

	pdf(paste0(ids[i],".pdf"),width=4000,height=2000)
	VlnPlot(sce, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
	dev.off()
	#plot1 <- FeatureScatter(sce, feature1 = "nCount_RNA", feature2 = "percent.mt")
	#plot2 <- FeatureScatter(sce, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
	#plot1 + plot2

	sce <- subset(sce, subset = nFeature_RNA > 10 & percent.mt < 10)
	sce <- NormalizeData(sce)

	sce <- FindVariableFeatures(sce, selection.method = "vst", nfeatures = 2000)

	# Identify the 10 most highly variable genes
	#top10 <- head(VariableFeatures(sce), 10)

	# plot variable features with and without labels
	#plot1 <- VariableFeaturePlot(sce)
	#plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
	#plot1 + plot2

	all.genes <- rownames(sce)
	sce <- ScaleData(sce, features = all.genes)

	sce <- RunPCA(sce, features = VariableFeatures(object = sce))
	#print(sce[["pca"]], dims = 1:5, nfeatures = 5)

	#VizDimLoadings(sce, dims = 1:2, reduction = "pca")
	#DimPlot(sce, reduction = "pca")
	#DimHeatmap(sce, dims = 1, cells = 500, balanced = TRUE)
	#DimHeatmap(sce, dims = 1:15, cells = 500, balanced = TRUE)

	sce <- JackStraw(sce, num.replicate = 100)
	sce <- ScoreJackStraw(sce, dims = 1:20)

	pdf(paste0(ids[i],".2.pdf"))
	JackStrawPlot(sce, dims = 1:15)
	dev.off()

	pdf(paste0(ids[i],".3.pdf"))
	ElbowPlot(sce,ndims=30)
	dev.off()

	sce <- FindNeighbors(sce, dims = 1:30)
	# resolution
	sce <- FindClusters(sce, resolution = 0.5)

	#head(Idents(sce), 5)

	sce <- RunUMAP(sce, dims = 1:30)
	#DimPlot(sce, reduction = "umap")

	sce=RunTSNE(sce, dims = 1:30)

	saveRDS(sce,paste0(ids[i],"/",ids[i],".rds"))

	pdf(paste0(ids[i],"/",ids[i],"tSNE1.pdf"),width=5000,height=2000)
	DimPlot(sce, reduction = "tsne", pt.size = 0.5,group.by="gfp_marker")
	dev.off()
	pdf(paste0(ids[i],"/",ids[i],"tSNE2.pdf"),width=5000,height=2000)
	DimPlot(sce, reduction = "tsne", pt.size = 0.5,split.by="gfp_marker")
	dev.off()

	#sce=readRDS(paste0(ids[i],"/",ids[i],".rds"))
		
	sce.markers <- FindAllMarkers(sce, only.pos = TRUE, min.pct = 0, logfc.threshold = 0)
	write.csv(sce.markers,paste0(ids[i],"/",ids[i],".DE_cluster.csv"),quote=FALSE,row.names=FALSE)

	sce\$clus <- Idents(sce)
		
	Idents(sce) <- sce[[]]\$gfp_marker
		
	sce.markers <- FindMarkers(sce, ident.1 = "G", ident.2 = "nG", only.pos = FALSE, min.pct = 0, logfc.threshold = 0)
	write.csv(sce.markers,paste0(ids[i],"/",ids[i],".DE_gfp.csv"),quote=FALSE)
	#Add column name for 1st column.

	# Not used from here
	#sce\$clus <- Idents(sce)

	#Idents(sce) <- sce[[]]\$gfp_marker

	#sce.markers <- FindAllMarkers(sce, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
	#sce.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_logFC)
	#write.csv(sce.markers,paste0(ids[i],"//",ids[i],".DE_gfp.csv"),quote=FALSE,row.names=FALSE)

	rm(sce)
	print(i)

}
```

### Integration

```
library(Seurat)

ids=c("Sample_P","Sample_KC","Sample_KC_12","Sample_12")

for (i in 1:length(ids)) {

	name=paste0(ids[i],"/outs/filtered_feature_bc_matrix")

	seurat_data <- Read10X(data.dir = name)
	seurat_obj <- CreateSeuratObject(counts = seurat_data, 
	                                     project = ids[i])    
    
	a=seurat_obj[["RNA"]]@counts
	a=as.matrix(a)
	a=a[rownames(a)=="GFP",]
	a[a>0]="G"
	a[a==0]="nG"

	seurat_obj <- AddMetaData(object = seurat_obj, metadata = a, col.name = "gfp_marker") 

	seurat_obj$tp <- ids[i]

	seurat_obj <- subset(seurat_obj, subset = nFeature_RNA > 200)
	seurat_obj <- NormalizeData(seurat_obj)
	seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst", nfeatures = 2000)

	assign(ids[i], seurat_obj)
}

data.anchors <- FindIntegrationAnchors(object.list = list(Sample_P,Sample_KC,Sample_KC_12,Sample_12), dims = 1:30)

saveRDS(data.anchors,"data.anchors.rds")

#data.anchors=readRDS("data.anchors.rds")
data.combined <- IntegrateData(anchorset = data.anchors, dims = 1:30)
saveRDS(data.combined,"data.combined.rds")
```

#### Clustering and plots

```
#data.combined=readRDS("data.combined.rds")
DefaultAssay(data.combined) <- "integrated"

data.combined <- ScaleData(data.combined, verbose = FALSE)
data.combined <- RunPCA(data.combined, npcs = 30, verbose = FALSE)
# t-SNE and Clustering
data.combined <- RunUMAP(data.combined, reduction = "pca", dims = 1:30)
data.combined <- FindNeighbors(data.combined, reduction = "pca", dims = 1:30)
data.combined <- FindClusters(data.combined, resolution = 0.5)
data.combined <- RunTSNE(data.combined, reduction = "pca", dims = 1:30)

saveRDS(data.combined,"data.combined2.rds")

pdf("data.combined.tSNE1.pdf",width=2000,height=2000)
DimPlot(data.combined, reduction = "tsne", pt.size = 0.5,group.by="gfp_marker")
dev.off()
pdf("data.combined.tSNE2.pdf",width=4000,height=2000)
DimPlot(data.combined, reduction = "tsne", pt.size = 0.5,split.by="gfp_marker")
dev.off()

new=paste(data.combined[[]]$tp,data.combined[[]]\$gfp_marker,sep="__")
data.combined <- AddMetaData(object = data.combined, metadata = new, col.name = "new") 

pdf("data.combined.tSNE3.pdf",width=3000,height=2000)
DimPlot(data.combined, reduction = "tsne", pt.size = 0.5,group.by="new",cols="glasbey")
dev.off()

pdf("data.combined.tSNE4.pdf",width=6000,height=2000)
DimPlot(data.combined, reduction = "tsne", pt.size = 0.5,group.by="gfp_marker",split.by="tp")
dev.off()

pdf("data.combined.tSNE5.pdf",width=6000,height=2000)
DimPlot(data.combined, reduction = "tsne", pt.size = 0.5,split.by="tp")
dev.off()

pdf("data.combined.tSNE6.pdf",width=6000,height=2000)
DimPlot(data.combined, reduction = "tsne", pt.size = 0.5,split.by="gfp_marker")
dev.off()
```

#### Write DE lists

```
#Use "RNA" assay instead of "integrated". See https://www.biostars.org/p/399789/ and https://github.com/satijalab/seurat/issues/1168.

data.combined=readRDS("data.combined2.rds")
DefaultAssay(data.combined) <- "RNA"
new=paste(data.combined[[]]\$tp,data.combined[[]]\$gfp_marker,sep="__")
data.combined <- AddMetaData(object = data.combined, metadata = new, col.name = "new") 

data.combined\$clus <- Idents(data.combined)
Idents(data.combined) <- data.combined[[]]\$new

data.markers <- FindMarkers(data.combined, ident.1 = "Sample_KC__G", ident.2 = "Sample_P__G", only.pos = FALSE, min.pct = 0, logfc.threshold = 0)
write.csv(data.markers,"data.combined.Sample_KC__G_to_Sample_P__G.csv",quote=FALSE)

data.markers <- FindMarkers(data.combined, ident.1 = "Sample_KC_12__G", ident.2 = "Sample_P__G", only.pos = FALSE, min.pct = 0, logfc.threshold = 0)
write.csv(data.markers,"data.combined.Sample_KC_12__G_to_Sample_P__G.csv",quote=FALSE)

data.markers <- FindMarkers(data.combined, ident.1 = "Sample_12__G", ident.2 = "Sample_P__G", only.pos = FALSE, min.pct = 0, logfc.threshold = 0)
write.csv(data.markers,"data.combined.Sample_12__G_to_Sample_P__G.csv",quote=FALSE)
```

## Annotate genes

Made dmel-coding-GFP-r6.35.annot.txt from dmel-coding-GFP-r6.35.gtf.

```
annot=read.delim("dmel-coding-GFP-r6.35.annot.txt", header=FALSE)

fls=list.files(".","gfp.csv")

for (fl in fls) {
	d1=read.csv(fl)
	fl=gsub("\\.csv","",fl)
	d1=merge(d1,annot,all.x=TRUE,by.x="gene",by.y="V1")
	write.csv(d1,paste0(fl,".annot.csv"),quote=FALSE,row.names=FALSE)
}


#

fls=list.files(".","data.combined.*csv")

for (fl in fls) {
	d1=read.csv(fl)
	fl=gsub("\\.csv","",fl)
	d1=merge(d1,annot,all.x=TRUE,by.x="id",by.y="V1")
	write.csv(d1,paste0(fl,".annot.csv"),quote=FALSE,row.names=FALSE)
}

#

fls=list.files(".","cluster.csv")

for (fl in fls) {
	d1=read.csv(fl)
	fl=gsub("\\.csv","",fl)
	d1=merge(d1,annot,all.x=TRUE,by.x="gene",by.y="V1")
	write.csv(d1,paste0(fl,".annot.csv"),quote=FALSE,row.names=FALSE)
}
```

Order by cluster and log fold change.

## Plots

FBgn0020307     dve
FBgn0259986     nab
FBgn0085424     nub
FBgn0267337     rn
FBgn0005613     Sox15
FBgn0004607     zfh2
FBgn0003866     tsh
FBgn0284084     wg
FBgn0000490     dpp
FBgn0261648     salm
FBgn0000179     omb
FBgn0024250     brk
FBgn0000157     Dll
FBgn0003975     vg

```
library(Seurat)

data.combined=readRDS("data.combined2.rds")

DefaultAssay(data.combined) <- "integrated"

ge="FBgn0261648"

a=data.combined[["RNA"]]@counts
#a=as.matrix(a)
a=a[rownames(a)==ge,]
a[a>0]="E"
a[a==0]="NE"

data.combined <- AddMetaData(object = data.combined, metadata = a, col.name = ge) 

new=paste(data.combined[[]]$tp,data.combined[[]]$gfp_marker,sep="__")
data.combined <- AddMetaData(object = data.combined, metadata = new, col.name = "new") 

png(paste0("data.combined.",ge,".tSNE.png"),width=10000,height=2000,res=300)
DimPlot(data.combined, reduction = "tsne", pt.size = 0.5,group.by=ge,split.by="new")
dev.off()

table(data.combined[[]]$new,data.combined[[]][[ge]])

###

library(Seurat)

ids=c("Sample_P","Sample_KC","Sample_KC_12","Sample_12")

for (i in 1:length(ids)) {

sce=readRDS(paste0(ids[i],"//",ids[i],".rds"))
	
sce$clus <- Idents(sce)


for (ge in c("FBgn0020307","FBgn0259986","FBgn0085424","FBgn0267337","FBgn0005613","FBgn0004607","FBgn0003866","FBgn0284084","FBgn0000490","FBgn0261648","FBgn0024250","FBgn0000157","FBgn0003975")) {
	a=sce[["RNA"]]@counts
	a=as.matrix(a)
	a=a[rownames(a)==ge,]
	a[a>0]="E"
	a[a==0]="NE"
	
	sce <- AddMetaData(object = sce, metadata = a, col.name = ge)

	Idents(sce) <- sce[[]][[ge]]
	
	png(paste0(ids[i],"//",ids[i],".",ge,".tSNE.png"),width=2000,height=2000,res=300)
	print(DimPlot(sce, reduction = "tsne", pt.size = 0.5,split.by="gfp_marker"))
	dev.off()
	}
}
```

#### Plots2

```
library(Seurat)

ids=c("Sample_P","Sample_KC","Sample_KC_12","Sample_12")
	
for (i in 1:length(ids)) {
	sce=readRDS(paste0(ids[i],"/",ids[i],".rds"))
	png(paste0(ids[i],"/",ids[i],"tSNE1.png"),width=2500,height=2000,res=300)
	print(DimPlot(sce, reduction = "tsne", pt.size = 0.5,group.by="gfp_marker"))
	dev.off()
	png(paste0(ids[i],"/",ids[i],"tSNE2.png"),width=4500,height=2000,res=300)
	print(DimPlot(sce, reduction = "tsne", pt.size = 0.5,split.by="gfp_marker"))
	dev.off()
}

#

data.combined=readRDS("data.combined2.rds")

png("data.combined.tSNE1.png",width=2500,height=2000,res=300)
print(DimPlot(data.combined, reduction = "tsne", pt.size = 0.5,group.by="gfp_marker"))
dev.off()
png("data.combined.tSNE2.png",width=4500,height=2000,res=300)
print(DimPlot(data.combined, reduction = "tsne", pt.size = 0.5,split.by="gfp_marker"))
dev.off()
	
new=paste(data.combined[[]]$tp,data.combined[[]]$gfp_marker,sep="__")
data.combined <- AddMetaData(object = data.combined, metadata = new, col.name = "new") 
	
png("data.combined.tSNE3.png",width=3000,height=2000,res=300)
print(DimPlot(data.combined, reduction = "tsne", pt.size = 0.5,group.by="new",cols="glasbey"))
dev.off()

png("data.combined.tSNE4.png",width=6000,height=2000,res=300)
print(DimPlot(data.combined, reduction = "tsne", pt.size = 0.5,group.by="gfp_marker",split.by="tp"))
dev.off()

png("data.combined.tSNE5.png",width=6000,height=2000,res=300)
print(DimPlot(data.combined, reduction = "tsne", pt.size = 0.5,split.by="tp"))
dev.off()

png("data.combined.tSNE6.png",width=6000,height=2000,res=300)
print(DimPlot(data.combined, reduction = "tsne", pt.size = 0.5,split.by="gfp_marker"))
dev.off()
```

## nub+ teashirt-

```
for (i in 1:length(ids)) {
sce=readRDS(paste0(ids[i],"/",ids[i],".rds"))

dta=subset(sce, subset = FBgn0085424 > 0)
#dta=subset(dta, subset = FBgn0003866 == 0)

dta <- RunPCA(dta, features = VariableFeatures(object = dta))
dta <- FindNeighbors(dta, dims = 1:30)
# resolution
dta <- FindClusters(dta, resolution = 0.4)

#head(Idents(dta), 5)

dta <- RunUMAP(dta, dims = 1:30)
#DimPlot(dta, reduction = "umap")

dta=RunTSNE(dta, dims = 1:30)

png(paste0(ids[i],".nub_p.tSNE.png"),width=2500,height=2000,res=300)
print(DimPlot(dta, reduction = "tsne", pt.size = 0.5,split.by="gfp_marker"))
dev.off()

saveRDS(dta,paste0(ids[i],".nub_p.rds"))

dta.markers <- FindAllMarkers(dta, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.25)
write.csv(dta.markers,paste0(ids[i],".nub_p.DE_cluster.csv"),quote=FALSE)

dta$clus <- Idents(dta)
Idents(dta) <- dta[[]]$gfp_marker

dta.markers <- FindAllMarkers(dta, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.25)
write.csv(dta.markers,paste0(ids[i],".nub_p.DE_gfp.csv"),quote=FALSE)

print(i)

}

#

for (i in 1:length(ids)) {
sce=readRDS(paste0(ids[i],"/",ids[i],".rds"))

dta=subset(sce, subset = FBgn0085424 > 0)
dta=subset(dta, subset = FBgn0003866 == 0)

dta <- RunPCA(dta, features = VariableFeatures(object = dta))
dta <- FindNeighbors(dta, dims = 1:30)
# resolution
dta <- FindClusters(dta, resolution = 0.4)

#head(Idents(dta), 5)

dta <- RunUMAP(dta, dims = 1:30)
#DimPlot(dta, reduction = "umap")

dta=RunTSNE(dta, dims = 1:30)

png(paste0(ids[i],".nub_p.teashirt_n.tSNE.png"),width=2500,height=2000,res=300)
print(DimPlot(dta, reduction = "tsne", pt.size = 0.5,split.by="gfp_marker"))
dev.off()

saveRDS(dta,paste0(ids[i],".nub_p.teashirt_n.rds"))

dta.markers <- FindAllMarkers(dta, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.25)
write.csv(dta.markers,paste0(ids[i],".nub_p.teashirt_n.DE_cluster.csv"),quote=FALSE)

dta$clus <- Idents(dta)
Idents(dta) <- dta[[]]$gfp_marker

dta.markers <- FindAllMarkers(dta, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.25)
write.csv(dta.markers,paste0(ids[i],".nub_p.teashirt_n.DE_gfp.csv"),quote=FALSE)

print(i)

}

#

annot=read.delim("dmel-coding-GFP-r6.35.annot.txt", header=FALSE)

fls=list.files(".","nub_p.*csv")

for (fl in fls) {
d1=read.csv(fl,row.names=1)
fl=gsub("\\.csv","",fl)
d1=merge(d1,annot,all.x=TRUE,by.x="gene",by.y="V1")
write.csv(d1,paste0(fl,".annot.csv"),quote=FALSE,row.names=FALSE)
}

###

library(Seurat)

setwd("D://scrna")

ids=c("Sample_P","Sample_KC","Sample_KC_12","Sample_12")

for (i in 1:length(ids)) {
dta=readRDS(paste0(ids[i],".nub_p.rds"))

Idents(dta) <- dta[[]]$gfp_marker

dta.markers <- FindAllMarkers(dta, only.pos = TRUE, min.pct = 0, logfc.threshold = 0)
write.csv(dta.markers,paste0(ids[i],".nub_p.DE_gfp0.csv"),quote=FALSE)

print(i)

}

annot=read.delim("dmel-coding-GFP-r6.35.annot.txt", header=FALSE)

fls=list.files(".","nub_p.DE_gfp0.csv")

for (fl in fls) {
d1=read.csv(fl,row.names=1)
fl=gsub("\\.csv","",fl)
d1=merge(d1,annot,all.x=TRUE,by.x="gene",by.y="V1")
write.csv(d1,paste0(fl,".annot.csv"),quote=FALSE,row.names=FALSE)
}
```

## Clusters in the combined dataset for nub+ teashirt-

```
library(Seurat)

sce=readRDS("data.combined2.rds")
a=sce[['RNA']]@counts
b=a[rownames(a)=="FBgn0085424",]
b[b>0]="Pos"
b[b==0]="Neg"
sce <- AddMetaData(object = sce, metadata = b, col.name = "nub") 
new=paste(sce[[]]$tp,sce[[]]$gfp_marker,sep="__")
sce <- AddMetaData(object = sce, metadata = new, col.name = "new")

dta=subset(sce, subset = (nub == "Pos"))
#dta=subset(dta, subset = FBgn0003866 == 0)

dta <- RunPCA(dta, features = VariableFeatures(object = dta))
dta <- FindNeighbors(dta, dims = 1:30)
# resolution
dta <- FindClusters(dta, resolution = 0.4)

#head(Idents(dta), 5)

dta <- RunUMAP(dta, dims = 1:30)
#DimPlot(dta, reduction = "umap")

dta=RunTSNE(dta, dims = 1:30)

saveRDS(dta,"combined.nub_p.rds")

dta.markers <- FindAllMarkers(dta, only.pos = TRUE, min.pct = 0, logfc.threshold = 0)
write.csv(dta.markers,"combined.nub_p.DE_cluster.csv",quote=FALSE)

#dta$clus <- Idents(dta)
#Idents(dta) <- dta[[]]$new

#png("combined.nub_p.tSNE.png",width=16000,height=2000,res=300)
#png("combined.nub_p.tSNE.png")
#print(DimPlot(dta, reduction = "tsne", pt.size = 0.5,split.by="new"))
#dev.off()
	
##

sce=readRDS("data.combined2.rds")
a=sce[['RNA']]@counts

b=a[rownames(a)=="FBgn0085424",]
b[b>0]="Pos"
b[b==0]="Neg"
sce <- AddMetaData(object = sce, metadata = b, col.name = "nub") 

b=a[rownames(a)=="FBgn0003866",]
b[b>0]="Pos"
b[b==0]="Neg"
sce <- AddMetaData(object = sce, metadata = b, col.name = "teashirt") 

new=paste(sce[[]]$tp,sce[[]]$gfp_marker,sep="__")
sce <- AddMetaData(object = sce, metadata = new, col.name = "new")

dta=subset(sce, subset = (nub == "Pos"))
dta=subset(dta, subset = (teashirt == "Neg"))

dta <- RunPCA(dta, features = VariableFeatures(object = dta))
dta <- FindNeighbors(dta, dims = 1:30)
# resolution
dta <- FindClusters(dta, resolution = 0.4)

#head(Idents(dta), 5)

dta <- RunUMAP(dta, dims = 1:30)
#DimPlot(dta, reduction = "umap")

dta=RunTSNE(dta, dims = 1:30)

saveRDS(dta,"combined.nub_p.ts_n.rds")

dta.markers <- FindAllMarkers(dta, only.pos = TRUE, min.pct = 0, logfc.threshold = 0)
write.csv(dta.markers,"combined.nub_p.ts_n.DE_cluster.csv",quote=FALSE)

#dta$clus <- Idents(dta)
#Idents(dta) <- dta[[]]$new

#png("combined.nub_p.ts_n.tSNE.png",width=16000,height=2000,res=300)
#png("combined.nub_p.ts_n.tSNE.png")
#print(DimPlot(dta, reduction = "tsne", pt.size = 0.5,split.by="new"))
#dev.off()
```

```
library(Seurat)
dta=readRDS("combined.nub_p.rds")
png("combined.nub_p.tSNE.png",width=16000,height=2000,res=300)
print(DimPlot(dta, reduction = "tsne", pt.size = 0.5,split.by="new"))
dev.off()

dta=readRDS("combined.nub_p.ts_n.rds")
png("combined.nub_p.ts_n.tSNE.png",width=16000,height=2000,res=300)
print(DimPlot(dta, reduction = "tsne", pt.size = 0.5,split.by="new"))
dev.off()
```

### More plots

FBgn0003117	pnr
FBgn0085424 nub
FBgn0003866	tsh
FBgn0001235	hth
FBgn0004607	zfh2
FBgn0010453	Wnt4
	
```
library(Seurat)

dta=readRDS("combined.nub_p.rds")

data.frame(table(dta[[]]$seurat_clusters))
data.frame(table(dta[[]]$new))

for (ge in c("FBgn0003117","FBgn0085424","FBgn0003866","FBgn0001235","FBgn0004607","FBgn0010453")) {

	a=dta[['RNA']]@counts
	b=a[rownames(a)==ge,]
	dta <- AddMetaData(object = dta, metadata = b, col.name = ge) 
	png(paste0(ge,".nub_p.png"),width=18000,height=2000,res=300)
	print(FeaturePlot(dta, features = ge,reduction = "tsne",split.by="new",cols=c("turquoise","blue")))
	dev.off()
	png(paste0(ge,".nub_p.legend.png"),width=2000,height=2000,res=300)
	print(FeaturePlot(dta, features = ge,reduction = "tsne",cols=c("turquoise","blue")))
	dev.off()
}

##

dta=readRDS("combined.nub_p.ts_n.rds")

data.frame(table(dta[[]]$seurat_clusters))
data.frame(table(dta[[]]$new))

for (ge in c("FBgn0003117","FBgn0085424","FBgn0003866","FBgn0001235","FBgn0004607","FBgn0010453")) {

	a=dta[['RNA']]@counts
	b=a[rownames(a)==ge,]
	dta <- AddMetaData(object = dta, metadata = b, col.name = ge) 
	png(paste0(ge,".nub_p.tsh_n.png"),width=18000,height=2000,res=300)
	print(FeaturePlot(dta, features = ge,reduction = "tsne",split.by="new",cols=c("turquoise","blue")))
	dev.off()
	png(paste0(ge,".nub_p.tsh_n.legend.png"),width=2000,height=2000,res=300)
	print(FeaturePlot(dta, features = ge,reduction = "tsne",cols=c("turquoise","blue")))
	dev.off()
}
```

### Annotation
	
```
annot=read.delim("dmel-coding-GFP-r6.35.annot.txt", header=FALSE)

#fl="combined.nub_p.ts_n.DE_cluster.csv"
fl="combined.nub_p.DE_cluster.csv"
d1=read.csv(fl)
fl=gsub("\\.csv","",fl)
d1=merge(d1,annot,all.x=TRUE,by.x="gene",by.y="V1")
write.csv(d1,paste0(fl,".annot.csv"),quote=FALSE,row.names=FALSE)


for (fl in c("data.combined.Sample_12__G_to_Sample_P__G.csv","data.combined.Sample_KC_12__G_to_Sample_P__G.csv","data.combined.Sample_KC__G_to_Sample_P__G.csv")) {
	d1=read.csv(fl)
	fl=gsub("\\.csv","",fl)
	d1=merge(d1,annot,all.x=TRUE,by.x="X",by.y="V1")
	write.csv(d1,paste0(fl,".annot.csv"),quote=FALSE,row.names=FALSE)
}
```

