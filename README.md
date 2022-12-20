# Hypergeo_SingleCell_Markers
R function for hypergeometric distribution test of cluster marker overlap.

## Overview
Single cell sequencing processing produces clusters of cells, from which marker genes can be defined. Datasets that use the same feature genes can be compared by investigating the overlap of markers from clusters. Using a hypergeometric distribution test this R function can take two Seurat objects, detect markers and compare in an all versus all cluster approach.

This tool maybe helpful for the detection of similar/different cell types/clusters between samples.

## R code with an example

To run the function you need Seurat objects that have clusters determined so that markers can be defined.

An example of how to use the below function.

```{r}
# load packages
library(Seurat)
library(SeuratData)
library(sctransform)
library(sctransform)

# get test Seurat data
SeuratData::InstallData("pbmcsca")

# SEURAT PROCESSING WITH SCTransform
dims   = 50
thresh = 1.3
res    = 1.5

#normalise
pbmcsca <- SCTransform(pbmcsca,
  method = "glmGamPoi", 
  do.scale = TRUE,
  variable.features.n = NULL, 
  variable.features.rv.th = thresh, 
  vars.to.regress = c("nCount_RNA","nFeature_RNA"),
  verbose = FALSE)  

#PCA
pbmcsca <- RunPCA(pbmcsca,
  verbose = FALSE,
  dims = 1:dims)

#UMAP
pbmcsca <- RunUMAP(pbmcsca,
  reduction = "pca",
  dims = 1:dims,
  umap.method="umap-learn",
  n.neighbors = 5L,
  min.dist = 0.3,
  metric = "correlation",
  seed.use = 42)

#Neighbors
pbmcsca <- FindNeighbors(pbmcsca, 
  dims = 1:dims, 
  verbose = FALSE, 
  n.trees = 50, 
  k.param = 200)

#clusters
pbmcsca <- FindClusters(pbmcsca,
  verbose = FALSE,
  algorithm=1,
  resolution = res)                         

pbmcsca <- PrepSCTFindMarkers(pbmcsca, 
  assay = "SCT",
  verbose=FALSE)

#Plot UMAP
DimPlot(pbmcp, 
  label=T)
```
UMAP plot for the pbmcsca data.

<img src="https://github.com/djleite/Hypergeo_SingleCell_Markers/blob/main/UMAP_pbmcsca.png" width="400">


The function can be run to compare that test dataset to itself. Clusters with similar markers have low pvalues.

```{r}
# input Seurat object, the percentage of cells and the pvalue threshold for marker detection
p <- HGD(pbmcsca,pbmcsca,0.1,1e-3)
p+annotate("text", y=0.99,x=1,fontface='plain',hjust=1, label = "log(p-value)")
```

<img src="https://github.com/djleite/Hypergeo_SingleCell_Markers/blob/main/HGD_pbmcsca.png" width="400">



## The R function

```{r}
# load packages
library(Seurat)
library(pheatmap)
library(ggplotify)
library(ggplot2)

# GET MARKERS FROM SEURAT OBJECT FOR HYPERGEOMETRIC/FISHERS-EXACT TEST
# <indata> needs to be a seurat object with neighbors and clusters
# the pct and return thresh hold are adjustable

seurat_markers_for_HGD <- function(indata,pct,rthresh) {
    indata.markers = FindAllMarkers(object = indata,
                                        test.use = "wilcox", 
                                        only.pos=TRUE, 
                                        min.pct = pct, 
                                        return.thresh = rthresh,
                                        verbose=T)
    
    # Return in unique marker list to test the overlap
    return(as.data.frame(cbind(as.numeric(indata.markers[["cluster"]])-1,indata.markers[["gene"]])))
}

# FISHER'S EXACT TEST FOR CLUSTER MARKERS
# <*.obj> needs to be a seurat object with neighbors and clusters
# the pct and return thresh hold are adjustable

HGD <- function(x.obj,y.obj,pct,rthresh){
    
    # get markers
    x.markers <- seurat_markers_for_HGD(x.obj,pct,rthresh)
    y.markers <- seurat_markers_for_HGD(y.obj,pct,rthresh)
    
    # turnoff scientific numbers
    options(scipen = 999)
    
    # prepare dataframe
    df = data.frame(matrix(ncol = length(unique(x.markers$V1)), nrow = length(unique(y.markers$V1))))
    colnames(df) <- seq(0,length(unique(x.markers$V1))-1)
    
    # for each x cluster
    for (x in unique(x.markers$V1)){
        
        # get markers for each x cluster
        clx.markers = x.markers[x.markers$V1 == x,][,2]
        dfx <- numeric()

        # for each y cluster
        for (y in unique(y.markers$V1)){
            
            # get markers for each y cluster
            cly.markers = y.markers[y.markers$V1 == y,][,2]
            
            intersect <- length(intersect(clx.markers,cly.markers)) #intersect between x and y markers
            xlen      <- length(clx.markers)                        #size of x markers
            ylen      <- length(cly.markers)                        #size of y markers
            tsize     <- length(unique(c(rownames(data.matrix(rowSums(x.obj[['RNA']]@counts) != 0)),rownames(data.matrix(rowSums(y.obj[['RNA']]@counts) != 0)))))-length(x.markers)

            xy <- sum(dhyper(intersect:xlen, ylen, tsize, xlen)) #hypergeometric distribution test
            dfx <- c(dfx,xy) 


        }
        
        #dfx.adjusted <- log(p.adjust(dfx, method='bonferroni'))
        df[x] <- dfx#.adjusted

    }
    
    df <- log(matrix(p.adjust(as.vector(as.matrix(df)), method='bonferroni'),ncol=length(unique(x.markers$V1))))

    # add cluster numbers
    rownames(df) <- seq(0,length(unique(y.markers$V1))-1)
    
    # dataframe adjustments
    df[sapply(df, is.infinite)] <- NA                 # replace Inf values with min of df 
    df[sapply(df, is.na)] <- min(df,na.rm = TRUE)     # replace Inf values with min of df 
    df[df == 0] <- NA                                 # replace zeros with NAs
    
    # USE THIS TO FILTER RESULTS 
    df[df > log(1e-5)] <- NA                          # replace values below p-value with NAs

    
    # return dataframe. You can comment out the plotting function and return the dataframe.
    #return(df)
    
    # turnon scientific numbers
    options(scipen = 0)
    
    #PLOT HEATMAP
    # colour palette
    paletteLength <- 50
    myColors      <- colorRampPalette(c("darkred","red","orange","yellow"),interpolate = c("spline"))(paletteLength)

    # labels for legend
    dfmin <- as.integer(min(df,na.rm = TRUE))
    dfmax <- as.integer(max(df,na.rm = TRUE)-1)
    dfq <- as.integer((min(df,na.rm = TRUE))/4)
    
    # heatmap plot
    gg <- as.ggplot(pheatmap(data.matrix(df),
             annotation_names_row=TRUE,
             annotation_names_col=TRUE,
             main=' ',
             cexRow=10, 
             cexCol=10, 
             cluster_cols=F, 
             cluster_rows=F,
             color=myColors,
             legend=T,
             legend_labels = c(dfmin,dfq,dfq*2,dfq*3,dfmax),
             annotation_legend=T,
             legend_breaks=c(dfmin,dfq,dfq*2,dfq*3,dfmax),
             na_col='grey'))
    
    # return the plot
    plot <- gg+annotate("text", y=0.99,x=1,fontface='plain',hjust=1, label = "log(p-value)")
    return(plot)
}
```
