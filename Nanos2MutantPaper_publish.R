# Rscript for generating data objects and figures in Denner et al. 2023 BioRxiV
# Setup ----
# load libraries 
library(easypackages)
libraries("Seurat", "Matrix", "readxl","RColorBrewer",'Rmagic',
          'patchwork','dplyr','viridis','ggplot2','pals','SeuratWrappers')

## load the data
load (file='nanos2.paper.revised.RData') #workspace with necessary objects.

# Sections ----
# two sections: one to analyze the data, one to generate the figures from the paper
generate_data= F
generate_figures = T

## Generate datasets ----
if (generate_data)
{
data1=Tissue.NanosMutant
levels(SetIdent(data1,value='orig.ident'))

### Run Seurat pipeline All Data----
{
  #### normalize the dataset (merged)
  data1 <- NormalizeData(data1, scale.factor = 5000) 

  #### calculate variable genes
  #select 500 variable features from each library and use this set.
  {
    list=  NULL
    vargenelist <- SplitObject(data1, split.by = "orig.ident")
    for (i in 1:length(vargenelist)) {
      vargenelist[[i]] <- FindVariableFeatures(vargenelist[[i]],nfeatures = 500, verbose = FALSE)
    }
    for (i in 1:length(vargenelist)) {
      x <- vargenelist[[i]]@assays$RNA@var.features
      list=c(list,x)}
    list=unique(list)
    length(list)
    
    data1@assays$RNA@var.features = list
    
  }
  data1 <- ScaleData(data1,split.by = 'orig.ident')#
  
  #### Dimensional reductions
  data1 <- RunPCA(data1, pcs.compute = 50)
  PCAPlot(data1,group.by='orig.ident')
  ElbowPlot(object = data1, ndims = 50)
  
  # set dimensions
  d=as.integer(which(data1@reductions$pca@stdev>2))
  d=1:30 #used in atlas, not in nanos.mutant merge
  
  #UMAP
  data1 <- RunUMAP(data1, dims = d,
                   reduction = 'pca',
                   n.neighbors = 20L,
                   spread =1, min.dist = 0.6, local.connectivity = 1,
                   seed.use = 0)
  
  ### Find Clusters
  data1 <- FindNeighbors(object = data1,reduction ="pca",dims = d,
                         annoy.metric = 'cosine',
                         k.param = 30)
  
  data1 <- FindClusters(object = data1,resolution = 0.2)
  data1 <- BuildClusterTree(object = data1, reorder = TRUE,reorder.numeric = T)
}

## label clusters ----
#run a semi-automated script for updating cluster names based on marker gene expression
{
  clusterNames<- read_excel("ClusterAnnotations_Nv2.xlsx") #From Nv2 Atlas paper: Supplemental data 2
  goi = clusterNames$Marker

    #assign cluster ID to the individual libraries
  data1<-ScaleData(data1,features = goi, split.by = 'orig.ident') #works with scaled gene values, so must scale all of these first (in case they are not within the VarFeatures set!)
  
  cl <-length(levels(data1@active.ident)) #how many clusters
  C.suffix <-seq(1:cl) #the numbers
  
  g=length(goi)
  clName = vector()
  m=matrix(0L,g,cl) #generate a matrix of values of each cluster for each gene:
  for (j in 1:cl) #for each cluster set
  {
    for (i in 1:g) #for each gene
      m[i,j]=mean(data1@assays$RNA@scale.data[goi[i],WhichCells(data1,idents = C.suffix[j])]) #average scaled value for cluster/gene
    clName[j]=as.integer(which.max(m[,j])) #choose the highest value
  }
  sort(clName) #this prints which IDs were selected; sanity check that it did what you wanted... should have no unwanted replicated IDs
  clusterNames$ID.Alldata[clName] #you can also look at the names to see where the problem might have been
  
  #use the desired order from the spreadsheet to re-order the clusters:
  #first order the identities..
  data1@active.ident = factor(data1@active.ident,
                              levels(data1@active.ident)[order(clName)])
  #set the names to your IDs; once the gene list works, it is really easy to then update your cluster names with your excel spreadsheet.
  levels(data1@active.ident) = clusterNames$ID.Alldata[clName][order(clName)]
  #save the IDs in metadata:
  data1@meta.data$IDs = data1@active.ident
}

Tissue.NanosMutant=data1

## Run Seurat pipeline PGC subset ----
data1=subset(Tissue.NanosMutant,idents = levels(Tissue.NanosMutant)[c(1,2)])
{
  data1@active.assay = 'RNA'
  data1 <- FindVariableFeatures(data1,nfeatures = 2000)#,
  
  coi = as.vector(data1@assays$RNA@counts@Dimnames[[2]])
  t=ScaleData(Tissue.NanosMutant,
              split.by = 'orig.ident', features = data1@assays$RNA@var.features)
  data1@assays$RNA@scale.data = t@assays$RNA@scale.data[,coi]
  
  # run different reduction algorythms:
  # #PCA
  data1 <- RunPCA(data1, pcs.compute = 50)
  ElbowPlot(object = data1, ndims = 50)
  
  # set dimensions
  d=as.integer(which(data1@reductions$pca@stdev>2))
  d# 
  data1 <- RunUMAP(data1, reduction ="pca", 
                   n.neighbors = 20L,spread = 1,
                   dims = d,reduction.name ='umap',
                   reduction.key ='umap',min.dist = 0.6,
                   local.connectivity = 1)#
  data1 <- FindNeighbors(object = data1,
                         reduction ="pca",dims = d,
                         nn.method = 'annoy',
                         annoy.metric = 'cosine',
                         k.param = 10)
  data1 <- FindClusters(object = data1,resolution = 0.1,random.seed = 0)
  # #look at relationship between clusters
  data1 <- BuildClusterTree(object = data1, reorder = TRUE,
                            reorder.numeric = TRUE)#, dims = c(1:d))
  data1$IDs.coarse = data1@active.ident
  
}#

##label the clusters: ----
names=c('primary oocytes','spermatozoa','maturing sperm 1','maturing sperm 2','spermatagonia','pSC','PGCs')
levels(data1@active.ident) = names
data1@active.ident = factor(data1@active.ident,
                            levels(data1@active.ident)[c(6,7,1,5,3,4,2)])
DimPlot(data1,label = T,cells.highlight =WhichCells(SetIdent(Tissue.NanosMutant,value='orig.ident'),ident='nanos.mutant') )
nanos.PGCs=data1

# add the geneotype information to the objects:
Tissue.NanosMutant$genotype = 'wildtype'
Tissue.NanosMutant$genotype[WhichCells(SetIdent(Tissue.NanosMutant,value='orig.ident'),ident='nanos.mutant')]='nanos2.mutant'

nanos.PGCs$genotype = 'wildtype'
nanos.PGCs$genotype[WhichCells(SetIdent(nanos.PGCs,value='orig.ident'),ident='nanos.mutant')]='nanos2.mutant'

#import the clustering of PGC/pSC into the full dataset:
coi=nanos.PGCs@assays$RNA@counts@Dimnames[[2]]
Tissue.NanosMutant <- SetIdent(Tissue.NanosMutant, cells = coi, value = nanos.PGCs@active.ident[coi])

# clean up the workspace
rm(data1,cluster.plot,library.plot,dist.clust2,d,g,coi,i,j,list,names,x,cl,C.suffix,clName,goi,t,vargenelist,m,ids.cluster.library,dist.clust,Alldata.tissue.Nv2,annotations,nanos.mutant,clust.cp,gene.cp.dark)
}

## Figures paper ----
if (generate_figures)
{
#set the colour palet
clust.cp=clust.cp.graded[2:(length(levels(Tissue.NanosMutant))+1)]

FigA=DimPlot(Tissue.NanosMutant, label = T,label.size = 5, repel = T,reduction = 'umap',
             cols = clust.cp,raster=T,pt.size = 2)+NoLegend()+NoAxes()+labs(tag='A',title='Adult cell populations')
FigA

data1=Tissue.NanosMutant
ids.cluster.library = as.data.frame(table(Idents(data1), data1@meta.data$genotype))
colnames(ids.cluster.library) = c('ID','Genotype','CellCount')

#barplot of cluster identities in each library:

FigB=
  ggplot(ids.cluster.library, aes(fill=ID, y= CellCount,
                                  x=Genotype)) +
  geom_bar(mapping =aes(fill=ID, y= (CellCount),
                        x=(Genotype)),
           position="fill", stat="identity", width = 0.5)+
  scale_fill_manual(values = clust.cp)+
  theme(axis.text.x = element_text(#face="bold", color="#993333", 
    size=8, angle=-45,hjust=0,vjust = 0.5))+
  geom_area(mapping =aes(fill=ID, y= (CellCount),
                         x=as.integer(Genotype)),
            position="fill", stat="identity",alpha=0.2 , size=.5, colour="white") +
  geom_bar(mapping =aes(fill=ID, y= (CellCount),#this re-plots the bars over the area
                        x=(Genotype)),
           position="fill", stat="identity", width = 0.5)+
  ggtitle("Distribution of cell types in time and space")+labs(tag='B')

FigB

FigC = DimPlot(Tissue.NanosMutant,cells.highlight = WhichCells(SetIdent(Tissue.NanosMutant,value='orig.ident'),ident='nanos.mutant'),raster=T,pt.size = 2)&NoAxes()&NoLegend()&labs(tag='C',title='Distribution of cells from mutant',subtitle = 'Nanos2 Mutant')

FigC

pdf(file = 'figureABCrevised.pdf',width = 24,height = 8)
FigA+FigB+FigC+FigD+plot_layout(ncol=4)
dev.off()

pdf(file = 'figureD.pdf',width = 24,height = 6,onefile = F)
DotPlot(Tissue.NanosMutant,'RNA',rev(c('Nanos2','Nanos1','PIWL1-like-5',vasa2,'SoxC','SoxB.2')),scale.by = 'size',split.by='genotype',cols=c('black','darkred'),col.min = 0 , idents = levels(Tissue.NanosMutant)[c(1:7,9:13,20)],dot.min = 0.05,dot.scale = 15)& theme(panel.grid = element_line(size=0.002,colour = 'grey90'))&RotatedAxis()&labs(tag='D')&coord_flip()&theme(legend.position = 'bottom')
dev.off()
}