## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
library(CoNI)

## ----install_dependencies-----------------------------------------------------
# dependencies<-c("igraph","doParallel","cocor","tidyverse","foreach","ggrepel","gplots","gridExtra","plyr","ppcor","tidyr","Hmisc")
# 
# `%notin%`<-Negate(`%in%`)
# for(package in dependencies){
#   if(package %notin% rownames(installed.packages())){
#     install.packages(package,dependencies = TRUE)
#   }
# }

## ---- eval=FALSE, echo=TRUE---------------------------------------------------
#  # if (!requireNamespace("BiocManager", quietly = TRUE))
#  #     install.packages("BiocManager")
#  # BiocManager::install("genefilter")

## ----Download_data------------------------------------------------------------
#Chow Data
data(Chow_MetaboliteData) #Metabolite data
data(Chow_GeneExpData) #Gene expression data

#HFD data
data(HFD_MetaboliteData) #Metabolite data
data(HFD_GeneExpData) #Gene expression data

## ----matchRows----------------------------------------------------------------
#Match rownames both omics data
rownames(Chow_MetaboliteData)<-rownames(Chow_GeneExpData)
rownames(HFD_MetaboliteData)<-rownames(HFD_GeneExpData)
#Shorten names
Chow_metabo<-Chow_MetaboliteData
Chow_gene<-Chow_GeneExpData
HFD_metabo<-HFD_MetaboliteData
HFD_gene<-HFD_GeneExpData

## ----CoNI_chow,eval=FALSE, echo=TRUE------------------------------------------
#  #Run for Chow
#  # CoNIResults_Chow<-CoNI(edgeD = Chow_gene,vertexD = Chow_metabo,
#  #                        saveRaw = FALSE,filter_highVarianceEdge = TRUE,correlateDFs=TRUE,
#  #                        padjustvertexD = FALSE, split_number = 200,
#  #                        outputDir = "./Chow/",outputName = "CoNIChow",
#  #                        splitedgeD = TRUE,numCores = 2,onlySgRes = TRUE)

## ----load_CoNIResults_chow----------------------------------------------------
#Load chow results
data(CoNIResults_Chow)

## ----CoNI_hfd,eval=FALSE, echo=TRUE-------------------------------------------
#  #Run for HFD
#  # CoNIResults_HFD<-CoNI(edgeD = HFD_gene,vertexD = HFD_metabo,
#  #                        saveRaw = FALSE,filter_highVarianceEdge = TRUE,
#  #                        padjustvertexD = FALSE, split_number = 200,correlateDFs=TRUE,
#  #                        outputDir = "./HFD/",outputName = "CoNIHFD",
#  #                        splitedgeD = TRUE,numCores = 2,onlySgRes = TRUE)

## ----load_CoNIResults_hfd-----------------------------------------------------
#Load HFD results
data(CoNIResults_HFD)

## ----read_metabolite_annotation-----------------------------------------------
#Read Annotation table
data(MetaboliteAnnotation)
MetaboliteAnnotation<-assign_colorsAnnotation(MetaboliteAnnotation,col="Class")

## ----network_chow-------------------------------------------------------------
#Generate Network
ChowNetwork<-generate_network(ResultsCoNI = CoNIResults_Chow,
                             colorVertexTable = MetaboliteAnnotation,
                             outputDir = "./",
                             outputFileName = "Chow")

## ----network_chow_class-------------------------------------------------------
#Generate Network Chow
ChowNetworkWithClass<-generate_network(ResultsCoNI = CoNIResults_Chow,
                             colorVertexTable = MetaboliteAnnotation,
                             outputDir = "./",
                             outputFileName = "Chow",
                             Class = MetaboliteAnnotation)


## ----network_hfd_class--------------------------------------------------------
#Generate Network HFD
HFDNetworkWithClass<-generate_network(ResultsCoNI = CoNIResults_HFD,
                             colorVertexTable = MetaboliteAnnotation,
                             outputDir = "./",
                             outputFileName = "HFD",
                             Class = MetaboliteAnnotation)

## ----network_stats_chow-------------------------------------------------------
library(knitr)
library(kableExtra)
kable(NetStats(Network = ChowNetworkWithClass),caption="Network statistics Chow") %>% kable_styling(full_width = FALSE)

## ----network_stats_hfd--------------------------------------------------------
kable(NetStats(Network = HFDNetworkWithClass),caption="Network statistics HFD")  %>% kable_styling(full_width = FALSE)

## -----------------------------------------------------------------------------
library(igraph)
coordinates = layout_with_fr(ChowNetworkWithClass) #define layout

## ----spectral-----------------------------------------------------------------
Spectral = cluster_leading_eigen(ChowNetworkWithClass)
#See membership for the nodes
Spectral$membership

## ----spectral_image,fig.show='hide', out.width = '100%', out.height = '100%', fig.width = 7, fig.heigt = 7----
#Plot the network
plot(ChowNetworkWithClass, vertex.color=membership(Spectral), layout=coordinates)

## ----greedy-------------------------------------------------------------------
greedy = cluster_fast_greedy(ChowNetworkWithClass)
#See membership for the nodes
greedy$membership

## ----greedy_image,fig.show='hide',out.width = '100%', out.height = '100%', fig.width = 7, fig.heigt = 7----
#Plot the network
plot(ChowNetworkWithClass, vertex.color=membership(greedy), layout=coordinates)

## ----betweenness--------------------------------------------------------------
betweenness = cluster_edge_betweenness(ChowNetworkWithClass,weights=NULL)
#See membership for the nodes
betweenness$membership

## ----betweenness_image,fig.show='hide',out.width = '100%', out.height = '100%', fig.width = 7, fig.heigt = 7----
#Plot the network
plot(ChowNetworkWithClass, vertex.color=membership(betweenness), layout=coordinates)

## ----local_controlling_genes_Chow---------------------------------------------
#Chow
#Get results binomial test 
LCGenes_ResultsBinomialTable_Chow<- find_localControllingFeatures(ResultsCoNI = CoNIResults_Chow,network = ChowNetworkWithClass)
#Get list local controlling genes
LCGenes_Chow<-as.character(unique(LCGenes_ResultsBinomialTable_Chow$edgeFeatures))

## ----local_controlling_genes_HFD----------------------------------------------
#HFD
#Get results binomial test
LCGenes_ResultsBinomialTable_HFD<- find_localControllingFeatures(ResultsCoNI = CoNIResults_HFD,network = HFDNetworkWithClass)
#Get list local controlling genes
LCGenes_HFD<-as.character(unique(LCGenes_ResultsBinomialTable_HFD$edgeFeatures))

## ----local_controlling_tables-------------------------------------------------
#Chow
#Generate table local controlling genes
TableLCFsChow<-tableLCFs_VFs(CoNIResults = CoNIResults_Chow,LCFs = LCGenes_Chow)
#Show first two rows
TableLCFChowk<-kable(TableLCFsChow[1:2,],caption="Table Local Controlling Genes") %>% kable_styling(full_width = FALSE)

#HFD
#Generate table local controlling genes
TableLCFsHFD<-tableLCFs_VFs(CoNIResults = CoNIResults_HFD,LCFs = LCGenes_HFD)
#Show first two rows
TableLCFHFDk<-kable(TableLCFsHFD[1:2,],caption="Table Local Controlling Genes") %>% kable_styling(full_width = FALSE)
TableLCFHFDk

## ----gene_Magnitude,warning=FALSE,message=FALSE-------------------------------
Top10Chow<-top_n_LF_byMagnitude(CoNIResults_Chow,topn = 10)
Top10HFD<-top_n_LF_byMagnitude(CoNIResults_HFD,topn = 10)
head(Top10HFD[,c(1:3,ncol(Top10HFD))])

## ----CorvsPcor_oneCombination,fig.show='hold',warning=FALSE,message=FALSE,out.width="50%",fig.align='center'----
plotPcorvsCor(ResultsCoNI = CoNIResults_HFD,edgeFeature = "Lpin2",
              vertexFeatures = c("PC.ae.C42.2","PC.ae.C42.0"),
              vertexD = HFD_metabo,edgeD = HFD_gene,
              label_edgeFeature = "Gene",plot_to_screen = TRUE,
              outputDir = "./")


## ----CorvsPcor_AllCombinations,fig.show='hold',warning=FALSE,message=FALSE,out.width="40%",fig.show='hide'----
plotPcorvsCor(ResultsCoNI = CoNIResults_HFD,edgeFeature = "Lpin2",
              vertexD = HFD_metabo,edgeD = HFD_gene,
              label_edgeFeature = "Gene",plot_to_screen = TRUE,
              outputDir = "./")


## ----bipartite_graphs---------------------------------------------------------
#Chow
ChowBipartiteGraph<-createBipartiteGraph("./TableForNetwork_Chow.csv",MetaboliteAnnotation)
#Save bipartite graph
write.graph(ChowBipartiteGraph,file="./Chow_bipartite.graphml",format="graphml")

#HFD
HFDBipartiteGraph<-createBipartiteGraph("./TableForNetwork_HFD.csv",MetaboliteAnnotation)
#Save bipartite graph
write.graph(HFDBipartiteGraph,file="./HFD_bipartite.graphml",format="graphml")

## -----------------------------------------------------------------------------
Chow_HypergraphIncidenceM<-createBipartiteGraph("./TableForNetwork_Chow.csv",MetaboliteAnnotation,incidenceMatrix = TRUE)

## ----triplet_comparison-------------------------------------------------------
Compare_Triplets(Treat1_path = "./TableForNetwork_Chow.csv",
                 Treat2_path = "./TableForNetwork_HFD.csv",
                 OutputName = "Shared_triplets_HFDvsChow.csv")


## -----------------------------------------------------------------------------
(LCP_sharedGene_HFDvsChow<-Compare_VertexClasses_sharedEdgeFeatures(
                          Treat1_path = "./TableForNetwork_HFD.csv",
                          Treat2_path = "./TableForNetwork_Chow.csv",
                          OutputName = "Comparison_LipidClassesPerGene_HFDvsChow.csv",
                          Treat1Name = "HFD",
                          Treat2Name = "Chow"))

## ----gene_metabolitePairProfile,fig.align='center',fig.align='center',fig.width = 6, fig.height = 5----
create_edgeFBarplot(CompTreatTable = LCP_sharedGene_HFDvsChow,
                    edgeF = "Gm4553",
                    treat1 = "HFD",
                    treat2 = "Chow",
                    factorOrder = c("HFD","Chow"),
                    col1="#E76BF3",
                    col2 = "#F8766D",EdgeFeatureType = "Gene")

## ----global_metabolitePairProfile,fig.align='center',fig.width = 6, fig.height = 5----
(HFDvsChow_GlobalLipidProfile<-create_GlobalBarplot(CompTreatTable = LCP_sharedGene_HFDvsChow,
                                                     treat1 = "HFD",
                                                     treat2 = "Chow",
                                                     factorOrder = c("HFD","Chow"),
                                                     col1="#E76BF3",
                                                     col2 = "#F8766D", 
                                                     maxpairs = 1,
                                                     szggrepel = 6,
                                                     szaxisTxt = 15,
                                                     szaxisTitle = 15,
                                                     xlb = "Metabolite-pair classes"))

## ----stacked_metabolitePairProfile,fig.align='center',fig.show='hold',fig.width = 6, fig.height = 5----
create_stackedGlobalBarplot_perTreatment(CompTreatTable = LCP_sharedGene_HFDvsChow,
                                         treat = "HFD",
                                         max_pairsLegend = 1,
                                         xlb = "Metabolite-class-pairs",
                                         szTitle = 20,
                                         szggrepel = 6,
                                         szaxisTxt = 15, 
                                         szaxisTitle = 15)

create_stackedGlobalBarplot_perTreatment(CompTreatTable = LCP_sharedGene_HFDvsChow,
                                         treat = "Chow",
                                         max_pairsLegend = 1,
                                         ylim = 3,
                                         xlb = "Metabolite-pair classes",
                                         szTitle = 20,
                                         szggrepel = 4,
                                         szaxisTxt = 15,
                                         szaxisTitle = 15)

## ----stacked_metabolitePairProfile_grid,fig.align='center',fig.width = 7, fig.height = 5----
HFDvsChow_StackedLipidProfile<-getstackedGlobalBarplot_and_Grid(
                                 CompTreatTable = LCP_sharedGene_HFDvsChow,
                                 xlb = "Metabolite-pair classes",
                                 Treat1 = "HFD",
                                 Treat2 = "Chow",
                                 szTitle = 20,
                                 szggrepel = 6,
                                 szaxisTxt = 15, szaxisTitle = 15)
plot(HFDvsChow_StackedLipidProfile)

## ----metabolite_classes_per_gene,fig.align='center',fig.width = 7, fig.height = 5----
HFD_vs_Chow_LCP_Gene<-getVertexsPerEdgeFeature_and_Grid(LCP_sharedGene_HFDvsChow,"HFD","Chow",
                                                        Annotation=MetaboliteAnnotation,
                                                        ggrep=FALSE,
                                                        small = TRUE,
                                                        szTitle = 20,
                                                        szaxisTxt = 15, 
                                                        szaxisTitle = 15)
plot(HFD_vs_Chow_LCP_Gene)

