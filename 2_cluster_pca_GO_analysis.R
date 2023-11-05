
#This script is for clustering RNASeq data,identifying differentially expressed genes and preforming GO enrichment analysis 
#Please run the 1_import_filter_norm.R script before this one to ensure all necessary variables are defined 

#Installing packages----
# install.packages('DT')
# install.packages('plotly')
# install.packages('gt')

# Loading the required libraries -----
library(tidyverse) # provides access to Hadley Wickhma's collection of R packages for data science 
library(DT) # for making interactive tables
library(plotly) # for making interactive plots
library(gt) # A layered 'grammar of tables' - like ggplot, but for tables

log2.cpm.filtered.norm.df

# ---- HIERARCHICAL CLUSTERING ----
distance_euc <- dist(t(log2.cpm.filtered.norm), method = "euclidean")  # compute distances based on the normalized and filtered data
clusters_euc <- hclust(distance_euc, method = "complete") # cluster the data based on distance
c_euc <- plot(clusters_euc, labels=sampleLabels) # plot the dendrogram        

# ---- PRINCIPAL COMPONENT ANALYSIS(PCA)----
# Identify variables of interest in study design file ----
targets
group <- targets$group
group <- factor(group)

pca.res <- prcomp(t(log2.cpm.filtered.norm), scale.=F, retx=T) #preform PCA 

#look at the PCA result (pca.res) that you just created
ls(pca.res)
summary(pca.res) # Prints variance summary for all principal components.
pca.res$rotation #$rotation shows you how much each gene influenced each PC (called 'scores')
pca.res$x # 'x' shows you how much each sample influenced each PC (called 'loadings')
#note that these have a magnitude and a direction (this is the basis for making a PCA plot)
screeplot(pca.res) # A screeplot is a standard way to view eigenvalues for each PCA         

pc.var<-pca.res$sdev^2 # sdev^2 captures these eigenvalues from the PCA result
pc.per<-round(pc.var/sum(pc.var)*100, 1) # we can then use these eigenvalues to calculate the percentage variance explained by each PC
pc.per

# Visualize your PCA result ------------------
#lets first plot any two PCs against each other
#We know how much each sample contributes to each PC (loadings), so let's plot
pca.res.df <- as_tibble(pca.res$x)
ggplot(pca.res.df) +
  aes(x=PC1, y=PC2, label=sampleLabels) + # color = group
  geom_point(size=4) +
  geom_label() +
  #stat_ellipse() +
  xlab(paste0("PC1 (",pc.per[1],"%",")")) + 
  ylab(paste0("PC2 (",pc.per[2],"%",")")) +
  labs(title="PCA plot",
       caption=paste0("produced on ", Sys.time())) +
  #coord_fixed() +
  theme_bw()                      

# Create a PCA 'small multiples' chart ----
pca.res.df <- pca.res$x[,1:4] %>% 
  as_tibble() %>%
  add_column(sample = sampleLabels,
             group = group)
  
pca.pivot <- pivot_longer(pca.res.df, # dataframe to be pivoted
                          cols = PC1:PC4, # column names to be stored as a SINGLE variable
                          names_to = "PC", # name of that new variable (column)
                          values_to = "loadings") # name of new variable (column) storing all the values (data)

ggplot(pca.pivot) +
  aes(x=sample, y=loadings, fill=group) +
  geom_bar(stat="identity") +
  facet_wrap(~PC) +
  labs(title="PCA 'small multiples' plot",
       caption=paste0("produced on ", Sys.time())) +
  theme_bw() +
  coord_flip()                

#---- GENE EXPRESSION AVERAGES
mydata.df <- log2.cpm.filtered.norm.df %>% 
  mutate(control.AVG = (C01 + C02 + C03)/3, # calculate averages for control 
         salt.AVG = (S01 + S02 + S03)/3, # calculate averages for treatment
         LogFC = (salt.AVG - control.AVG)) %>% # create columns comparing each of the averages above
  mutate_if(is.numeric, round, 2)

#now look at this modified data table
mydata.df

mydata.sort <- mydata.df %>%
  dplyr::arrange(desc(LogFC)) %>% 
  dplyr::select(geneID, LogFC)

mydata.filter <- mydata.df %>%
  dplyr::filter(geneID=="Solyc06g043170.3.1" | geneID=="Solyc08g076580.2.1" | geneID=="Solyc04g011650.3.1" | geneID=="Solyc02g087935.1.1" | geneID=="Solyc02g085370.3.1"
                | geneID=="Solyc05g044580.3.1" | geneID=="Solyc02g065220.3.1" | geneID=="Solyc07g025390.3.1" | geneID=="Solyc08g016080.3.1" ) %>%
  dplyr::select(geneID, control.AVG, salt.AVG, LogFC) %>%
  dplyr::arrange(desc(LogFC))

# you can also filter based on any regular expression
mydata.grep <- mydata.df %>%
  dplyr::filter(grepl('CXCL|IFI', geneID)) %>%
  dplyr::select(geneID, control.AVG, salt.AVG, LogFC) %>%
  dplyr::arrange(desc(geneID))


# Produce publication-quality tables using the gt package ----
gt(mydata.filter)                                     

# Make an interactive table using the DT package ----
datatable(mydata.df[,c(1,10:8)], #view as a table 
          extensions = c('KeyTable', "FixedHeader"), 
          filter = 'top',
          options = list(keys = TRUE, 
                         searchHighlight = TRUE, 
                         pageLength = 10, 
                         lengthMenu = c("10", "25", "50", "100")))   #10.000

# Make an interactive scatter plot with plotly -----
# begin by storing your ggplot object
myplot <- ggplot(mydata.df) + 
  aes(x=control.AVG, y=salt.AVG) +
  geom_point(shape=16, size=1) +
  ggtitle("salt vs. control") +
  theme_bw()

#now use the ggplotly function from the plotly package to convert this ggplot object into an interactive plot
ggplotly(myplot)                                 

#let's customize this graphic by adding a more informative mouseover tooltip
myplot <- ggplot(mydata.df) +
  aes(x=control.AVG, y=salt.AVG, 
      text = paste("Symbol:", geneID)) +
  geom_point(shape=16, size=1) +
  ggtitle("salt vs. control") +
  theme_bw()

ggplotly(myplot)  

# ---- IDENTIFICATION of DIFFERENTIALLY EXPRESSED GENES (DEGs)----

# Install new packages -----
# install.packages("reshape2")
# install.packages('heatmaply')

# Load packages -----
library(tidyverse) 
library(limma) # package for differential gene expression using linear modeling
library(edgeR)
library(gt)
library(DT)
library(plotly)
library(reshape2)
library(heatmaply)

# Set up your design matrix ----
group <- factor(targets$group)
design <- model.matrix(~0 + group)
colnames(design) <- levels(group)

v.DEGList.filtered.norm <- voom(myDGEList.filtered.norm, design, plot = TRUE) ## Use VOOM function from Limma package to model the mean-variance relationship

fit <- lmFit(v.DEGList.filtered.norm, design) # fit a linear model to the data 

contrast.matrix <- makeContrasts(infection = salt - control,  # create a contrast matrix
                                 levels=design)

fits <- contrasts.fit(fit, contrast.matrix) #extract the linear model fit 

ebFit <- eBayes(fits) #get bayesian stats for your linear model fit

#---VIEW DEGs---
myTopHits <- topTable(ebFit, adjust ="BH", coef=1, number=40000, sort.by="logFC") #number is set to capture all DEGs

# convert to a tibble
myTopHits.df <- myTopHits %>%
  as_tibble(rownames = "geneID")

gt(myTopHits.df)          

#----Create a volcano plot of the DEGs----
vplot = ggplot(myTopHits.df) + # vplot <- ggplot(myTopHits.df) +
  aes(y=-log10(adj.P.Val), x=logFC, text = paste("Symbol:", geneID)) +
  geom_point(size=2) +
  geom_hline(yintercept = -log10(0.01), linetype="longdash", colour="grey", linewidth=1) +
  geom_vline(xintercept = 1, linetype="longdash", colour="#BE684D", linewidth=1) +
  geom_vline(xintercept = -1, linetype="longdash", colour="#2C467A", linewidth=1) +
  annotate("rect", xmin = 1, xmax = 12, ymin = -log10(0.01), ymax = 7.5, alpha=.2, fill="#BE684D") +
  annotate("rect", xmin = -1, xmax = -12, ymin = -log10(0.01), ymax = 7.5, alpha=.2, fill="#2C467A") +
  labs(title="Volcano plot",
       subtitle = "Solanum Lycopersicum",
       caption=paste0("produced on ", Sys.time())) +
  theme_bw()

# Now make the volcano plot above interactive with plotly
ggplotly(vplot)       


# decideTests to pull out the DEGs and make Venn Diagram ----
results <- decideTests(ebFit, method="global", adjust.method="BH", p.value=0.05, lfc=4) # output -1 or 1: t-statistic for gene is classified as significant(0 = not significant)

# take a look at what the results of decideTests looks like
head(results)
summary(results)
vennDiagram(results, include="both")

# retrieve expression data for your DEGs ----
head(v.DEGList.filtered.norm$E)
colnames(v.DEGList.filtered.norm$E) <- sampleLabels


diffGenes <- v.DEGList.filtered.norm$E[results[,1] !=0,] # E:numeric matrix of normalised expression values on the log2 scale 
head(diffGenes)
dim(diffGenes)
#convert your DEGs to a dataframe using as_tibble
diffGenes.df <- as_tibble(diffGenes, rownames = "geneID")


# create interactive tables to display your DEGs ----
datatable(diffGenes.df,
          extensions = c('KeyTable', "FixedHeader"),
          caption = 'Table 1: DEGs in Solanum lycopersicum',
          options = list(keys = TRUE, searchHighlight = TRUE, pageLength = 10, lengthMenu = c("10", "25", "50", "100"))) %>%
  formatRound(columns=c(1:6), digits=2)


#write your DEGs to a file
write_tsv(diffGenes.df,"DiffGenes.txt")

#---- Create a heatmap of differentially expressed genes ----

heatmaply(diffGenes.df[2:7], 
          #dendrogram = "row",
          xlab = "Samples", ylab = "DEGs", 
          main = "DEGs in Solanum Lycopersicum",
          scale = "column",
          margins = c(60,100,40,20),
          grid_color = "white",
          grid_width = 0.0000001,
          titleX = T,
          titleY = T,
          hide_colorbar = FALSE,
          branches_lwd = 0.1,
          label_names = c("Gene", "Sample:", "Value"),
          fontsize_row = 5, fontsize_col = 5,
          labCol = colnames(diffGenes.df)[2:7],
          labRow = diffGenes.df$geneID,
          heatmap_layers = theme(axis.line=element_blank())
)



#---- GENE ONTOLOGY (GO) ENRICHMENT USING gProfiler2 ---

# install.packages("gprofiler2")
# Load packages ----
library(gprofiler2) #tools for accessing the GO enrichment results using g:Profiler web resources

# use topTable result to pick the top genes for carrying out a Gene Ontology (GO) enrichment analysis
myTopHits <- topTable(ebFit, adjust ="BH", coef=1, number=50, sort.by="logFC") # number:50 pick the top 50 genes by logFC value for carrying out GO enrichment analysis 
# use the 'gost' function from the gprofiler2 package to run GO enrichment analysis
gost.res <- gost(rownames(myTopHits), organism = "slycopersicum", correction_method = "fdr")
gost.res
# produce an interactive manhattan plot of enriched GO terms
mygostplot = gostplot(gost.res, interactive = F, capped = F) # set interactive = F to publish the plot  
# produce a publication quality static manhattan plot with specific GO terms highlighted
publish_gostplot(
  mygostplot, #your static gostplot from above
  highlight_terms = c("GO:0016491","GO:0001101","GO:0009414","GO:0044036","GO:1901700","GO:1902074"), # higHlight lowest  p values ones 
  filename = NULL,
  width = NA,
  height = NA)




