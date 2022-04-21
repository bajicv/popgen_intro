# Multivariate analyses

_How can we display a set of pairwise genetic distances in a comprehensible manner?_ 
If we have n populations, we require n dimensions to fully display their pairwise genetic distances on a graph. However, we often have more than four populations or molecules that we want to compare, yet we cannot conceive of, or represent, the four or more dimensions required to display these data. __Multivariate analyses__ allow us to reduce these multiple dimensions to the two or three dimensions we can draw and comprehend, while minimizing the inevitable loss of information. 

![How many dimensions are needed to display population relationships?](/pics/multivariate_analyses.jpeg){width=100%}

Figure from [Jobling et al. 2014](https://www.amazon.com/Human-Evolutionary-Genetics-Mark-Jobling/dp/081534148)

Distance matrices are shown relating genetic distances between increasing numbers of populations. As the number of populations increases, the number of dimensions needed to display these distances visually such that the lengths of the lines represent the genetic distances between them quantitatively also increases. A single line of a given length can relate two populations. Three populations can be related by a triangle, whose sides have lengths proportional to the genetic distance. Four populations can be represented by a pyramid. With five or more populations it is no longer possible to represent population relationships in three-dimensional space.


****************************************************************************************************

## MDS

__Multi-dimensional scaling (MDS)__ is one of these methods, where the loss of information is represented as a __stress statistic__, with the lower the stress statistic, the better the MDS fit to the data and therefore the less information lost. 

We can use plink to construct MDS for `.genome` file that we created
```{bash, echo=T, eval=F}
# First we can create 4 dimensions for MDS
plink1.9 --bfile hgdp_newFIDs_modernPops_MindGeno --read-genome hgdp_newFIDs_modernPops_MindGeno_genome.genome --cluster --mds-plot 4 --out mds_4D_hgdp


# And then we can calculate MDS with 2 dimensions
plink1.9 --bfile hgdp_newFIDs_modernPops_MindGeno --read-genome hgdp_newFIDs_modernPops_MindGeno_genome.genome --cluster --mds-plot 2 --out mds_2D_hgdp
```



```{r, MDS_plot, eval=T, echo=T, fig.width=15, fig.align='center'}

# MDS plot

# Loading necessary packages
library(data.table)
library(dplyr)
library("ggplot2")

# Specifying path to directory and loading files
setwd("~/Nowick_Lab/Master_course_Bioinformatics_for_Biologists/SoSe2020/popge_intro/")
indINFO <- fread("HGDP_metainformation.txt")
mds <- fread("mds_2D_hgdp.mds")
#mds <- fread("mds_4D_hgdp.mds") # try with 4 dimensions

# to vlookup to indINFO table 
mds_inf <- left_join(mds, indINFO, by="IID")

# Deffine shape for each population within each Region
symbol_col <- unique(mds_inf[,c("FID","Region")])
symbol_col$Symbol <- NA
symbol_col$Color <- NA

# assign shape/symbol
for (i in unique(symbol_col$Region)){
  symbol_col[which(symbol_col$Region == i),"Symbol"] <- 1:nrow(symbol_col[which(symbol_col$Region == i),])
}

# assign color
n = 0
for (i in unique(symbol_col$Region)){
  n = n + 1
  symbol_col[which(symbol_col$Region == i),"Color"] <- n
}

# Merge with mds dataframe
mds_inf_shape <- left_join(mds_inf, symbol_col[,c("FID","Symbol","Color")], by="FID")
mds_inf_shape$Region_Population <- paste(mds_inf_shape$Region, mds_inf_shape$Population, sep = " - ")

# Order dataset before plotting
ordered_mds_inf_shape <- mds_inf_shape[order(mds_inf_shape$Region_Population),]
# MDS plot with colors according to Regions and shapes as Populations
ggplot(mds_inf_shape, aes(x=C1, y=C2)) +
  geom_point(aes(col=Region_Population, shape=Region_Population)) +
  scale_colour_manual(values = c("#a6cee3","#1f78b4", "#b2df8a", "#33a02c","#fb9a99", "#e31a1c", "#fdbf6f")[as.factor(unique(ordered_mds_inf_shape[,c("Color","Region_Population")])$Color)]) +
  scale_shape_manual(values=unique(ordered_mds_inf_shape[,c("Symbol","Region_Population")])$Symbol) +
  theme_bw()

```



****************************************************************************************************

## PCA

It is also possible to use multivariate analysis to generate two- or three-dimensional graphical representations of distances between populations using the _raw data of allele frequencies_, rather than _genetic distances_. __Principal component analysis (PCA)__ is a commonly used example of this approach, and represents one of the most useful techniques to visualise genetic diversity in a dataset. The methodology is not restricted to genetic data, but in general allows breaking down high-dimensional datasets to two or more dimensions for visualization in a two-dimensional space. Individual axes, known as __principal components (PCs)__ or __eigenvectors__, are extracted sequentially, with each PC being independent and encapsulating as much of the remaining variation as possible. Using PCA it is possible to estimate the proportion of the total variance in the dataset that has been summarized within these reduced dimensions. PCA can reveal global relationships from a set of _pairwise distances between populations_.

PCA is also used to display genetic distances between a __set of individuals__, rather than a __set of populations__. Each data point represents an individual person, plotted according to the eigenvalues of the first two PCs of the data calculated from the __covariance matrix__ of a genome-wide SNP dataset. Plotting values for two PCs for every individual allows the spatial separation of the different groups.

For the plots of genome-wide data, the first two PCs are normally shown, but any two PCs can be plotted to search for useful information. It should be remembered that, for genome-wide analysis, although the first two PCs will represent the largest amount of variation of all PCs, they may still represent a small amount of the total variation in the dataset. Principal components can also be rotated, a statistical procedure that is sometimes used to emphasize the relationship between the genetic structure and geography.

![Population structure within Europe](/pics/PCA.jpeg){width=100%}

Population structure within Europe revealed by genome-wide SNPs. Principal component analysis of data on 197,146 SNPs in 1387 Europeans. Small colored labels represent individuals and large colored points represent median PC1 and PC2 values for each country. The inset map provides a key to the labels. The PC axes are rotated to emphasize the similarity to the geographic map of Europe (from [Novembre J et al. 2008](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2735096/)).

Soon we will see that PC plots based on genome-wide SNP analyses show that world populations from the HGDP-CEPH collection cluster according to the continent, and as shown here European populations cluster approximately according to country. This indicates that our genomes contain detailed information about geographical ancestry, a finding with important implications for __persona genomics__ and __forensic genetics__.

****************************************************************************************************


### PCA in plink

With the modifier `--pca` PLINK extracts the top 20 principal components of the variance-standardized relationship matrix. The results consist in a `.eigenvec` file with the coordinates for each individual in rows and eigenvectors in columns, and a `.eigenval` which explains how much variance there is in the data for that given vector.

```{bash plink_PCA, echo=T, eval=F}
plink1.9 --bfile hgdp_newFIDs_modernPops_MindGeno --pca --out hgdp_newFIDs_modernPops_MindGeno_PCA
```


Now we can visualize PCA in R:
```{r PCA_plot, echo=T, eval=T, fig.width=10, fig.align='center'}

# PCA

# Loading necessary packages
library(data.table)
library(dplyr)
library("ggplot2")

# Specifying path to directory and loading files
setwd("~/Nowick_Lab/Master_course_Bioinformatics_for_Biologists/SoSe2020/popge_intro/")

indINFO <- fread("HGDP_metainformation.txt")

eigenvec <- fread("hgdp_newFIDs_modernPops_MindGeno_PCA.eigenvec")
names(eigenvec) <- c("FID","IID", paste("PCA",1:20, sep = "_"))

eigenval <- fread("hgdp_newFIDs_modernPops_MindGeno_PCA.eigenval")
eigneval <- data.frame("PCA" = 1:20, "Eigenval" = eigenval$V1)

# Plot eigenval
ggplot(data=eigneval, aes(x=PCA, y=Eigenval_prop)) +
  geom_bar(stat="identity")

# Or we can calculate proportion for each eigenval from total eigenval and plot it
eigneval$Eigenval_prop <- round(eigneval$Eigenval / sum(eigneval$Eigenval)*100,2)

ggplot(data=eigneval, aes(x=PCA, y=Eigenval_prop)) +
  geom_bar(stat="identity")


# Preparing to plot PCA
# to vlookup to indINFO table 
pca_inf <- left_join(eigenvec, indINFO, by="IID")

# Deffine shape for each population within each Region
symbol_col <- unique(pca_inf[,c("FID","Region")])
symbol_col$Symbol <- NA
symbol_col$Color <- NA

# assign shape/symbol
for (i in unique(symbol_col$Region)){
  symbol_col[which(symbol_col$Region == i),"Symbol"] <- 1:nrow(symbol_col[which(symbol_col$Region == i),])
}

# assign color
n = 0
for (i in unique(symbol_col$Region)){
  n = n + 1
  symbol_col[which(symbol_col$Region == i),"Color"] <- n
}

# Merge with mds dataframe
pca_inf_shape <- left_join(pca_inf, symbol_col[,c("FID","Symbol","Color")], by="FID")
pca_inf_shape$Region_Population <- paste(pca_inf_shape$Region, pca_inf_shape$Population, sep = " - ")

# Order dataset before plotting
ordered_pca_inf_shape <- pca_inf_shape[order(pca_inf_shape$Region_Population),]
# PCA plot with colors according to Regions and shapes as Populations
ggplot(pca_inf_shape, aes(x=PCA_1, y=PCA_2)) +
  geom_point(aes(col=Region_Population, shape=Region_Population)) +
  scale_colour_manual(values = c("#a6cee3","#1f78b4", "#b2df8a", "#33a02c","#fb9a99", "#e31a1c", "#fdbf6f")[as.factor(unique(ordered_pca_inf_shape[,c("Color","Region_Population")])$Color)]) +
  scale_shape_manual(values=unique(ordered_pca_inf_shape[,c("Symbol","Region_Population")])$Symbol) +
  theme_bw()
```

Do you notice any outliers? Let's check other PCA's and figure out which individuals should be removed. 
As this is well established dataset, there is already indication in `~/popgen_intro/sample_all_snp.txt`


!!! warning "Homework"
    Make a new PCA after removing the outliers.

    1. Make file with individuals to be removed;
      
    2. Use this file to subset your original file;
      
    3. Run PCA on new file that doesn't contain outlier individuals; 
      
    4. Visualize PCA in R.
    
    _Do you see any differences?_

    *****************

    ??? warning "Click for Answer"
        <br />1. Make file with individuals to be removed;
      
        ``` bash
        cd ~/popgen_intro/plink_exercise
        grep discover ~/popgen_intro/sample_all_snp.txt > ind_exclude_tmp.txt
        join -j 1 -o 2.2,1.1 <(sort -k1 ind_exclude.txt) <(sort -k1 ../HGDP_metainformation.txt) > ind_exclude.txt
        ```
            
        <br />2. Use this file to subset your original file;
            
        ``` bash
        plink1.9 --bfile hgdp_newFIDs_modernPops_MindGeno --remove ind_exclude.txt --make-bed --out hgdp_outliers_removed
        ```
            
        <br />3. Run PCA on new file that doesn't contain outlier individuals; 
            
        ``` bash
        plink1.9 --bfile hgdp_outliers_removed --pca --out hgdp_outliers_removed_PCA
        ```
            
        <br />4. Download necessary ouput files to your local machine and visualize PCA in R using the script abowe. Don't forget to exchange input files for `eigenval` and `eigenvec`.