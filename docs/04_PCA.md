# PCA

__Principal component analysis (PCA)__ is a commonly used multivariate analysis in popgen, and represents one of the most useful techniques to visualise genetic diversity in a dataset. The methodology is not restricted to genetic data, but in general allows breaking down high-dimensional datasets to two or more dimensions for visualization in a two-dimensional space. Individual axes, known as __principal components (PCs)__ or __eigenvectors__, are extracted sequentially, with each PC being independent and encapsulating as much of the remaining variation as possible. Using PCA it is possible to estimate the proportion of the total variance in the dataset that has been summarized within these reduced dimensions.

PCA is used to display genetic distances between a __set of individuals__, rather than a __set of populations__. Each data point represents an individual person, plotted according to the eigenvalues of the first two PCs of the data calculated from the __covariance matrix__ of a genome-wide SNP dataset. Plotting values for two PCs for every individual allows the spatial separation of the different groups.

For the plots of genome-wide data, the first two PCs are normally shown, but any two PCs can be plotted to search for useful information. It should be remembered that, for genome-wide analysis, although the first two PCs will represent the largest amount of variation of all PCs, they may still represent a small amount of the total variation in the dataset.

![Population structure within Europe](https://www.researchgate.net/profile/John-Novembre/publication/264245693/figure/fig2/AS:295962284183553@1447574315970/A-statistical-summary-of-genetic-data-from-1387-Europeans-based-on-the-two-components_W640.jpg){width=100%}

Population structure within Europe revealed by genome-wide SNPs. Principal component analysis of data on 197,146 SNPs in 1387 Europeans. Small colored labels represent individuals and large colored points represent median PC1 and PC2 values for each country. The inset map provides a key to the labels. The PC axes are rotated to emphasize the similarity to the geographic map of Europe _(from [Novembre J et al. 2008](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2735096/))_.

Soon we will see that PC plots based on genome-wide SNP analyses show that world populations from the HGDP-CEPH collection cluster according to the continent, and as shown here European populations cluster approximately according to country. This indicates that our genomes contain detailed information about geographical ancestry, a finding with important implications for __persona genomics__ and __forensic genetics__.

****************************************************************************************************


### Creating PCA in plink

With the modifier `--pca` PLINK extracts the top 20 principal components of the variance-standardized relationship matrix. The results consist in a `.eigenvec` file with the coordinates for each individual in rows and eigenvectors in columns, and a `.eigenval` which explains how much variance there is in the data for that given vector.

``` bash
plink1.9 --bfile hgdp_newFIDs_modernPops_MindGeno --pca --out hgdp_newFIDs_modernPops_MindGeno_PCA
```

****************************************************************************************************

### Visualizing PCA in R

Now we can visualize PCA in R:

``` R
# Visualizing PCA in R

# Load necessary packages
library("tidyverse")

indINFO   <- read_table("HGDP_metainformation.txt", col_types = "cffffddffffffiii")

eigenvec <- read_table("hgdp_newFIDs_modernPops_MindGeno_PCA.eigenvec", 
                       col_names = c("FID","IID", paste("PC",1:20, sep = "")),
                       col_types = "fcdddddddddddddddddddd")

eigenval <- read_table("hgdp_newFIDs_modernPops_MindGeno_PCA.eigenval", col_names = "Eigenval") %>% 
  mutate(PC = 1:20)

PC_x <- 1
PC_y <- 3

# Plot eigenval
eigenval %>% 
  ggplot(aes(x=PC, y=Eigenval)) +
  geom_bar(stat="identity")

# Or we can calculate proportion for each eigenval from total eigenval and plot it
eigenval <- eigenval %>% 
  mutate(Prop = round(Eigenval/sum(Eigenval)*100,1)) 

eigenval %>% 
  ggplot(aes(x=PC, y=Prop)) +
  geom_bar(stat="identity")

# merge with indINFO table 
pca_inf <- left_join(eigenvec, indINFO, by="IID")

# Deffine shape for each population within each Region
shape_color <- unique(pca_inf[,c("FID","Region")]) %>% arrange(Region)
shape_color$Region_Population <- as.factor(paste(shape_color$Region, shape_color$FID, sep = " - "))
shape_color$Color <- c("#a6cee3","#1f78b4", "#b2df8a", "#33a02c","#fb9a99", "#e31a1c", "#fdbf6f")[shape_color$Region]
shape_color$Shape <- NA
for (i in unique(shape_color$Region)){
  shape_color[which(shape_color$Region == i),"Shape"] <- 1:nrow(shape_color[which(shape_color$Region == i),])
}


# Merge with PCA dataframe
pca_inf_shape_color <- left_join(pca_inf, shape_color[,-2], by="FID")

# PCA plot with colors according to Regions and shapes as Populations
pca_inf_shape_color %>% 
  ggplot(aes(x=unlist(.[paste0("PC", PC_x)]), y=unlist(.[paste0("PC", PC_y)]))) +
  xlab(paste0("PC", PC_x, " (", eigenval$Prop[PC_x], "%)")) + 
  ylab(paste0("PC", PC_y, " (", eigenval$Prop[PC_y], "%)")) +
  geom_point(aes(col=Region_Population, shape=Region_Population)) +
  scale_colour_manual(values = deframe(shape_color[,c("Region_Population","Color")])) +
  scale_shape_manual(values = deframe(shape_color[,c("Region_Population","Shape")])) +
  theme_bw()
```

Do you notice any outliers? Let's check other PCA's and figure out which individuals should be removed. 
As this is well established dataset, there is already indication in `~/popgen_intro/sample_all_snp.txt`. Search for "discover" in name of populations to find those that should be filter out.


!!! warning "Optional homework"
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