### IBD and PI_HAT

https://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1002287

<center>
![Consanguinity - degrees of relationships](/pics/Consanguinity.jpeg){width=100%}
_Consanguinity - degrees of relationships_<br />
_Consanguinity ("blood relation", from Latin consanguinitas) is the characteristic of having a kinship with another person (being descended from a common ancestor)._
</center> 

The HGDP data set is already filtered for relatives according to [Rosenberg et al. 2006](https://www.ncbi.nlm.nih.gov/pubmed/17044859). There are different ways to determine related individuals. Here we will mention plink's widely used __IBD__ estimate __PI_HAT__ ($\hat{π}$) which can be used to estimate the level of relatedness between a pair of individuals. PI_HAT is a measure of overall IBD alleles. It measures the fraction of the genome shared (IBD) between pairs of the individuals (P(IBD=2)+0.5*P(IBD=1)). 

Since our data is already filtered for relatives using more sophisticated methodologies (see [Rosenberg et al. 2006](https://www.ncbi.nlm.nih.gov/pubmed/17044859)) we will only use PI_HAT to search for individuals that still have very high PI_HAT. This will illustrate to us how isolation, consanguinity, bottleneck, drift... can violate simple principles on which PI_HAT relies on (e.g. populations from Oceania and Americas).

``` bash
plink1.9 --bfile hgdp_newFIDs_modernPops_MindGeno --genome --out hgdp_newFIDs_modernPops_MindGeno_genome
```

Check plink webpage for more information about the [output of the `--genome` option](https://www.cog-genomics.org/plink/1.9/formats#genome).

![Degrees of relationships](https://upload.wikimedia.org/wikipedia/commons/thumb/5/51/Pedigree%2C_recombination_and_resulting_IBD_segments%2C_schematic_representation_modified.png/1280px-Pedigree%2C_recombination_and_resulting_IBD_segments%2C_schematic_representation_modified.png){width=100%}

Columns Z0, Z1, and Z2 indicate the probabilities of having IBD of 0, 1, or 2 over the loci, which gives us a way of discriminating between relationship types. 
Ideal parent-offspring has (Z0, Z1, Z2) = (0, 1, 0), i.e. all loci have one allele identical by descent; 
Ideal full sibling = (1/4, 1/2, 1/4), i.e. 25% of loci have 0 alleles IBD, 50% have 1 allele IBD, 25% have 2 alleles IBD; etc. So, we can use Z0, Z1, and Z2 to distinguishing between relationship types. 

![IBD 0 1 2](https://upload.wikimedia.org/wikipedia/commons/0/0a/IBD_0_1_2.png)

IBD=0, IBD=1, and IBD=2 are the genome proportions shared on 0, 1, and 2 chromosomes, respectively, between two individuals. Many relationships share the same expected mean IBD proportions; however, for full-sibling, second-degree, and third-degree relationships, a variance around the expected mean is due to the random nature of recombination events. Genotyping and other technical errors can contribute to this variance.

The table below shows the expected mean IBD proportions for the outbred familial relationship categories

<center>

| Familial Relationship |	IBD0 |	IBD1 | IBD2 |
|-----------------------|------|-------|------|
|Parental               |	0    |	1    |	0   |
|Full-sibling           |	0.25 |	0.5  |	0.25|
|Half-sibling, avuncular, and grandparental |	0.5 |	0.5 |	0 |
|First-cousin, great-grandparental, great-avuncular, and half-avuncular |	0.75 |	0.25 |	0 |
|Distantly related      |	varies |	varies |	0 |
|Unrelated (includes relationships beyond the third degree) |	1 |	0 |	0 |

</center>

In theory, PI_HAT values that are around: 
<br />0.125 are __third degree__, 
<br />0.25 are __second degree__, and 
<br />0.5 are __first degree__. 
<br />


PI_HAT assumes homogeneous population with known population allele frequencies, however if the population is not homogeneous, the values are biased compared to the expectations under the known relationship categories. 

In other words, plink's PI_HAT is based on the assumption that the study sample is drawn from a single, homogeneous, randomly mating population. This assumption is violated if pedigree founders are drawn from multiple populations or include admixed individuals. In the presence of population structure, the method of moments estimator has an inflated variance and can be biased because it relies on sample-based allele frequency estimates. For plink's estimator which truncates genome-wide sharing estimates at zero and one to generate biologically interpretable results, the bias is most often towards over-estimation of relatedness between ancestrally similar individuals. 

Thus the PI_HAT/Z0/Z1/Z2 values produced by `plink --genome` should not be taken as solid evidence to justify relationship inference. Instead it is better to use [KING](http://people.virginia.edu/~wc9c/KING/) or [PC-Relate](https://www.rdocumentation.org/packages/GENESIS/versions/2.2.2/topics/pcrelate).


*******************************************************

### Heterozygosity - inbreeding - consanguinity
https://link.springer.com/article/10.1007/s00439-019-02045-1

A __method-of-moments estimator__ (__F__) developed by Ritland (1996) can be used to estimate the distribution of __inbreeding__ among and __relatedness__ between individuals of a natural population.

Plink's `--het` option computes observed and expected autosomal homozygous genotype counts for each sample, and reports method-of-moments F coefficient estimates (i.e. ([observed hom. count] - [expected count]) / ([total observations] - [expected count])) to `plink.het`.

```
plink1.9 --bfile hgdp_newFIDs_modernPops_MindGeno --het --out hgdp_newFIDs_modernPops_MindGeno_het
```

We can visualize the difference in heterozygosity within populations and we can visualize it on the map like this:
```{r, echo=T, eval=T, fig.align='center'}
# Loading necessary packages
library(data.table)
library(dplyr)
library("ggplot2")

# Specifying path to directory and loading files
setwd("~/Nowick_Lab/Master_course_Bioinformatics_for_Biologists/SoSe2020/popge_intro/")
indINFO <- fread("HGDP_metainformation.txt")
het <- fread("hgdp_newFIDs_modernPops_MindGeno_het.het")


#to vlookup to indINFO table 
het <- left_join(het, indINFO, by="IID")

# Preparing to plot Missingness boxplots per population
# First ordering the populations based on region and color so that they appear like that on the boxplot
ordered_het <- het[order(het$Region, het$Country, het$FID),]
het$FID <- factor(het$FID, levels = unique(ordered_het$FID))

# ploting boxplots per population colored by region

ggplot(het, aes(x=FID, y=F, fill=Region)) + 
  geom_boxplot() + 
  theme(axis.text.x = element_text(angle = 90))


# Calculate mean F per population
pop_mean_F =  aggregate(F ~ FID, het, mean )

# Extract metainformation for each population and merge it with pop_mean_F 
uniqMetaPop <- unique(het[,c("FID", "Region", "Latitude", "Longitude")])
pop_mean_F_extra = left_join(pop_mean_F, uniqMetaPop, by="FID")


library("maps")
world_map <- map_data(map="world")

ggplot() + 
  theme(panel.background = element_rect(fill = "white")) +
  geom_map(data=world_map, map=world_map, 
           aes(map_id=region), fill="black", colour="gray20", size=0.15) +
  coord_quickmap(ylim=c(-35,70), xlim=c(-115,150)) +
  geom_point(data=pop_mean_F_extra, aes(x=Longitude, y=Latitude, color=F), alpha =0.5, size=5 )+
  scale_color_gradientn(colours = rainbow(10)) +
  ggtitle("Average method-of-moments [F] coefficient estimates for HGDP populations")
```

In which populations F has lowest and highest values? How would you interpret these results?
