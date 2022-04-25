# Plink

We will be working with `plink`, a free, open-source whole-genome association analysis toolset, designed to perform a range of basic, large-scale analyses in a computationally efficient manner. Even though the focus of `plink` is  on analysis of genotype/phenotype data, it is widely used in popgen as it has many features for data manipulation, it offers basic statistics, and many popgen tools assume input files to be in plink format (e.g. `fastStructure`, `ADMIXTURE`, etc.).

[Plink](https://www.cog-genomics.org/plink/1.9/) is a tool created for genome-wide association studies ([GWAS](https://www.genome.gov/genetics-glossary/Genome-Wide-Association-Studies#:~:text=A%20genome%2Dwide%20association%20study,the%20presence%20of%20a%20disease.)) and research in population genetics. `Plink` parses each command line as a collection of flags (each of which starts with two dashes), plus parameters (which immediately follow a flag). Because Plink was developed for GWAS medical studies, many statistics and pieces of information from plink files will be not used in our analysis, such as pedigree or phenotype.

*********************************************************************

## Plink file formats

`plink` can either read 
<br />1. __text-format files__ (`.ped` + `.map`) or 
<br />2. __binary files__ (`.bed` + `.bim` + `.fam`). 

Because reading large text files can be time-consuming, it is recommended to use binary files. 

![plink data formats (Marees et al. 2017)](https://www.researchgate.net/publication/323424714/figure/fig3/AS:667766705098757@1536219397189/Overview-of-various-commonly-used-PLINK-files-SNP-single-nucleotide-polymorphism_W640.jpg){width=80%}

*********************************************************************

### Text-format plink files

Text `plink` data consist of two files (which should have matching names and they should be stored together): 
<br />1. `.ped` contains information on the __individuals and their genotypes__; 
<br />2. `.map` contains information on the __genetic markers__.

*********************************************************************

### Binary plink files

Binary `plink` data consist of three files (one binary and two text files which should have matching names and which should be stored together):
<br />1. `.bed`  contains individual identifiers (IDs) and genotypes, 
<br />2. `.fam` contain information on the individuals, 
<br />3. `.bim` contains information on the genetic markers. 

Analysis using __covariates__ often requires the fourth file, containing the values of these covariates for each individual.

*********************************************************************

## Basic statistics with plink

We can generate some simple summary statistics using `plink` e.g. rates of missing data in the file, diversity within and between the samples, etc.

Before we start our exercises let's make binary files out of our text-format plink files and output them in `~/popgen_intro/plink_exercise/`:

``` bash
cd ~/popgen_intro
mkdir plink_exercise
plink1.9 --file all_snp --make-bed --out ~/popgen_intro/plink_exercise/hgdp
```

Now go to `plink_exercise` and check which files are there:

``` bash
$ cd ~/popgen_intro/plink_exercise/
$ ls

hgdp.bed  hgdp.bim  hgdp.fam  hgdp.log 
```

As expected we got 3 files that go together in binary file format in `plink` (i.e. `bed` + `bim` + `fam`), but also `log`. Let's take a look at it:

``` bash
$ cat hgdp.log

PLINK v1.90b6.6 64-bit (12 Oct 2018)
Options in effect:
  --file all_snp
  --make-bed
  --out /home/bajiv90/popgen_intro/plink_exercise/hgdp

Hostname: evop-login
Working directory: /home/bajiv90/popgen_intro
Start time: Mon Apr 25 18:19:33 2022

Random number seed: 1650903573
16042 MB RAM detected; reserving 8021 MB for main workspace.
Scanning .ped file... done.
Performing single-pass .bed write (62679 variants, 942 people).
--file: /home/bajiv90/popgen_intro/plink_exercise/hgdp-temporary.bed +
/home/bajiv90/popgen_intro/plink_exercise/hgdp-temporary.bim +
/home/bajiv90/popgen_intro/plink_exercise/hgdp-temporary.fam written.
62679 variants loaded from .bim file.
942 people (621 males, 321 females) loaded from .fam.
Using 1 thread (no multithreaded calculations invoked).
Before main variant filters, 942 founders and 0 nonfounders present.
Calculating allele frequencies... done.
Warning: Nonmissing nonmale Y chromosome genotype(s) present; many commands
treat these as missing.
Total genotyping rate is 0.995297.
62679 variants and 942 people pass filters and QC.
Note: No phenotypes present.
--make-bed to /home/bajiv90/popgen_intro/plink_exercise/hgdp.bed +
/home/bajiv90/popgen_intro/plink_exercise/hgdp.bim +
/home/bajiv90/popgen_intro/plink_exercise/hgdp.fam ... done.

End time: Mon Apr 25 18:19:35 2022
```

`plink` by default creats `.log` file with information on how linked `bed`, `bim`, `fam` or `map` and `ped` files were created. This is very useful when trying to understand your own work after some time, or understanding the steps of your collaborators.

We can include population names as family ID (FID) so we can have an idea to which population each individual belongs. We can get this information from the metainformation file that we have in `~/popgen_intro/HGDP_metainformation.txt`. To do so we should first make a file with old and new IID and FID, and then make new plink files with updated IDs (`--update-ids`):

``` bash
cd ~/popgen_intro/plink_exercise/

join -j 1 -o 1.1,1.2,2.2,2.1 <(sort -k1 hgdp.fam) <(sort -k1 ../HGDP_metainformation.txt) > hgdp_updateIDs.txt

plink1.9 --bfile hgdp --update-ids hgdp_updateIDs.txt --make-bed --out hgdp_newFIDs
```

!!! tip "__Tip: Process Substitution__"
    __Piping__ (`|`) the stdout of a command into the stdin of another is a powerful technique. But, what if we need to pipe the stdout of multiple commands? This is where __process substitution__ (`<(command)` or `>(command)`) comes in.

    Process substitution feeds the output of a process (or processes) into the stdin of another process. In other words process substitution treats the output (or input) of commands as files. 
    <br /> `>(command list)`  will treat it as an input
    <br /> `<(command list)`  will treat it as an output

    The syntax for process substitution is: `<(list)` or `>(list)`
    <br /> e.g. `diff -u <( ls dir1 | sort)  <( ls dir2 | sort )`



***************************************************************************

### Missing data

We can use `--missing` to produces __sample-based__ (`.imiss`) and __variant-based__ (`.lmiss`) files with missing data reports.

``` bash
cd ~/popgen_intro/plink_exercise/

plink1.9 --bfile hgdp_newFIDs --missing --out hgdp_newFIDs_miss
```

Now, we can visualize the data in __R__. Let's download newly created files together with the metainformation file to our local machine, and then run the following R commands to plot them:

``` R
# load necessary packages
library("tidyverse")

# load files
indINFO   <- read_table("HGDP_metainformation.txt", col_types = "cffffddffffffiii")
imiss     <- read_table("hgdp_newFIDs_miss.imiss", col_types = "fcfiid")
lmiss     <- read_table("hgdp_newFIDs_miss.lmiss", col_types = "fciid")

# add metainformation to imiss
imiss <- left_join(imiss, indINFO, by="IID")

# order levels for FID to be ordered by Region, Country and FID
imiss$FID <- factor(imiss$FID, levels = unique(imiss[order(imiss$Region, imiss$Country, imiss$FID),]$FID))

# ploting boxplots per population colored by region
imiss %>% 
  ggplot(aes(x=FID, y=F_MISS, fill=Region)) + 
  geom_boxplot() + 
  theme(axis.text.x = element_text(angle = 90))

# Let's keep only modern human populations
# and plot it again so we can see better which of modern populations has a lot of individuals that are missing a lot of data
imiss %>% 
  filter(!FID %in% c("Href","Clint","Denisova", "Gorilla","Macaque","Marmoset","Orang","Vindija")) %>% 
  ggplot(aes(x=FID, y=F_MISS, fill=Region)) + 
  geom_boxplot() + 
  theme(axis.text.x = element_text(angle = 90))

# Plot histograms of Missingness per Chromosome
lmiss %>% 
  ggplot(aes(x=CHR, y=F_MISS, fill=CHR)) + 
  geom_boxplot() +
  theme(legend.position = "none")

# Now let's see only those that have F_MISS < 0.1
lmiss %>% 
  filter(F_MISS < 0.1) %>% 
  ggplot(aes(x=CHR, y=F_MISS, fill=CHR)) + 
  geom_boxplot() +
  theme(legend.position = "none")
```


Now, let's keep only individuals from modern populations, and remove SNPs and individuals with more than 10% missingness. To do so, we will need to specify which individuals we want to keep (with option `--keep`) or remove (with option `--remove`). One way to make list of individuals which we want to keep is to find all individuals that contain "HGDP" in their name.

``` bash
grep HGDP hgdp_newFIDs.fam > ind_keep.txt
plink1.9 --bfile hgdp_newFIDs --make-bed --mind 0.1 --geno 0.1 --keep ind_keep.txt --out hgdp_newFIDs_modernPops_MindGeno
```
!!! question
    How many individuals and how many SNPs were removed?


***************************************************************************


### IBD and PI_HAT

![Consanguinity - degrees of relationships](/pics/Consanguinity.jpeg){width=100%}

The HGDP data set is already filtered for relatives according to [Rosenberg et al. 2006](https://www.ncbi.nlm.nih.gov/pubmed/17044859). There are different ways to determine related individuals. Here we will mention plink's widely used __IBD__ estimate __PI_HAT__ ($\hat{\pi}$) which can be used to estimate the level of relatedness between a pair of individuals. PI_HAT is a measure of overall IBD alleles. It measures the fraction of the genome shared (IBD) between pairs of the individuals (P(IBD=2)+0.5*P(IBD=1)). 

Since our data is already filtered for relatives using more sophisticated methodologies (see Rosenberg et al. 2006) we will only use PI_HAT to search for individuals that still have very high PI_HAT. This will illustrate to us how isolation, consanguinity, bottleneck, drift can violate simple principles on which PI_HAT relies on (e.g. populations from Oceania and Americas).

``` bash
plink1.9 --bfile hgdp_newFIDs_modernPops_MindGeno --genome --out hgdp_newFIDs_modernPops_MindGeno_genome
```

Check plink webpage for more information about the [output of the `--genome` option](https://www.cog-genomics.org/plink/1.9/formats#genome).

![Degrees of relationships](https://upload.wikimedia.org/wikipedia/commons/thumb/5/51/Pedigree%2C_recombination_and_resulting_IBD_segments%2C_schematic_representation_modified.png/1280px-Pedigree%2C_recombination_and_resulting_IBD_segments%2C_schematic_representation_modified.png){width=100%}

Columns Z0, Z1, and Z2 indicate the probabilities of having IBD of 0, 1, or 2 over the loci, which gives us a way of discriminating between relationship types. 
Ideal parent-offspring has (Z0, Z1, Z2) = (0, 1, 0), i.e. all loci have one allele identical by descent; 
Ideal full sibling = (1/4, 1/2, 1/4), i.e. 25% of loci have 0 alleles IBD, 50% have 1 allele IBD, 25% have 2 alleles IBD; etc. So, we can use Z0, Z1, and Z2 to distinguishing between relationship types. 

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
