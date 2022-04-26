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
!!! note "Task"
    1. Make directory called `plink_exercise` in `~/popgen_intro`.
    2. Make __binary files__ from `all_snp` __text-format__ plink files and output them in newly created `~/popgen_intro/plink_exercise`.

    ??? note "Tip"
        Check plink manual on how to make [BED file format](https://zzz.bwh.harvard.edu/plink/binary.shtml).

    ??? note "Click for answer"
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

As expected we got 3 files that go together in binary file format in `plink` (i.e. `bed` + `bim` + `fam`), but also `log`. 

Let's take a look at it:

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

If we take a look at `.fam` file, we will see that FID (column 1) and IID (column 2) are the same. However, for popgen studies, FID could be used to store the information about the population from which the individual was sampled. We can find population information for each individual in the metainformation file that we have in `~/popgen_intro/HGDP_metainformation.txt`. 

In order to [update sample information](https://www.cog-genomics.org/plink2/data#update_indiv) in our plink files we should first make a text file with old and new FID and IID, and then provide it to `--update-ids` when making new plink files with updated IDs. 

`--update-ids` expects input with the following four fields:

1. Old FID
2. Old IID
3. New FID
4. New IID

To create this file we can use the command `join` which will join lines of two files on a common field.

!!! tip "__Tip: Process Substitution__"
    __Piping__ (`|`) the stdout of a command into the stdin of another is a powerful technique. But, what if we need to pipe the stdout of multiple commands? This is where __process substitution__ (`<(command)` or `>(command)`) comes in.

    Process substitution feeds the output of a process (or processes) into the stdin of another process. In other words process substitution treats the output (or input) of commands as files. The syntax for process substitution is: `<(list)` or `>(list)` where each list is a command or a pipeline of commands. The effect of process substitution is to make each list act like a file.
    <br /> To substitute a command pipeline for an input file the syntax is: `command ... <(list) ...`
    <br /> To substitute a command pipeline for an output file the syntax is:`command ... >(list) ...`
    <br /> E.g. `diff -u <( ls dir1 | sort)  <( ls dir2 | sort )`


``` bash
cd ~/popgen_intro/plink_exercise/

join -j 1 -o 1.1,1.2,2.2,2.1 <(sort -k1 hgdp.fam) <(sort -k1 ../HGDP_metainformation.txt) > hgdp_updateIDs.txt

plink1.9 --bfile hgdp --update-ids hgdp_updateIDs.txt --make-bed --out hgdp_newFIDs
```


*********************************************************************

## Filtering with plink

We can generate some simple summary statistics using `plink` e.g. rates of missing data in the file, diversity within and between the samples, etc., and use some of those statistics to filter our dataset to include only high-quality entries.

We can use `--missing` to produces __sample-based__ (`.imiss`) and __variant-based__ (`.lmiss`) files with missing data reports.

``` bash
cd ~/popgen_intro/plink_exercise/

plink1.9 --bfile hgdp_newFIDs --missing --out hgdp_newFIDs_miss
```

We can visualize the missing data reports in __R__. Let's download newly created files together with the metainformation file to our local machine, and then run the following R commands to plot them _(remember to specify correct path to files)_:

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

!!! task
    Let's create new plink files in which we will keep only individuals from modern populations, and in which we will remove SNPs and individuals with more than 10% missingness. To do so, we will need to specify which individuals we want to keep (with option `--keep`) or remove (with option `--remove`). One way to make list of individuals which we want to keep is to find all individuals that contain "HGDP" in their name.

    Steps:<br />
    1. Create a new txt file that will contain only lines with individuals that contain "HGDP" in IID column of `hgdp_newFIDs.fam`<br />
    2. Use plink to create new binary files that will have only individuals specified in the txt file you created in the previous step, and remove SNPs and individuals with more than 10% missingness.
    
    ??? task "Click for answers"
    
        ``` bash
        grep HGDP hgdp_newFIDs.fam > ind_keep.txt
        plink1.9 --bfile hgdp_newFIDs --make-bed --mind 0.1 --geno 0.1 --keep ind_keep.txt --out hgdp_newFIDs_modernPops_MindGeno
        ```
!!! question
    How many individuals and how many SNPs were removed?
    ??? question "Click for answers"
        `--keep`: 934 out of 942 people remaining.<br />
        `--mind`: 0 people removed due to missing genotype data.<br />
        `--geno`: 8 variants removed due to missing genotype data.
