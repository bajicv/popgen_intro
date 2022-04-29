# SNP arrays

[SNP array](https://en.wikipedia.org/wiki/SNP_array) is a type of DNA microarray containing designed probes harboring the SNP positions, which is hybridized with fragmented DNA to determine the specific alleles of all SNPs on the array for the hybridized DNA sample. 
SNP arrays can be used to detect polymorphisms within a population. 

In order to create an SNP array, it is necessary to know where are the variable sites in the human genome. This knowledge comes from several big projects that sequenced entire human genomes from different human populations (e.g. HapMap, 1K Genome Project, HGDP, etc.) . Once SNPs and sequences around them were known it was possible to design SNP arrays that would capture only positions in which humans differ from each other. SNP genotyping is much cheaper than sequencing the entire human genomes, and thus it allows us to look into more individuals. 

Many companies nowadays offer different genomewide SNP chips. Here we will get familiar with [Illumina microarray technology](https://www.youtube.com/embed/lVG04dAAyvY) (also known as BeadArray technology) which uses silica microbeads. On the surface of each array, or BeadChip, hundreds of thousands to millions of genotypes for a single individual can be assayed at once. These tiny silica beads are housed in carefully etched microwells and coated with multiple copies of an oligonucleotide probe targeting a specific locus in the genome.

![Illumina beads](https://www.illumina.com/content/dam/illumina-marketing/images/science/v2/web-graphic/multi-sample-array-formats-web-graphic.jpg)

As DNA fragments of a sample pass over the BeadChip, each probe binds to a complementary sequence in the sample DNA, stopping one base before the locus of interest. Allele specificity is conferred by a single base extension that incorporates one of four labeled nucleotides. When excited by a laser, the nucleotide label emits a signal. The intensity of that signal conveys information about the allelic ratio at that locus.

![Illumina genotypes](https://www.illumina.com/content/dam/illumina-marketing/images/technology/microarray/how-microarrays-work-web-graphic.jpg)

!!! note "_What is a rs or Reference SNP number?_"
    __dbSNP Reference SNP (rs or RefSNP) number__ is a locus accession for a variant type assigned by [dbSNP](https://www.ncbi.nlm.nih.gov/snp/). The RefSNP catalog is a non-redundant collection of submitted variants which were clustered, integrated and annotated. RefSNP number is the stable accession regardless of the differences in genomic assemblies. RefSNP numbers facilitate large-scale studies in association genetics, medical genetics, functional and pharmaco-genomics, population genetics and evolutionary biology, personal genomics, and precision medicine. They provide a stable variant notation for mutation and polymorphism analysis, annotation, reporting, data mining, and data integration.


*********************************************************************

# Ascertainment bias

Although the SNP chips have been widely used in human population genetics studies, this was not what they were designed for. The SNP arrays are a useful tool for studying slight variations between whole genomes, and their most important applications are for determining disease susceptibility and for measuring the efficacy of drug therapies designed specifically for individuals. In research, SNP arrays are most frequently used for genome-wide association studies.

The ascertainment and choice of the SNPs included in a given SNP chip are often poorly documented. Using SNPs known from one population to study the genetic diversity of another population is not ideal for popgen studies as such an approach could omit numerous SNPs that might exist in the second population, and which were not taken into account when the SNP chip was designed. Effect of the ascertainment bias has been documented for most of the commercially available SNP array platforms. 

_In the figure below each line represents a chromosome, and each dot represents a variant. Using a set of markers optimized to detect the genetic diversity of population 1 (B) on population 2 underestimates the higher diversity of population 2 (A)._ 

<center>
![Ascertainment bias](https://www.researchgate.net/profile/Luca-Pagani/publication/303388742/figure/fig1/AS:614371868811287@1523489076589/Effect-of-the-ascertainment-bias-documented-for-most-of-the-commercially-available-SNP_W640.jpg){width=70%}

_Figure taken from [Luca Pagani's Research gate](https://www.researchgate.net/publication/303388742_Through_the_layers_of_the_Ethiopian_genome_a_survey_of_human_genetic_variation_based_on_genome-wide_genotyping_and_re-sequencing_data/figures?lo=1)._
</center>

*********************************************************************

# Exercise dataset

The aim of this exercise is to get familiar with a few widely used analyses in human popgen as well as to understand  genetic variation in human populations. For this exercise we are going to work with [HGDP-CEPH Dataset 11](http://www.cephb.fr/en/hgdp_panel.php) submitted by Harvard Genetic Department.

<center>
![HGDP](https://upload.wikimedia.org/wikipedia/commons/a/a9/Worldwide_human_populations_-_frappe_results.png){width=70%}
</center>

This dataset contains __629,443 SNPs__ that were obtained by genotyping __934 unrelated HGDP-CEPH individuals from different human populations across the globe__. The dataset we chose is genotyped on the [Affymetrix Axiom™ Human Origins Array](https://www.thermofisher.com/order/catalog/product/901853#:~:text=Axiom%20Genome%2DWide%20Human%20Origins,%2C%20migration%2C%20and%20natural%20selection.) which was designed to reduce problems of [ascertainment bias](https://www.ncbi.nlm.nih.gov/books/NBK9792/). Read the [detailed technical document](ftp://ftp.cephb.fr/hgdp_supp10/8_12_2011_Technical_Array_Design_Document.pdf) before using these data in order to avoid pitfalls. 


!!! note "__Affymetrix Axiom™ Human Origins Array__"
    To address ascertainment bias problem, the Affymetrix Axiom™ Human Origins Array has been designed. This array analyzes __~630,000 SNPs__ ascertained by resequencing in __13 population-specific panels__, including San, Yoruba, Mbuti Pygmy, French, Sardinian, Han, Cambodian, Mongolian, Karitiana, Papuan, and Bougainville populations from the CEPH-HGDP panel, as well as variants identified from the Neanderthal Genome Project. The array also contains 87,000 SNPs that include sets from mtDNA and the Y chromosome, and variants from standard chips for data comparison purposes.

All the data that we need for the exercises are located here:
`/opt/evop/public/BIOINFORMATICS/popgen_intro/`

We will copy them to our home directory instead of downloading them so we can save some time.

``` bash
mkdir ~/popgen_intro
cd /opt/evop/public/BIOINFORMATICS/popgen_intro/
cp -r HGDP_metainformation.txt Harvard_HGDP-CEPH/all_snp.* Harvard_HGDP-CEPH/sample_all_snp.txt ~/popgen_intro/
cd ~/popgen_intro
```

`Harvard_HGDP-CEPH` directory contains several files explained in `8_12_2011_Harvard_HGDP_readme.txt`. 


| File(s) | Description |
|---|---|
| `all_snp.ped` and <br />`all_snp.map` | These two files are in the `plink` format (used by software such as `fastStructure`, `ADMIXTURE`, `TreeMix`, and many others). <br />They are _"connected"_ and they should be always processed together. <br />So be sure that you keep them in the same directory and that they have the same name with only difference in their file extension name. <br />These files contain only 10% of randomly chosen SNPs (with option `--thin 0.1` in plink) from the original file that was downloaded from HGDP-CEPH webpage. <br />This was done to ensure that we are not using too much memory on the evop-login server and to make sure that we will finish exercises on time ;) |
| `sample_all_snp.txt` | Contains a list of individuals their sex and populations to which they belong. |
| `HGDP_metainformation.txt` | A file in which I collected the meta-information on each of the samples which might be useful for population genetic research, <br />data interpretation, and/or data visualization (e.g. country or origin, longitude, latitude, language...). |

!!! note

    If you would like to download it yourself, you can do it like this:

    ``` bash
    wget ftp://ftp.cephb.fr/hgdp_supp10/*
    tar -xzvf Harvard_HGDP-CEPH.tgz
    ```
    `wget` is the non-interactive network downloader which is used to download files from the server even when the user has not logged on to the system and it can work in the background without hindering the current process.

    A **tarball** is a group of files that are kept together using the `tar` command. Tarballs are common file formats on Linux operating systems, and they are often used for distribution of software/media or backup purposes. Typically they have `.tar` extension, while compressed `.tar` files have `.tgz` or `.tar.gz` extension.

For all data manipulation, such as extracting certain SNPs, chromosomes, individuals, merging datasets, removing SNPs in high LD, filtering out low-quality genotypes, and many others we will use [plink](https://www.cog-genomics.org/plink2). 

