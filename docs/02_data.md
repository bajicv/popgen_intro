# Exercise dataset

For this exercise we are going to work with [HGDP-CEPH Dataset 11](http://www.cephb.fr/en/hgdp_panel.php) submitted by Harvard Genetic Department.

<center>
![HGDP](https://upload.wikimedia.org/wikipedia/commons/a/a9/Worldwide_human_populations_-_frappe_results.png){width=70%}
</center>

This dataset contains SNPs from different human populations across the globe. The aim of this exercise is to get familiar with a few widely used analyses in human popgen as well as to understand  genetic variation in human populations. The dataset we chose is genotyped on the Affymetrix Axiom Human Origins Array which was designed to reduce problems of [ascertainment bias](https://www.ncbi.nlm.nih.gov/books/NBK9792/). Effect of the ascertainment bias has been documented for most of the commercially available SNP array platform. 

_In the figure below each line represents a chromosome, and each dot represents a variant. The usage on population 2 of a set of markers optimized to detect the genetic diversity of population 1 (B) underestimates the higher diversity of population 2 (A)._ 

<center>
![Ascertainment bias](https://www.researchgate.net/profile/Luca-Pagani/publication/303388742/figure/fig1/AS:614371868811287@1523489076589/Effect-of-the-ascertainment-bias-documented-for-most-of-the-commercially-available-SNP_W640.jpg){width=70%}

Figure taken from [Luca Pagani Research gate](https://www.researchgate.net/publication/303388742_Through_the_layers_of_the_Ethiopian_genome_a_survey_of_human_genetic_variation_based_on_genome-wide_genotyping_and_re-sequencing_data/figures?lo=1)
</center>

*********************************************************************

## Affymetrix Axiom Human Origins Array

Although the __SNP arrays__ have been widely used in human population genetics studies, this was not what they were designed for. The SNP arrays are a useful tool for studying slight variations between whole genomes, and their most important clinical applications are for determining disease susceptibility and for measuring the efficacy of drug therapies designed specifically for individuals. In research, SNP arrays are most frequently used for genome-wide association studies.

The ascertainment and choice of the SNPs are not well documented, and the tagging design is not ideal for population studies (__ascertainment bias__). To address this problem, the __Affymetrix Axiom Human Origins Array__ has been designed. This analyzes ~630,000 SNPs ascertained by resequencing in 13 population-specific panels, including San, Yoruba, Mbuti Pygmy, French, Sardinian, Han, Cambodian, Mongolian, Karitiana, Papuan, and Bougainville populations from the CEPH-HGDP panel, as well as variants identified from the Neanderthal Genome Project. The array also contains 87,000 SNPs that include sets from mtDNA and the Y chromosome, and variants from standard chips for data comparison purposes.

Dataset 11 contains __629,443 SNPs__ that were obtained by genotyping __934 unrelated HGDP-CEPH individuals__ with Affymetrix Axiom Human Origins Array Plate, and merging the genotypes with data from _Neandertal_, _Denisova_ and _chimpanzee_. 

The SNP data are divided among 14 partially overlapping datasets, 13 of which are of value for analysis of different population genetics scenarios. 

For each of __datasets 1-12__, SNPs were discovered as heterozygotes by whole-genome shotgun sequencing of a different HGDP-CEPH individual of known ancestry, as per [Keinan et al. 2007](https://www.nature.com/articles/ng2116). 

__Dataset 13__ contains heterozygote SNPs for each of which a random _Denisovan_ allele matches that of _chimpanzee_, and the random _San_ allele is derived. 

__Dataset 14__, which is valuable for studying population structure using a maximum number of SNPs, _does not allow demographic modeling_. This dataset combines __all SNPs__ along with an __additional 87,044 SNPs chosen to allow haplotype inference at mitochondrial DNA and the Y chromosome__, and to provide overlap with previous Affymetrix and Illumina genotyping arrays so that users can merge the data available here with previously published datasets. 

Read the [detailed technical document](ftp://ftp.cephb.fr/hgdp_supp10/8_12_2011_Technical_Array_Design_Document.pdf) before using these data in order to avoid pitfalls. The array was developed by David Reich and colleagues in collaboration with Affymetrix for the purpose of generating data with clearly documented ascertainment.

All the data that we need for the exercises are located here:
`/opt/evop/public/BIOINFORMATICS/popgen_intro/`

We will copy them to our home directory instead of downloading them so we can save some time.
``` bash
mkdir ~/popgen_intro
cd /opt/evop/public/BIOINFORMATICS/popgen_intro/
cp -r HGDP_metainformation.txt Harvard_HGDP-CEPH/all_snp.* Harvard_HGDP-CEPH/sample_all_snp.txt ~/popgen_intro/
cd ~/popgen_intro
```

If you would like to download it yourself, you can do it like this:<br />
`wget ftp://ftp.cephb.fr/hgdp_supp10/*`

A **tarball** is a group of files that are kept together using the `tar` command. Tarballs are common file formats on Linux operating systems, and they are often used for distribution of software/media or backup purposes. Typically they have `.tar` extension, while compressed `.tar` files have `.tgz` or `.tar.gz` extension.

``` bash
tar -xzvf Harvard_HGDP-CEPH.tgz
```

`Harvard_HGDP-CEPH` directory contains several files explained in `8_12_2011_Harvard_HGDP_readme.txt`. 

- `all_snp.ped` and `all_snp.map`<br />
  are in the `plink` format (used by software such as `ADMIXTURE`, `TreeMix`, and many others). These two files are _"connected"_ and they should be always processed together. So be sure that you keep them in the same directory and that they have the same name with only difference in their file extension name.

- `sample_all_snp.txt`<br /> 
  contains a list of individuals their sex and populations to which they belong.

- `HGDP_metainformation.txt` 
  is a file in which I collected the meta-information on each of the samples which might be useful for population genetic research, data interpretation, and/or data visualization (e.g. country or origin, longitude, latitude, language...).

For all data manipulation, such as extracting certain SNPs, chromosomes, individuals, merging datasets, removing SNPs in high LD, filtering out low quality genotypes, and many others we will use [plink](https://www.cog-genomics.org/plink2).
