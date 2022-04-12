# ADMIXTURE



[ADMIXTURE](https://www.genetics.ucla.edu/software/admixture/) is a software that we can use to identify ancestry components shared between individuals in a set of populations. It takes `plink` format input files.


*********************************************************************

### Pruning

Before running ADMIXTURE, it is advised to prune the dataset. We will use `plink` to prune the dataset by excluding the SNPs which are in Linkage. The resulting file will have fewer SNPs, and the computation will be faster. We will also only keep SNPs that are on the autosomes (`--autosome`)

<br />To prune for LD we will use these settings previously recommended for HO array `--indep-pairwise 200 25 0.4`:
<br />window size in SNPs = 200;
<br />number of SNPs to shift the window at each step = 25;
<br />r^2 threshold (pairwise SNP-SNP metric) = 0.4;

<br />Using these parameters specified above we will 
<br />a) consider a window of 200 SNPs,
<br />b) calculate LD between each pair of SNPs in the window,
<br />b) remove one of a pair of SNPs if the LD is greater than 0.4,
<br />c) shift the window 25 SNPs forward and repeat the procedure.


```{bash, echo=T, eval=F}
plink1.9 --bfile hgdp_newFIDs_modernPops_MindGeno --indep-pairwise 200 25 0.4 --out prune_hgdp
plink1.9 --bfile hgdp_newFIDs_modernPops_MindGeno --extract prune_hgdp.prune.in --recode12 --make-bed --autosome --out hgdp_pruned_autosome
```

How many SNPs are left after pruning?

*****************************************************

ADMIXTURE demands a bit more computational power and time than previous analysis we did, and we will not be able to run it on `evop-login` server. So, instead I will explain how to run ADMIXTURE and then we will focus on visualizing the output of the ADMIXTURE in R. 


*************************************************************

### Running ADMIXTURE

We can use the following commands to run ADMIXTURE for each K. One of the K values will have better supported than the others and we can find it by searching for the one K which has the lowest __cross-validation (cv) error__. And finaly, we will use the best K as the best representation of the actual data variation and visualize it.

```
for K in 6 7 8 9 
do
  admixture -s $RANDOM --cv hgdp_pruned_autosome.ped $K | tee log.K${K}.out;
done
```

This will create three output files for each K: 
<br />1.`.out` 
<br />1.`.P`
<br />1.`.Q`

*************************************************************

### Cross-validation error to identify the best value of K

We can find the best supported K by looking at the `.out` files in a section where is CV error reported. We can take a look at the file with command `less` and then we will see that we can simply grep pattern "CV" from all our log files and store in new textfile like this:

```{bash, eval=F, echo=T}
grep -h CV log*out | awk -F'[(|=|)|\t| ]' '{print $5 "\t" $7}' > CV.txt
```

Now we can download thsi file to our local machines and plot it in R.
```{r, CV_plot, eval=T, echo=T, fig.align='center'}
library(data.table)
library(dplyr)
library("ggplot2")

setwd("~/Nowick_Lab/Master_course_Bioinformatics_for_Biologists/SoSe2020/popge_intro/")
cv <- fread("CV.txt", col.names = c("K","cv"))   

ggplot(cv, aes(group=1, x=K, y=cv)) + 
  geom_point() +
  geom_line() +
  ggtitle("CV error for each K")
```

When performing ADMIXTURE we should also do several runs for each K. For simplicity reasons in our example we did only one run for each K. When multiple runs are present for each K we can determine which one of them has the highest likelihood, and then we can use it as the representative run for the given K.

*************************************************************

### Plotting ADMIXTURE results

For simplicity reasons, we will plot the results of the admixture on chromosome 6. In reality, one would run admixture on all autosomal chromosomes for multiple K's and with multiple runs for each K.

First copy files, and download them to your local machine, and then visualize results in R.
```{bash, eval=F, echo=T}
cp -r /opt/evop/public/BIOINFORMATICS/popgen_intro/files_4_admixture ~/popgen_intro/
```

And now once you downloaded the files you can plot the results in R:

```{r ADMIXTURE_plot, eval=T, echo=T, fig.width=30, fig.align='center'}
library(data.table)
library(dplyr)
library(ggplot2)
library(tidyr)

setwd("~/Nowick_Lab/Master_course_Bioinformatics_for_Biologists/SoSe2020/popge_intro/")

# read in the necessary files
q_file <- read.table("hgdp_pruned_chr6.9.Q")
fam_file <- read.table("hgdp_pruned_chr6.fam", col.names = c("Population", "IID", "PID","MID","Sex","P"))
popGroups <- read.table("HGDP_metainformation.txt", header = T)

# Prepare the data by adding new columns that contain information abut Region and Population
q_file_ID <- cbind(fam_file[,1:2], q_file)
q_file_ID_Region <- left_join(q_file_ID, popGroups[,c("Population","Region")], by="Population")

# Tranforme the data before plotting
plot_data <- q_file_ID_Region %>% 
  mutate(id = row_number()) %>% 
  gather('pop', 'prob', V1:V9) %>% 
  group_by(id) %>% 
  mutate(likely_assignment = pop[which.max(prob)],
         assingment_prob = max(prob)) %>% 
  arrange(likely_assignment, desc(assingment_prob)) %>% 
  ungroup() %>% 
  mutate(id = forcats::fct_inorder(factor(id)))

# Make a plot
p <- ggplot(plot_data, aes(id, prob, fill = pop)) +
  geom_col() +
  theme(strip.text.y = element_text(angle=90),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        strip.text.x = element_text(angle=90),
        axis.ticks.x = element_blank()) +
  facet_grid(~Population, scales = 'free', space = 'free')

p

```

Type `ggsave(plot = p, width = 30, height = 4,filename = "admixture_plot_HGDP_chr6_ordered_ancestry.pdf")` to save the plot as a pdf.

Look at patterns across populations. Do they follow a geographic structure? Is there a sign of Admixture?

Plotting admixture results with R might be tricky when many K's and Runs are present, but we can use [PONG](https://github.com/ramachandran-lab/pong) to visualize our ADMIXTURE runs and to summarize them. This program will take care of the colors and order of individuals and populations for us. It's great! If you ever have to use ADMIXTURE I highly recommend it.
