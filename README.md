# new_rna_seq_analysis
tutorial on rna seq analysis, focusing on quality control, differential expression and gene set testing


### Read data into R

small_counts <- read.table('C:/Users/Tony Isebe/Desktop/RNASEQANALYSIS/small_counts.txt', header = TRUE)
print(small_counts)
dim(small_counts)
##getting sample1 only

small_counts$Sample_1

small_counts[,1,2]

small_counts[,c('Sample_1','Sample_2')]

log(small_counts)
sum(small_counts$Sample_1)

###sum the ounts for each sample

sample_sums=apply(small_counts, MARGIN = 1, sum)
print(sample_sums)

??MARGIN
sample_sums=apply(small_counts, MARGIN = 3, sum)
print(sample_sums)


ResultsTable_small <- read.table('C:/Users/Tony Isebe/Desktop/RNASEQANALYSIS/ResultsTable_small.txt', header = TRUE)
print(ResultsTable_small)


##structure of ResultsTable_small

str(ResultsTable_small)


##sort data from largest to smallest

sort(ResultsTable_small$logFC)

##sort from decreasing

sort(ResultsTable_small$logFC, decreasing=TRUE)


### load all databases needed

library(edgeR)
library(limma)
library(Glimma)
library(gplots)
library(org.Mm.eg.db)
library(RColorBrewer)

if(!requireNamespace('BiocManager', quietly = TRUE))
  install.packages('BiocManager')
##else
library(BiocManager)
BiocManager::install('Glimma')


###read the data to R

Seqdata <- read.delim('C:/Users/Tony Isebe/Desktop/RNASEQANALYSIS/GSE60450_Lactation-GenewiseCounts.txt', stringsAsFactors = FALSE)

###read sample information to R

sampleinfo <- read.delim('C:/Users/Tony Isebe/Desktop/RNASEQANALYSIS/SampleInfo.txt')


### visualize the data

head(Seqdata)

## rows and columns
dim(Seqdata)
View(Seqdata)
sampleinfo

##create a data object with only 12 samples

countdata <- Seqdata[, -(1:2)] ##removes the first two columns
countdata
head(countdata)

### store Entrez geneID as rownames

rownames(countdata) <- Seqdata[,1]
head(countdata)

#### check the colnames

colnames(countdata)

##shorten the names using substr to have 7 characters

colnames(countdata) <- substr(colnames(countdata), start = 1, stop = 7)
head(countdata)

table(colnames(countdata)==sampleinfo$SampleName)

###Filter lowly expressed genes

###start with obtaiing counts per million(cpm)

myCPM <- cpm(countdata)
head(myCPM)

##find values that are greater than 0.5

thresh <- myCPM>0.5
head(thresh)
##summary of how many TRUEs are there

table(rowSums(thresh))

###keep genes that have at least 2 TRUEs in the 12 samples

keep <- rowSums(thresh) >=2

##subset the data to keep the highly expressed genes

counts.keep <- countdata[keep, ]
summary(keep)
dim(counts.keep)
###A CPM of 0.5 is used as it corresponds to a count of 10-15 for the library sizes in this data set. If the count is any smaller, it is considered to be very low, indicating that the associated gene is not expressed in that sample. A requirement for expression in two or more libraries is used as each group contains two replicates. This ensures that a gene will be retained if it is only expressed in one group. Smaller CPM thresholds are usually appropriate for larger libraries. As a general rule, a good threshold can be chosen by identifying the CPM that corresponds to a count of 10, which in this case is about 0.5. You should filter with CPMs rather than filtering on the counts directly, as the latter does not account for differences in library sizes between samples.

# Let's have a look and see whether our threshold of 0.5 does indeed correspond to a count of about 10-15
# We will look at the first sample
library(ggplot2)
plot(myCPM[,2], countdata[,2]) + abline(v=median(2000),col='red')+abline(h=median(5000),col='blue')

## store count data in DEG

y <- DGEList(counts.keep)
y
names(y)
y$samples

###check number of libraries for each read

y$samples$lib.size


##plot library sizes as barplots to check for consistency
# The names argument tells the barplot to use the sample names on the x-axis
# The las argument rotates the axis names

barplot(y$samples$lib.size,names=colnames(y),las=2)

##add title

title('Barplot of library sizes')

###examine distribution of counts using the log function
###get log counts per million
logcounts <- cpm(y, log = TRUE)

##check distribution of counts using the boxplot
boxplot(logcounts,xlab='',ylab='Log2 Counts per million',las=2)

###add a blue line across the median

abline(v=median(logcounts),col='blue')

##Multidimensional Scaling Plots which helps in visualization of principal component analysis and identifies sources of variation ina  given dataset

plotMDS(y)

###specify option to let us plot two plots side by side

par(mfrow=c(1,2))
##setup color scheme for the celltype
##how many celltypes and their corresponding format of storage

levels(sampleinfo$CellType)
## we choose purple for basal and orange for luminal

col.cell <- c('purple','orange')[sampleinfo$CellType]
data.frame(sampleinfo$CellType,col.cell)

##redo the MDS with the celltype coloring

plotMDS(y,col=col.cell)
##add a legend to know which color corresponds to which type

legend("topright",fill=c("purple","orange"),legend=levels(sampleinfo$CellType))
##add title
title('celltype')

##similar for status,
levels(sampleinfo$Status)
col.status <- c('blue','red','green')[sampleinfo$Status]
data.frame(sampleinfo$Status,col.status)
plotMDS(y,col=col.status)

##add a legend
legend('bottomright',fill = c('green','blue','red'),legend = levels(sampleinfo$Status))
title('status')

##Hierarchical clustering using Heatmaps

##Identify 500 most variable genes

##estimate variance for each row in logcounts matrix

var_genes <- apply(logcounts,1,var)
head(var_genes)
##get gene names for for top 500 most variable genes

select_var <- names(sort(var_genes, decreasing = TRUE))[1:500]
head(select_var)

##subset logcounts matrix

highly_variable_cpm <- logcounts[select_var,]
dim(highly_variable_cpm)
head(highly_variable_cpm)

##getting nice colors

mypalette <- brewer.pal(11,'RdYlBu')
morecols <- colorRampPalette(mypalette)

##setup color vector for celltype variable

col.cell <- c('purple','orange')[sampleinfo$CellType]

###plot the heatmap
heatmap.2(highly_variable_cpm,col = rev(morecols(50)),trace = 'none',main = 'Top 500 Most Variable Genes across Samples',ColSideColors = col.cell,scale = 'row')

##save the heatmap

png(file='High_var_genes.heatmap.png')
heatmap.2(highly_variable_cpm,col=rev(morecols(50)),trace = 'none',main = 'Top 500 Most Variable Genes across Samples',ColSideColors = col.cell,scale = 'row')
dev.off()


###change color scheme to PiYG
mypalette <- brewer.pal(11,'PiYG')
morecols <- colorRampPalette(mypalette)

##setup color vector for celltype variable

col.cell <- c('purple','orange')[sampleinfo$CellType]

###plot the heatmap
heatmap.2(highly_variable_cpm,col = rev(morecols(50)),trace = 'none',main = 'Top 500 Most Variable Genes across Samples',ColSideColors = col.cell,scale = 'row')

##save the heatmap

png(file='High_var_genes.heatmap.png')
heatmap.2(highly_variable_cpm,col=rev(morecols(50)),trace = 'none',main = 'Top 500 Most Variable Genes across Samples',ColSideColors = col.cell,scale = 'row')
dev.off()

##Identify 500 least variable genes

##estimate variance for each row in logcounts matrix

var_genes <- apply(logcounts,1,var)
head(var_genes)
##get gene names for for top 500 most variable genes

select_var <- names(sort(var_genes, decreasing  = FALSE))[1:500]
head(select_var)

##subset logcounts matrix

least_variable_cpm <- logcounts[select_var,]
dim(least_variable_cpm)
head(least_variable_cpm)

##getting nice colors

mypalette <- brewer.pal(11,'RdYlBu')
morecols <- colorRampPalette(mypalette)

##setup color vector for celltype variable

col.status <- c('blue','green','red')[sampleinfo$Status]

###plot the heatmap
heatmap.2(least_variable_cpm,col = rev(morecols(50)),trace = 'none',main = 'Top 500 Least Variable Genes across Samples',ColSideColors = col.status,scale = 'row')

##save the heatmap

png(file='Least_var_genes.heatmap.png')
heatmap.2(least_variable_cpm,col=rev(morecols(50)),trace = 'none',main = 'Top 500 Least Variable Genes across Samples',ColSideColors = col.cell,scale = 'row')
dev.off()

###NORMALIZATION FOR ANY FORM OF BIAS

###done to eliminate bias in the composition of libraries

##apply normalization to DEGlist
y <- calcNormFactors(y) ##updates normalization factors

y$samples
par(mfrow=c(1,2))
plotMD(logcounts,column = 7)
abline(h=0,col='blue')
plotMD(logcounts,column = 11)
abline(h=0,col='red')
##mean difference plots show average expression (x-axis) against logfold changes (y-axis)

save(group,y,logcounts,sampleinfo,file="day1objects.Rdata")

##DIFFERENTIAL EXPRESSION

load(day1objects.Rdata)
group
sampleinfo$CellType
design <- model.matrix(~ 0 + group)

y <- calcNormFactors(y)
y$samples
par(mfrow=c(1,2))
plotMD(logcounts,column = 7)
abline(h=0,col="grey")
plotMD(logcounts,column = 11)
abline(h=0,col="grey")


par(mfrow=c(1,2))
plotMD(y,column = 7)
abline(h=0,col="grey")
plotMD(y,column = 11)
abline(h=0,col="grey")

# save(group,y,logcounts,sampleinfo,file="day1objects.Rdata")
# load("day1objects.Rdata")
# objects()

# Look at group variable again
group

# We want to highlight the significant genes. We can get this from decideTests.
par(mfrow=c(1,2))
plotMD(fit.cont,coef=1,status=summa.fit[,"B.PregVsLac"], values = c(-1, 1))

# save(group,y,logcounts,sampleinfo,file="day1objects.Rdata")
group

