#-----Load libraries-----
library(NOISeq)                # Version: 2.48.0
library(SummarizedExperiment)  # Version: 1.34.0
library(TCGAbiolinks)          # Version: 2.32.0
library(dplyr)                 # Version: 1.1.4
library(biomaRt)               # Version: 2.60.1
library(edgeR)                 # Version: 4.2.1
library(ggplot2)               # Version: 3.5.1
library(EDASeq)                # Version: 2.38.0

#-----Load data-----
samples_data <- read.csv("samples_data.tsv", sep = "\t", header = T) # .tsv file created in the Get_data.R script

#------Get reads data-----
exp_data <- GDCquery(project = "TCGA-LIHC",                          # Liver hepatocellular carcinoma
                     data.category = "Transcriptome Profiling",      # Self explanatory 
                     data.type = "Gene Expression Quantification",   # Self explanatory
                     workflow.type="STAR - Counts",                  # Read counts
                     barcode = samples_data$Experiment_id)           # Barcodes used to filter the download files

# DO NOT RUN THIS LINE IF YOU'VE ALREADY DOWNLOADED YOUR DATA
GDCdownload(exp_data)                                              # Download data (will create GDCdata folder unless directory is set)
                                                                   # This folder will contain all the experiments' reads and genes (228 as 26-08-2024)

# Downloaded data has unstranded, stranded first, stranded second, tpm, fpkm and fpkm-uq-untranded reads for every single sample
exp_lvls <- GDCprepare(exp_data, summarizedExperiment = F)         # Reads the downloaded data and prepares it into an R-readable object
dim(exp_lvls[which(duplicated(exp_lvls$gene_name) ==T),])          # There are 1,236 genes with the same name but different gene id (09-08-2024)
dim(exp_lvls[which(duplicated(exp_lvls$gene_id) ==T),]) 

# Format expression files and save the in a tsv file
exp_lvls <- as.matrix(exp_lvls[!which(exp_lvls$gene_name == ""),]) # Remove genes with no name and convert into matrix (removed 4)
rownames(exp_lvls) <- exp_lvls[,1]                                 # Add genes' ENSEMBL ids as row names
exp_lvls <- exp_lvls[,-1]                                          # Remove the ENSEMBL ids column
exp_lvls <- exp_lvls[,1:230]                                       # Only keep unstranded reads for every gene in every sample (there should be 230 columns as there were 228 samples)
write.table(exp_lvls,"RNAseq.tsv",sep='\t',quote=F)                # Save the matrix as a tsv file

#-----Format dataframes and matrices-----
design_exp <- samples_data                                           # Identical object to modify
design_exp$barcode <- strsplit(colnames(exp_lvls[,3:230]), "_") %>%  # Assign the barcode and remove the "unstranded" part of the name
  unlist(.) %>% .[c(F, T)]                     
rownames(exp_lvls) <- sapply(strsplit(rownames(exp_lvls),            # Format the genes' ids and set them as row names
                                      ".",fixed=T),function(x) x[1])

design_exp[which(design_exp$Subtype == "NA"),]$Subtype <- "No specific subtype" # Make subtypes==NA into "No specific subtype"
design_exp <- design_exp[which(!design_exp$Subtype == "No specific subtype"),]  # Remove the "No specific subtype" samples, 45 remain (27-08-2024)

#-----GC content, length and biotype-----
bmart_db <- useEnsembl("ensembl", dataset="hsapiens_gene_ensembl")   # BioMart selection of ENSEMBL database. Also selected homo sapiens dataset

myannot <- getBM(attributes = c("ensembl_gene_id",                   # Gene id
                                "percentage_gene_gc_content",        # Percentage of GC in each gene
                                "gene_biotype",                      # Gene biotype
                                "start_position",                    # Gene start position
                                "end_position",                      # Gene end position
                                "hgnc_id",                           # Unique and approved HGNC gene id
                                "hgnc_symbol"),                      # Unique and approved HGNC gene symbol
                filters = "ensembl_gene_id",                         # Parameter of data retrieved used as filter
                values = rownames(exp_lvls),                         # Filters' values
                mart=bmart_db)                                       # Object of class Mart

myannot$length <- abs(myannot$end_position-myannot$start_position)   # Add gene bp length
dim(myannot)                                                         # 60,280 genes were retrieved out of 60,660
exp_lvls <- exp_lvls[rownames(exp_lvls) %in%                         # Remove the expression of genes that are not in my annotation
                       myannot$ensembl_gene_id,]


# Filter transcripts without annotation
myannot <- myannot[myannot$gene_biotype=="protein_coding",]              # Only keep genes/transcripts that contains an open reading frame (19,925 at 19-08-2024)
myannot <- myannot[!duplicated(myannot$ensembl_gene_id),]                # Remove genes with duplicated ENSEMBL gene id (none removed as expected)
exprots_hgnc <- exp_lvls[rownames(exp_lvls)%in%myannot$ensembl_gene_id,] # Only keep expression data from filtered genes
exprots_hgnc <- exprots_hgnc[!duplicated(rownames(exprots_hgnc)),]       # Remove duplicated rownames
exprots_hgnc <- exprots_hgnc[,which(strsplit(colnames(exprots_hgnc),     # Remove samples with no specific subtype and the first two columns with gene name and gene type
                                             "_") %>% unlist(.) %>%
                                      .[c(F, T)] %in%          
                                      design_exp$barcode)]
dim(exprots_hgnc)                                                        # 19,925 genes expressed in 45 samples (27-08-2024)
exprots_hgnc <- matrix(as.numeric(exprots_hgnc),                         # Make the counts into numeric instead of characters
                       ncol = 45,
                       nrow = 19925,
                       dimnames = list(rownames(exprots_hgnc),
                                       colnames(exprots_hgnc)))             
exprots_hgnc <- exprots_hgnc[which(rowSums(exprots_hgnc) != 0),]         # Filter rows/genes with 0 reads in every single sample (19,164 genes remanining 27-08-2024)
myannot <- myannot[which(rownames(exprots_hgnc) %in%                     # Remove genes that are not present in the annotation object (None were removed 27-08-2024)
                           myannot$ensembl_gene_id),] 

#-----Check for biases-----
# Make object of class Expressionset 
noiseqData <- NOISeq::readData(data = exprots_hgnc,                                       # Matrix/dataframe with all the gene counts for all samples
                       gc = myannot[,c("ensembl_gene_id", "percentage_gene_gc_content")], # Percentage of GC content in each gene
                       biotype = myannot[,c("ensembl_gene_id", "gene_biotype")],          # Biotype of the genes
                       factors = design_exp,                                              # Sort of metadata that contains the conditions (subtypes)
                       length = myannot[,c("ensembl_gene_id", "length")])                 # Bp length of each feature

mycountsbio <- dat(noiseqData, type = "countsbio", factor = "Subtype")                    # Generate expression data with feature annotation object

# Check expression bias values
explo.plot(mycountsbio, plottype = "boxplot", samples = 1:8)      # Exploratory boxplot demonstrating the expression values for each subtype

# Check for low count genes
explo.plot(mycountsbio, plottype = "barplot", samples = 1:8)      # Exploratory barplot showing CPM per subtype

# Mean of log CPM
ggplot(data = exprots_hgnc,                                       # Exploratory histogram portrays mean log CPM to view low CPM values overall
       aes(x=rowMeans(cpm(exprots_hgnc, log = T)))) + 
  geom_histogram(color="black", fill="white", binwidth=.5) + 
  labs(x="Mean of log CPM") +
  geom_vline(xintercept = 0)

# Transcript composition bias
mycd <- dat(noiseqData, type = "cd", norm = FALSE) 
mycd                                                      # View object for better understanding
table(mycd@dat$DiagnosticTest[,"Diagnostic Test"])        # 43 failed, 1 passed (27-08-2024)
explo.plot(mycd,samples=sample(1:ncol(exprots_hgnc),10))  # Data normalization is needed

# GC bias
GCcontent <- dat(noiseqData,  k = 0, type = "GCbias", norm = F, factor = "Subtype")
GCcontent                                                                            # View object for better understanding
par(mfrow=c(2,4))                                                                    # Show the plots for the 8 subtypes
sapply(1:8,function(x) explo.plot(GCcontent, samples = x))                           # Show the GC content for the 8 subtypes

# Length bias
lenbias <- dat(noiseqData, k = 0, type = "lengthbias", factor = "Subtype")
par(mfrow=c(2,4))                                                                    # Show the plots for the 8 subtypes
sapply(1:8,function(x) explo.plot(lenbias, samples = x))                             # Show the length bias for the 8 subtypes

par(mfrow=c(1,1))
# Check for batch effect
my_pca <- dat(noiseqData, type = "PCA", norm = F, logtransf = F)                     # Make PCA
explo.plot(my_pca, samples = c(1,2), plottype = "scores", factor = "Subtype")        # Plot to check if samples group by subtype

#-----Reduce biases-----
# Filter genes with low read count
count_matrix_filtered <- filtered.data(exprots_hgnc, factor = "Subtype", norm = F,   
                                       depth = NULL, method = 1, cpm = 0, p.adj = "fdr")      # 10,153 features are to be kept (29-08-2024)
myannot <- myannot[myannot$ensembl_gene_id %in% rownames(count_matrix_filtered),]             # Keep only the filtered genes (10,049 genes 29-08-2024)
count_matrix_filtered <- count_matrix_filtered[rownames(count_matrix_filtered) %in%           # Remove the genes that are not in my annotation
                                                 myannot$ensembl_gene_id,]
colnames(count_matrix_filtered) <- strsplit(colnames(count_matrix_filtered), "_") %>%         # Format column names so they match with the PhenoData
  unlist(.) %>% .[c(F, T)]

# Create an object containing the new RNA expression set
new_exp_data <- newSeqExpressionSet(                                                          # Create expression set function
  counts = as.matrix(count_matrix_filtered),                                                  # Gene count matrix
  featureData = data.frame(myannot, row.names=myannot$ensembl_gene_id),                       # Features of data, such as GC content, lenght, etc.
  phenoData = data.frame(design_exp, row.names=design_exp$barcode))                           # Sample information assigned

# Correct GC and length biases
gcFull <- withinLaneNormalization(new_exp_data, "percentage_gene_gc_content", which = "full") # Correct for GC composition of reads
lFull <- withinLaneNormalization(gcFull, "length", which = "full")                            # Correct based on length

# Normalize data
normalized_count_matrix <- tmm(normCounts(lFull),                                             # Corrected read counts matrix
                               long = 1000,                                                   # No length correction is applied, already did it
                               lc = 0)                                                        # Same as previous argument

# Create new expression set based on the corrected normalized counts
noiseqData <- NOISeq::readData(data = normalized_count_matrix, factors = design_exp)
mycd <- dat(noiseqData, type="cd", norm=TRUE)                                                 # Generate results of pre-processed data for exploration        
table(mycd@dat$DiagnosticTest[,  "Diagnostic Test"])                                          # 44 PASSED (06-09-2024)

# Check for batch effect now that data has been pre-processed
myPCA <- dat(noiseqData, type = "PCA", norm = T, logtransf = F)                               # Perform PCA
explo.plot(myPCA, samples = c(1,2), plottype = "scores", factor = "Subtype")                  # Visualize PCA (Only seems to cluster clear cell and steatohepatic subtypes)

# Remove possible non identified batch effects 
full_data <- ARSyNseq(noiseqData,                                                             # Data to filter noise from
                      factor = "Subtype",                                                     # Factor of interest
                      batch = FALSE,                                                          # The factor SHALL NOT be used as the batch
                      norm = "n",                                                             # Data has already been normalized
                      logtransf = F)                                                          # False if we want to perform log-transformation of data

myPCA <- dat(full_data,                                                                       # Data to explore
             type = "PCA",                                                                    # PCA exploration
             norm = T,                                                                        # Data has already been normalized
             logtransf = T)                                                                   # Avoid performing a second log-transformation

explo.plot(myPCA, factor = "Subtype")                                                         # Visualize batch removal effect on PCA

#-----Final quality check-----
noiseqData <- NOISeq::readData(data = exprs(full_data),                                       # Expression data
                               gc = myannot[,1:2],                                            # GC content
                               biotype = myannot[,c(1,3)],                                    # Biotype
                               factor = design_exp,                                           # Factors
                               length = myannot[,c(1,8)])                                     # Gene length

# Counts by subtype
mycountsbio <- dat(noiseqData, type = "countsbio", factor = "Subtype", norm=T)                # Gene count by subtype
explo.plot(mycountsbio, plottype = "boxplot",samples=1:5)                                     # Visualization

# GC content bias by subtype
myGCcontent <- dat(noiseqData, k = 0, type = "GCbias", factor = "Subtype", norm=T)            # Explore data based on GC bias
par(mfrow=c(2,4))                                                                             # Showcase every plot
sapply(1:8,function(x) explo.plot(myGCcontent, samples = x))                                  # Visualize all plots

# Length bias by subtype
mylenBias <- dat(noiseqData, k = 0, type = "lengthbias", factor = "Subtype",norm=T)           # Explore data based on length bias
sapply(1:8,function(x) explo.plot(mylenBias, samples = x))                                    # Visualize plots

par(mfrow=c(1,1))                                                                             # Set one plot at the tome again

#-----Duplicates and save data-----
sum(duplicated(design_exp$Patient_id))                                                        # 4 patients have multiple samples (29-08-2024)
duplicates <- design_exp$Patient_id[duplicated(design_exp$Patient_id)]                        # Extract patients' ids
duplicates <- lapply(duplicates,function(x) design_exp$barcode[design_exp$Patient_id==x])     # Get the duplicates' barcodes

final  <- exprs(full_data)                                                                    # Get expression values from all samples, including duplicates
duplis <- final[, colnames(final) %in% unlist(duplicates)]                                    # Retrieve expression values only from duplicates
prefi  <- final[, !colnames(final) %in% unlist(duplicates)]                                   # Retrieve expression values from non-duplicates

duplicates <- do.call(cbind,lapply(duplicates,function(x)                                     # Get the average expression for the duplicates and set it as the expression value
  rowMeans(duplis[,colnames(duplis)%in%x])))

colnames(duplicates) <- design_exp$Patient_id[duplicated(design_exp$Patient_id)]              # Set the column names 
colnames(prefi) <- substr(colnames(prefi),1,19)                                               # Set column names

final <- cbind(prefi,duplicates)                                                              # Combine the averaged duplicates with the rest of the samples
dim(final)                                                                                    # 10,049 genes in 41 unique samples (06-09-2024)
final <- final[,order(match(colnames(final),samples_data$Experiment_id))]                     # Sort the expression matrix

write.table(final, "RNAseqnormalized.tsv", sep='\t', quote=F)                                 # Write table
