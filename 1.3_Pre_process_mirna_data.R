#-----Load libraries-----
library(TCGAbiolinks)          # Version: 2.32.0
library(biomaRt)               # Version: 2.60.1
library(NOISeq)                # Version: 2.48.0
library(edgeR)                 # Version: 4.2.1
library(EDASeq)                # Version: 2.38.0
library(tidyverse)             # Version: 2.0.0

#-----Load data-----
samples_data <- read.csv("samples_data.tsv", sep = "\t", header = T)   # .tsv file created in the Get_data.R script

#------Get miRNA data-----
mirna_data <- GDCquery(project = "TCGA-LIHC",                          # Liver hepatocellular carcinoma
                       data.category = "Transcriptome Profiling",      # Self explanatory 
                       data.type = "miRNA Expression Quantification",  # Self explanatory
                       barcode = samples_data$Experiment_id)           # Barcodes used to filter the download files

# DO NOT RUN THIS LINE IF YOU'VE ALREADY DOWNLOADED YOUR DATA
GDCdownload(mirna_data)                                                # Download data (will create GDCdata folder unless directory is set)
                                                                       # This folder will contain all the experiments' reads and genes (228 as 26-08-2024)

# Downloaded data has the read count, reads per million, and if they are cross mapped or not
mir_data <- GDCprepare(mirna_data)                                     # Reads the downloaded data and prepares it into an R-readable object
rownames(mir_data) <- mir_data$miRNA_ID                                # Set the miRNA id as the row names
mir_data <- mir_data[,grep("read_count",colnames(mir_data))]           # Only keep the "read_count" columns for each sample
colnames(mir_data) <- gsub("read_count_", "", colnames(mir_data))      # Format column names so we only keep the sample name without the "read_count" part
mir_data <- mir_data[which(!rowSums(mir_data) == 0),]                  # Remove miRNAs that were not detected in any single sample
dim(mir_data)                                                          # 1,490 miRNA to 228 samples (06-09-2024)

write.table(mir_data, "miRNAseq.tsv", sep='\t', quote=F)               # Write table with reads

#-----Format dataframes and matrices-----
reduced_id <- substr(colnames(mir_data),1,19)                          # Reduced id data up to portion and analyte 
duplicated_id <- reduced_id[duplicated(reduced_id)]                    # Check for duplicated id's. None were detected (06-09-2024)

design_exp <- samples_data                                             # New object for manipulation
design_exp <- designExp[order(match(design_exp$Experiment_id,          # Sort object by the experiment id
                                    substr(colnames(mirna_data),
                                           1,19))),]
design_exp$barcode <- colnames(mir_data)                                        # Add the reduced experiment id as a barcode
design_exp[which(is.na(design_exp$Subtype)),]$Subtype <- "No specific subtype"  # Make subtypes==NA into "No specific subtype"
design_exp <- design_exp[which(!design_exp$Subtype == "No specific subtype"),]  # Remove the "No specific subtype" samples, 45 remain (19-09-2024)
mir_data <- mir_data[,which(colnames(mir_data) %in% design_exp$barcode)]        # Remove samples from expression matrix
mir_data <- mir_data[which(rowSums(mir_data) > 0),]                             # Remove features with 0 reads across all samples

#-----GC content, length and biotype-----
bmart_db <- useEnsembl("ensembl", dataset="hsapiens_gene_ensembl")     # BioMart selection of ENSEMBL database. Also selected homo sapiens dataset
myannot <- getBM(attributes = c("ensembl_gene_id",                     # Get ENSEMBL gene id
                                "percentage_gene_gc_content",          # Get GC percentage for each gene
                                "mirbase_id",                          # Get the id from miRBase (most of the do not have one)
                                "start_position",                      # miRNA start position
                                "end_position"),                       # miRNA end position
                 mart = bmart_db)                                      # Object of class Mart

myannot <- myannot[myannot$mirbase_id%in%rownames(mir_data),]          # Only keep the annotation of the genes that are in the miRNA expression data
myannot$length <- myannot$end_position-myannot$start_position          # Get the transcript length

# Dealing with mirna duplicates is kind of confusing
# We will remove those that have the same mirbase_id and GC content just to make sure they are the same 
dim(myannot)                                                           # 1,433 genes (19-09-2024)
myannot <- myannot[!duplicated(myannot[,2:3]),]                        # 1,344 genes remain (19-09-2024)
myannot <- myannot[!duplicated(myannot$mirbase_id),]                   # Remove samples that have the same mirbase id (15 were removed 19-09-2024)
myannot <- myannot[!duplicated(myannot[,2:3]) &                        # Remove genes that have same ensembl_id and gc content but different mirbase_id
                     !duplicated(myannot$ensembl_gene_id),]

#-----Check for biases-----
# There should be little or no length bias due to the nature of micro RNAs, however, it is necessary to check either way
# I am not so sure about GC content bias and features with lower counts. Also, I need to check what does miRNA means in this case,
# some literature dumps every small RNA molecule in the same bin while others segregate based on origin and function

noiseqData <- readData(data = mir_data,                                # Create object based on counts matrix
                       factor = design_exp,                            # Assign conditions of interest
                       gc = myannot[,c(3,2)],                          # GC content for each feature
                       length = myannot[,c(3,6)])                      # Length of every feature

mycountsbio <- dat(noiseqData, type = "countsbio", factor = "Subtype") # Generate expression data with feature annotation object

# Check expression values
explo.plot(mycountsbio, plottype = "boxplot", samples = 1:8)           # Exploratory boxplot demonstrating the expression values for each subtype

# Check for low count genes
explo.plot(mycountsbio, plottype = "barplot", samples = 1:8)           # Exploratory barplot showing CPM per subtype

# Mean of log CPM
ggplot(data = mir_data,                                                # Exploratory histogram portrays mean log CPM to view low CPM values overall
       aes(x=rowMeans(cpm(mir_data, log = T)))) + 
  geom_histogram(color="black", fill="white", binwidth=.5) + 
  labs(x="Mean of log CPM") +
  geom_vline(xintercept = 0)

# Check GC bias
GCcontent <- dat(noiseqData, k = 0, type = "GCbias", factor = "Subtype")
par(mfrow=c(2,4))                                                      # Show the plots for the 8 subtypes
sapply(1:8,function(x) explo.plot(GCcontent, samples = x))             # Show the GC content for the 8 subtypes
# Results show R2 values among every subtype and their mean expression values, however 
# the p-values are only significant (p<0.05) in the sarcomatoid subtype (19-09-2024)

# Check length bias
mylenBias <- dat(noiseqData, k = 0, type = "lengthbias", factor = "Subtype",norm=T) # Explore data based on length bias
sapply(1:8,function(x) explo.plot(mylenBias, samples = x))                          # Visualize plots
# Results show very high R2 values, but only three subtypes' p-values are significant (p<0.05) (19-09-2024)

par(mfrow=c(1,1))                                                                   # Set one plot at the tome again

# Check for batch effect
my_pca <- dat(noiseqData, type = "PCA", norm = F, logtransf = F)                    # Make PCA
explo.plot(my_pca, samples = c(1,2), plottype = "scores", factor = "Subtype")       # Plot to check if samples group by subtype

# Transcript composition bias
mycd <- dat(noiseqData, type = "cd", norm = FALSE)
table(mycd@dat$DiagnosticTest[,"Diagnostic Test"])                     # 44 failed (19-09-2024)
explo.plot(mycd,samples=sample(1:ncol(mir_data),10))                   # Data normalization is needed

#-----Reduce biases-----
# Filter genes with low read counts
count_matrix_filtered <- filtered.data(mir_data, factor = "Subtype", norm = F,       # 251 features out of 1,353 are to be kept (20-09-2024)   
                                       depth = NULL, method = 1, cpm = 0)
myannot <- myannot[myannot$mirbase_id %in% rownames(count_matrix_filtered),]         # Keep only the filtered genes (246 genes 20-09-2024)
count_matrix_filtered <- count_matrix_filtered[rownames(count_matrix_filtered) %in%  # Remove the genes that are not in my annotation
                                                 myannot$mirbase_id,]

# Create an object containing the new RNA expression set
new_exp_data <- newSeqExpressionSet(                                   # Create expression set function
  counts = as.matrix(count_matrix_filtered),                           # Gene count matrix
  featureData = data.frame(myannot, row.names=myannot$mirbase_id),     # Features of data, such as GC content, length, etc.
  phenoData = data.frame(design_exp, row.names=design_exp$barcode))    # Sample information assigned







