#-----Load libraries-----
library(SummarizedExperiment)  # Version: 1.34.0
library(TCGAbiolinks)          # Version: 2.32.0
library(dplyr)                 # Version: 1.1.4

#-----Expression data-----
exp_data <- GDCquery(project = "TCGA-LIHC",                        # Liver hepatocellular carcinoma
                     data.category = "Transcriptome Profiling",    # Self explanatory 
                     data.type = "Gene Expression Quantification", # Self explanatory
                     workflow.type="STAR - Counts")                # Read counts

exp_data <- getResults(exp_data)                                   # Get results from query (must be a GDCquery object)
exp_ids <- substr(exp_data$cases,1,19)                             # Extract all sample ids. Remove institute code
length(exp_ids)                                                    # There are 424 samples  (06-08-2024)

#-----Methylation data-----
meth_data <-  GDCquery(project = "TCGA-LIHC",                      # Liver hepatocellular carcinoma
                     data.category = "DNA Methylation",            # Self explanatory
                     platform="Illumina Human Methylation 450")    # Self explanatory

meth_data <- getResults(meth_data)                                 # Get results from query (must be a GDCquery object)
meth_ids <- substr(meth_data$cases,1,19)                           # Extract all sample ids. Remove institute code
length(meth_ids)                                                   # There are 1,290 samples  (26-08-2024)

#-----Micro RNA data-----
mirna_data <- GDCquery(project = "TCGA-LIHC",                         # Liver hepatocellular carcinoma
                       data.category = "Transcriptome Profiling",     # Self explanatory
                       data.type = "miRNA Expression Quantification") # Self explanatory

mirna_data <- getResults(mirna_data)                                  # Get results from query (must be a GDCquery object)
mirna_ids <- substr(mirna_data$cases,1,19)                            # Extract all sample ids. Remove institute code
length(mirna_ids)                                                     # There are 425 samples  (26-08-2024)

#-----Extract shared data sample ids-----
samples_ids <- intersect(exp_ids, meth_ids) %>%                       # Intersect ids between the expression and methylation data ids
  intersect(., mirna_ids)                                             # Intersect the ids with the miRNA data ids
length(samples_ids)                                                   # There are 410 samples (26-08-2024)

#-----Format into dataframe-----
tissues <- c()                                                        # Make empty data frame
for(sample in samples_ids) {                                          # Iterate over every single sample id
  tissues <- append(tissues, meth_data[which(                         # Extract and append each sample type
    sample %in% 
      substr(meth_data$cases,1,19)),]$sample_type)
}

samples_data <- data.frame(cbind(samples_ids,                         # Make dataframe with first column being the whole experiment id
                                 substr(samples_ids,1,12),            # Second column has up to the patient id
                                 tissues),                            # Third column is sample types   
                           row.names = samples_ids)                   # The row names are the same as the whole experiment id

colnames(samples_data) <- c("Experiment_id", 
                            "Patient_id","Sample_type")               # Change column names

samples_data <- samples_data[which(samples_data$Sample_type == "Solid Tissue Normal"),] # None were removed (26-08-2024)                                  # Remove samples tagged as "Solid Tissue Normal"
 

#-----HCC variants-----
sample_subtypes <- TCGAquery_subtype(tumor = "LIHC")                                 # Table with hepatocellular carcinoma variants
sum(samples_data$Patient_id %in% sample_subtypes$patient)                            # 228 samples (26-08-2024)
samples_data <- samples_data[samples_data$Patient_id %in% sample_subtypes$patient, ] # Keep only annotated patient samples
table(sample_subtypes$`HCC subtypes`[sample_subtypes$patient %in%                    # View HCC subtypes
                                       samples_data$Patient_id])

# 26-08-2024
# Cirrhotomimetic hepatocellular carcinoma      Clear cell hepatocellular carcinoma           Fibrolamellar carcinoma
# 3                                             10                                            4

# Lymphocyte rich hepatocellular carcinoma      Myxoid hepatocellular carcinoma               NA
# 3                                             1                                             13

# No specific subtype                           Sarcomatoid hepatocellular carcinoma          Scirrhous hepatocellular carcinoma 
# 135                                           2                                             8

# Steatohepatitic 
# 10 

samples_data$Subtype <- sapply(samples_data$Patient_id,function(x)   # Assign each sample subtype as metadata column
  sample_subtypes$`HCC subtypes`[sample_subtypes$patient==x])

# Check for multiple samples in patients
duplicates <- samples_data |>
  add_count(Patient_id, Subtype) |>
  filter(n > 1) |>
  distinct()
table(duplicates$Subtype)             # There are 80 duplicates (26-08-2024)

# Cirrhotomimetic hepatocellular carcinoma  Clear cell hepatocellular carcinoma  Fibrolamellar carcinoma  
# 2                                         4                                    2

# NA                                        No specific subtype 
# 8                                         62

#-----Clinical data-----
# Extract clinical data
cli_data <- GDCquery_clinic("TCGA-LIHC","clinical")

#Keep specific clinical data
cli_data <- cli_data[, c("bcr_patient_barcode", "gender", "race", "vital_status", "ajcc_pathologic_t",
                         "ajcc_pathologic_stage", "primary_diagnosis", "prior_malignancy",
                         "ajcc_pathologic_m", "ajcc_pathologic_n",
                         "treatments_pharmaceutical_treatment_or_therapy",
                         "treatments_radiation_treatment_or_therapy")]

# Combine clinical data to expression sample data
samples_data <- cbind(samples_data,t(sapply(samples_data$Patient_id,function(x) 
  cli_data[cli_data$bcr_patient_barcode==x,2:ncol(cli_data)])))

# Basic data description
table(samples_data$Subtype)
table(as.character(samples_data$gender))
table(as.character(samples_data$race))
table(as.character(samples_data$vital_status))
table(as.character(samples_data$primary_diagnosis))
table(as.character(samples_data$ajcc_pathologic_t))
table(as.character(samples_data$ajcc_pathologic_stage))
table(as.character(samples_data$prior_malignancy))
table(as.character(samples_data$treatments_pharmaceutical_treatment_or_therapy))
table(as.character(samples_data$treatments_radiation_treatment_or_therapy))

#-----Make object into a tsv table for later pre-processing-----
samples_data <- apply(samples_data,2,as.character)
write.table(samples_data,"samples_data.tsv",sep='\t', row.names = F, quote = F)
