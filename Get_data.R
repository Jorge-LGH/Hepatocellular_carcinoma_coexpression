# Load libraries
library(SummarizedExperiment)  # Version: 1.34.0
library(TCGAbiolinks)          # Version: 2.32.0
library(VennDiagram)           # Version: 1.7.3

#-----Sample IDs for expression data-----
exp_data <- GDCquery(project = "TCGA-LIHC",                        # Liver hepatocellular carcinoma
                     data.category = "Transcriptome Profiling",    # Self explanatory 
                     data.type = "Gene Expression Quantification", # Self explanatory
                     workflow.type="STAR - Counts")                # Read counts

exp_data <- getResults(exp_data)         # Get results from query (must be a GDCquery object)
exp_ids <- substr(exp_data$cases,1,19)   # Extract all sample ids. Remove institute code
length(exp_ids)                          # There are 424 samples  (06-08-2024)

#-----Features of interest-----
# Make dataframe
exp_samples <- data.frame(cbind(exp_ids,                  # Column with the whole identifier
                                substr(exp_ids, 1, 12)))  # Column with patient identifier
colnames(exp_samples)[2] <- "patient_id"                  # Change second column name

# Get HCC variants and metadata
sample_subtypes <- TCGAquery_subtype(tumor = "LIHC")                              # Table with hepatocellular carcinoma variants
sum(exp_samples$patient_id %in% sample_subtypes$patient)                          # 240 samples (06-08-2024)
exp_samples <- exp_samples[exp_samples$patient_id %in% sample_subtypes$patient, ] # Keep only annotated patient samples

# Duplicated patient ids
sum(exp_samples$patient_id %in% sample_subtypes$patient)          # There are 47 duplicated patient ids (06-08-2024)
exp_samples <- exp_samples[!duplicated(exp_samples$patient_id), ] # Removed duplicated patient ids

# Remove samples with HCC variant set as "NA"
for(ids in exp_samples$patient_id){
  if(sample_subtypes[sample_subtypes$patient == ids,]$`HCC subtypes` != "NA"){
    exp_samples <- exp_samples[!exp_samples$exp_ids == ids, ]
  }
}

# Assign HCC variant type to the dataframe




