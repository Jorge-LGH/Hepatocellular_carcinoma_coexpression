# Load libraries
library(SummarizedExperiment)  # Version: 1.34.0
library(TCGAbiolinks)          # Version: 2.32.0
library(VennDiagram)           # Version: 1.7.3
library(dplyr)                 # Version: 1.1.4

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
samples_exp <- data.frame(cbind(exp_ids,                 # Column with the whole identifier
                                substr(exp_ids, 1, 12))) # Column with patient identifier
colnames(samples_exp) <- c("experiment_id","patient_id") # Change column names

# Get HCC variants and metadata
sample_subtypes <- TCGAquery_subtype(tumor = "LIHC")                              # Table with hepatocellular carcinoma variants
sum(samples_exp$patient_id %in% sample_subtypes$patient)                          # 240 samples (06-08-2024)
samples_exp <- samples_exp[samples_exp$patient_id %in% sample_subtypes$patient, ] # Keep only annotated patient samples

# Assign HCC subtypes 
subtypes <- c()                                             # Empty vector to extract HCC subtypes                                            
for(ids in samples_exp$patient_id){
  subtypes <- c(subtypes, sample_subtypes[which(            # Extract the subtypes
    sample_subtypes$patient == ids),]["HCC subtypes"][[1]])
}
samples_exp$subtype <- subtypes                             # Adding the subtype column

# Check for multiple samples for single patients
duplicates <- samples_exp |>
  add_count(patient_id, subtype) |>
  filter(n > 1) |>
  distinct()
table(duplicates$subtype)             # There are 94 duplicates (08-08-2024)

#Cirrhotomimetic hepatocellular carcinoma  |  Clear cell hepatocellular carcinoma  |  Fibrolamellar carcinoma  
#2                                           4                                        2
#NA  |  No specific subtype  |  Scirrhous hepatocellular carcinoma  |  Steatohepatitic 
#8      74                      2                                      2

#-----Clinical data-----
# Extract clinical data
cli_data <- GDCquery_clinic("TCGA-LIHC","clinical")

#Keep specific clinical data
cli_data <- select(cli_data, c("bcr_patient_barcode", "gender", "race", "vital_status", "ajcc_pathologic_t",
                               "ajcc_pathologic_stage", "primary_diagnosis", "prior_malignancy",
                               "ajcc_staging_system_edition", "ajcc_pathologic_m", "ajcc_pathologic_n",
                               "treatments_pharmaceutical_treatment_or_therapy",
                               "treatments_radiation_treatment_or_therapy"))

# Combine clinical data to expression sample data
samples_exp <- cbind(samples_exp,t(sapply(samples_exp$patient_id,function(x) 
  cli_data[cli_data$bcr_patient_barcode==x,-1])))


