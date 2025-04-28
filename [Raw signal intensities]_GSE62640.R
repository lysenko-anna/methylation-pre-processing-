# Read the entire file to inspect column names
library(data.table)
library(tidyr)
library(dplyr)

full_data <- fread("D:/Lab/Atlas/Pancreas/GSE62640/GSE62640_signal_intensities.txt", header = TRUE, sep = "\t")

# Manually select methylated and unmethylated columns based on position
methylated_cols <- names(full_data)[seq(2, ncol(full_data), by = 3)]  # Methylated columns (1st, 4th, 7th, ...)
unmethylated_cols <- names(full_data)[seq(3, ncol(full_data), by = 3)]  # Unmethylated columns (2nd, 5th, 8th, ...)


# Assuming you've already loaded full_data

# Extract the ID_ref column
id_ref <- full_data$ID_REF

# Create a dataframe for methylated columns (ID_ref + methylated columns)
methylated_data <- full_data[, c("ID_REF", methylated_cols), with = FALSE]

# Create a dataframe for unmethylated columns (ID_ref + unmethylated columns)
unmethylated_data <- full_data[, c("ID_REF", unmethylated_cols), with = FALSE]

# Remove "methylated" from column names in methylated_data
colnames(methylated_data) <- gsub("methylated", "", colnames(methylated_data))

# Remove "unmethylated" from column names in unmethylated_data
colnames(unmethylated_data) <- gsub("unmethylated", "", colnames(unmethylated_data))

# Preview the modified column names
head(methylated_data)
head(unmethylated_data)


library(minfi)

# Reshape methylated_long and unmethylated_long data

# Remove the ID_REF column and make sure they align by CpG
methylated_matrix <- as.matrix(methylated_data[, -1])  # Remove the first column (ID_REF)
unmethylated_matrix <- as.matrix(unmethylated_data[, -1])  # Remove the first column (ID_REF)

# Ensure that row names (ID_REF) match across both datasets
rownames(methylated_matrix) <- methylated_data$ID_REF
rownames(unmethylated_matrix) <- unmethylated_data$ID_REF

rownames(methylated_matrix) <- sub("^cg", "", rownames(methylated_matrix))
rownames(unmethylated_matrix) <- sub("^cg", "", rownames(unmethylated_matrix))

# Create the RGChannelSet
rgSet <- RGChannelSet(Green = methylated_matrix, Red = unmethylated_matrix)
# Install the annotation package for 450K:
BiocManager::install("IlluminaHumanMethylation450kanno.ilmn12.hg19")

# Also, make sure minfiData is installed:
BiocManager::install("minfiData")



RGset <- preprocessRaw(rgSet)
library(manifest)

MSet <- MSet <- MSet <- preprocessRaw(rgSet) # Potentially something went wrong here 

MSet_manual <- MethylSet(Meth = methylated_matrix, Unmeth = unmethylated_matrix, 
                         annotation = "IlluminaHumanMethylation450k", 
                         preprocessMethod = "raw")


RSet <- ratioConvert(MSet_manual, what = "both", keepCN = TRUE)

##### Saving the B-values #####
# Define the RSet file and output directory
output_file <- "D:/Lab/Atlas/Pancreas/GSE62640/GSE62640_BetaValues.csv"

# Step 1: Extract Beta values from RSet
BetaValues <- assay(RSet, "Beta")  # Extract Beta values

# Step 2: Save Beta values to a CSV file
write.csv(BetaValues, file = output_file, row.names = TRUE)

# Confirmation message
cat("Beta values have been successfully saved to:", output_file, "\n")



MethylSet(Meth = methylated_matrix, Unmeth = unmethylated_matrix,
          annotation = "", preprocessMethod = "")

Mset <- MethylSet(Meth = methylated_matrix, Unmeth = unmethylated_matrix,
                  annotation = c(array = "IlluminaHumanMethylation450k", annotation = "ilmn12.hg19"))

# Create MethylSet object manually
Mset2 <- MethylSet(Meth = methylated_matrix, Unmeth = unmethylated_matrix,
                  annotation = "IlluminaHumanMethylation450kanno.ilmn12.hg19",
                  preprocessMethod = "manual")
Mset
Mset2


# Map to genome
GMset <- mapToGenome(Mset)
GMset

GMset2 <- mapToGenome(Mset2) #this one doesn't work since Error in getAnnotationObject(object) :cannot load annotation package IlluminaHumanMethylation450kanno.ilmn12.hg19anno

predictedSex <- getSex(GMset, cutoff = -2)$predictedSex


# Get the genomic locations of the probes
probeLocations <- rowRanges(GMset)

# Filter for probes on the X and Y chromosomes
sexChromosomes <- probeLocations[seqnames(probeLocations) %in% c("chrX", "chrY")]

# Extract the methylation data for these probes
sexData <- assay(GMset)[names(sexChromosomes), ]

# Check variability of the data
variability <- apply(sexData, 1, var)
print(variability)

# If variability looks good, run getSex function
predictedSex <- getSex(GMset, cutoff = -2)$predictedSex
print(predictedSex)


sum(is.na(getBeta(GMset)))  # Number of missing values
beta_matrix_clean <- getBeta(GMset)
beta_matrix_clean <- beta_matrix_clean[complete.cases(beta_matrix_clean), ]  # Remove rows with NAs

pca <- prcomp(t(beta_matrix_clean))
summary(pca)  # View the PCA results
plot(pca)  # Visualize the PCA results



# Replace NA values in GMset$beta with the smallest possible value
GMset$beta[is.na(GMset$beta)] <- min(GMset$beta, na.rm = TRUE) - 1e-6




GRset <- mapToGenome(RSet,mergeManifest = FALSE)
predictedSex <- getSex(GRset, cutoff = -2)$predictedSex
print(predictedSex)


# Check what these variables hold
print(methylated_cols)
print(unmethylated_cols)

# Now use these column names in the `readGEORawFile` function
result <- readGEORawFile(
  filename = "D:/Lab/Atlas/Pancreas/GSE62640/GSE62640_signal_intensities.txt",
  sep = "\t",
  Uname = unmethylated_cols,  # Columns for unmethylated values
  Mname = methylated_cols,    # Columns for methylated values
  row.names = 1,
  pData = pheno_full,
  array = "450K"
)











# Load necessary libraries
library(tidyverse)
library(dplyr)
library(data.table)
library(minfi)
library(ChAMP)
library(sva)
library(IlluminaHumanMethylation450kmanifest)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
library(R.utils)

# Define your directory
directory <- "D:/Lab/Atlas/Pancreas/GSE62640"

# Set working directory and read phenotype data
setwd(directory)
pheno_full <- read_csv("D:/Lab/Atlas/Pancreas/GSE62640/GSE62640_phenotypes_full.csv")
print(colnames(pheno_full))

# Select the specified columns from the pheno dataframe
pheno <- pheno[, c("Sample_title", 
                              "Sample_geo_accession", 
                              "Sample_description", 
                              "Donor ID", 
                              "Sex", 
                              "Tissue", 
                              "Bisulfite Conversion Batch", 
                              "Illumina Chip Number")]

# Print the selected columns
print(selected_columns)



##### Already done, skip for now ######
# Preprocess phenotype data by separating and cleaning columns
pheno <- pheno %>%
  mutate(Sentrix_ID = sub("_.*", "", Title),  # Extract part before the underscore
         Sentrix_Position = Position) %>%
  select(-Position) %>%  # Remove original Position column
  mutate(Cellularity = ifelse(Cellularity == Age, "adjacent", Cellularity))  # Update 'Cellularity' column

# Save cleaned phenotype data
write.csv(pheno, "GSE49149_phenotypes.csv")

# Extract and unzip the raw data
file_path <- "D:/Lab/Atlas/Pancreas/GSE62640/GSE62640_RAW.tar"
dest_dir <- "D:/Lab/Atlas/Pancreas/GSE62640/Idats"
untar(file_path, exdir = dest_dir)

# List and decompress .idat.gz files
idat_files <- list.files(dest_dir, pattern = "\\.idat\\.gz$", full.names = TRUE)
lapply(idat_files, gunzip, remove = TRUE)  # Decompress files

# Load data into a methylation object and preprocess


# Read the data from the file (treating ! as a comment)
data = read.table("D:/Lab/Atlas/Pancreas/GSE62640/GSE62640_signal_intensities.txt", header = TRUE, comment.char = "!", sep = "\t", quote = "\"")
print(head(data))

# Read the data, excluding comment lines starting with '!'
data2 <- read.table("D:/Lab/Atlas/Pancreas/GSE62640/GSE62640_processed_data.txt", header = TRUE, comment.char = "!", sep = "\t", quote = "\"")

library(dplyr)

# Define the ilogit2 function
ilogit2 <- function(M) {
  2 ^ M / (1 + 2 ^ M)
}

# Assuming data2 is your preprocessed data frame with M values
# Print the first few rows of data2
head(data2)

# Select the columns that are not detection p-values
data_filtered <- data2 %>%
  select(-contains("detection.p.values"))

# Convert M values to B values using the ilogit2 function
data_b_values <- data_filtered %>%
  mutate(across(-ID_REF, ilogit2))

# Print the first few rows of the data with B values
print(head(data_b_values))

# Save the data with B values to a new CSV file
fwrite(data_b_values, file = "D:/Lab/Atlas/Pancreas/GSE62640/GSE62640_data_b_values.txt", row.names = FALSE)
library(dplyr)

# Save data_b_values as a CSV file
write.csv(data_b_values, file = "data_b_values.csv", row.names = TRUE)
print("data_b_values has been saved to data_b_values.csv")

# Transpose data_b_values so that Sample IDs become a new column and CpGs become column names
data_b_values_transposed <- as.data.frame(t(data_b_values))
colnames(data_b_values_transposed) <- data_b_values_transposed[1,]
data_b_values_transposed <- data_b_values_transposed[-1,]

data_b_values_transposed$Sample_description <- rownames(data_b_values_transposed)

# Read the phenotype table
phenotype_table <- read.csv("path/to/your/phenotype_table.csv", stringsAsFactors = FALSE)

# Merge the transposed data with the phenotype table
merged_data <- merge(pheno_full, data_b_values_transposed, by = "Sample_description")

# Print the first few rows of the merged data
print(colnames(merged_data))




print(head(data))
print(head(data2))

# Extract the header row (line 75 in your case, assuming it's the header with GSM IDs)
header_line <- data[75, ]  # Adjust this if necessary (this might depend on how your file is structured)

# Convert this line to a dataframe
header_df <- as.data.frame(t(header_line))

# Optionally, clean the column names (remove any extra quotes or special characters)
colnames(header_df) <- gsub("\"", "", colnames(header_df))

# Print the resulting dataframe
print(header_df)


# Save the data back to the same directory using fwrite
fwrite(pheno_data, "D:/Lab/Atlas/Pancreas/GSE62640/GSE62640_series_matrix_cleaned.txt")


##### Check for sex####

# Convert Beta matrix to GenomicRatioSet
GRSet <- makeGenomicRatioSetFromMatrix(as.matrix(merged_data))
# Check for CpG rownames
head(rownames(GRSet))

GRSet <- makeGenomicRatioSetFromMatrix(as.matrix(data_b_values_filtered),
                                       pData = pheno_full)

# Print the first few rows of the phenotype data
head(pData(GRSet))



####
library(minfi)
library(minfiData)
library(dplyr)

# Load Beta matrix and phenotype table
beta_matrix <- read.csv("path/to/data_b_values.csv", row.names = 1)
phenotype_table <- read.csv("path/to/your/phenotype_table.csv", stringsAsFactors = FALSE)

# Load the Illumina 450k annotation data
data(IlluminaHumanMethylation450kmanifest)
official_probe_names <- getManifestInfo(IlluminaHumanMethylation450kmanifest, "locusNames")

# Check and filter row names
data_b_values_filtered <- data_b_values[data_b_values$ID_REF %in% official_probe_names, ]

# Assuming the first column of data_b_values_filtered contains the CpG IDs
# Set the first column as row names and remove it from the dataframe
rownames(data_b_values_filtered) <- data_b_values_filtered[, 1]
data_b_values_filtered <- data_b_values_filtered[, -1]

# Print the first few rows to verify
head(data_b_values_filtered)

glimpse(pheno_full)



library(dplyr)
library(chammi)  # Alternative package for sex prediction

# Load Beta matrix and phenotype table
beta_matrix <- read.csv("path/to/data_b_values.csv", row.names = 1)
pheno_full <- read.csv("path/to/your/phenotype_table.csv", stringsAsFactors = FALSE)

# Ensure CpG IDs are correctly set as row names in the beta matrix
beta_matrix <- as.data.frame(t(beta_matrix))
colnames(beta_matrix) <- gsub("\\.0$", "", colnames(beta_matrix))  # Adjust column names if needed

# Match Beta matrix sample IDs with phenotype table
colnames(beta_matrix) <- pheno_full$Sample_description




# Convert the filtered Beta matrix to GenomicRatioSet
GRSet <- makeGenomicRatioSetFromMatrix(as.matrix(data_b_values_filtered))

# Predict sex
predictedSex <- getSex(GRSet, cutoff = -2)$predictedSex

# Add predicted sex to the phenotype table
phenotype_table <- phenotype_table %>%
  mutate(predictedSex = predictedSex,
         Sex = ifelse(Sex == "Male", "M", "F"))  # Standardize Sex column

# Check for discrepancies
discrepancy <- which(phenotype_table$Sex != phenotype_table$predictedSex)
if (length(discrepancy) == 0) {
  print("All good! Sex matches for all samples.")
} else {
  print("Discrepancies found in the following samples:")
  print(phenotype_table[discrepancy, ])
}

# Print the first few rows of the updated phenotype table
head(phenotype_table)

###### Read in Unmethylated and Methylated signals from a GEO raw file ######




dest_dir <- "D:/Lab/Atlas/Pancreas/GSE49149/Idats"
RGSet <- read.metharray.exp(dest_dir, targets = "D:/Lab/Atlas/Pancreas/GSE62640/GSE62640_signal_intensities.txt") #changable to B-matrix 
MSet <- preprocessRaw(RGSet)
RSet <- ratioConvert(MSet, what = "both", keepCN = TRUE)

# Map to genome
GRset <- mapToGenome(RSet)
predictedSex <- getSex(GRset, cutoff = -2)$predictedSex

# Standardize predictedSex and check for mismatches with 'Sex' column in phenotype data
pheno <- pheno %>%
  mutate(predictedSex = recode(predictedSex, "M" = "Male", "F" = "Female"))
mismatched_samples <- pheno$Sex != pheno$predictedSex
print(pheno[mismatched_samples, ])





myLoad <- champ.load(directory = "D:/Lab/Atlas/Pancreas/GSE49149/Idats",
                     filterXY=TRUE, #Remove sex chr
                     filterSNPs=TRUE,
                     arraytype="450K")

library(ChAMP)

pd = read.csv("D:/Lab/Atlas/Pancreas/GSE62640/GSE62640_phenotypes_full.csv")
head(pd)  # Inspect the top rows of the phenotype data
# Rename the column in the phenotype data (pd)
colnames(pd)[colnames(pd) == "Sample_description"] <- "Sample_Name"
# Verify the change
head(pd)
all(pd$Sample_Name %in% colnames(BetaValues))  # False
# Samples in pd but not in BetaValues
mismatched_in_pd <- setdiff(pd$Sample_Name, colnames(BetaValues))
print(mismatched_in_pd)

# Samples in BetaValues but not in pd
mismatched_in_beta <- setdiff(colnames(BetaValues), pd$Sample_Name)
print(mismatched_in_beta)
head(pd$Sample_Name)
head(colnames(BetaValues))
colnames(BetaValues) <- trimws(colnames(BetaValues))
all(pd$Sample_Name %in% colnames(BetaValues))  # Should return TRUE
myLoad <- champ.filter(beta = BetaValues, 
                       pd = pd, 
                       filterDetP = FALSE, 
                       filterBeads = FALSE, 
                       filterXY = TRUE, 
                       filterSNPs = TRUE, 
                       arraytype = "450K")
colnames(myLoad$beta) <- pd$Sample_Name
library(readxl)
sheets <- excel_sheets("D:/Lab/Annotations/Cross reactive and SNP probes from Pidsley et al.xlsx")
subbeta <- myLoad$beta
library(tidyverse)
for (s in sheets)
{
  qcprobes <- read_excel("D:/Lab/Annotations/Cross reactive and SNP probes from Pidsley et al.xlsx",
                         sheet=s)%>%
    pull(ProbeID)
  subbeta <- subbeta[setdiff(rownames(subbeta),qcprobes),]
}


# Set the output directory
output_dir <- "D:/Lab/Atlas/Pancreas/GSE62640"

# Specify the file name for the PDF output
pdf_filename <- file.path(output_dir, "QC_Graphs.pdf")

# Open a PDF graphical device
pdf(pdf_filename, width = 10, height = 8)

# Run champ.QC to produce QC graphs
champ.QC(beta = myLoad$beta,
         pheno = myLoad$pd$Sex,
         dendrogram = FALSE)

# Close the graphical device to save the plot
dev.off()

boxplot(myLoad$beta)
summary(myLoad$beta)
na_rows <- which(rowSums(is.na(myLoad$beta)) > 0)
myLoad$beta[na_rows, ]
na_cols <- which(colSums(is.na(myLoad$beta)) > 0)
myLoad$beta[, na_cols]




# Normalize Type I and Type II probes
# Normalize methylation data with specified method
myNorm <- champ.norm(beta = myLoad$beta, 
                     arraytype = "450K")

fwrite(myNorm, "D:/Lab/Atlas/Pancreas/GSE62640/GSE62640_norm.txt", quote = FALSE, row.names = TRUE, col.names = TRUE, sep = "\t")

# Check consistency of sample names between phenotype and normalized data
myNorm = as.data.frame(myNorm)
pheno_filtered <- pheno_filtered %>%
  filter(Sample %in% colnames(myNorm)) %>%
  arrange(match(Sample, colnames(myNorm)))  # Ensure order matches

# Perform Singular Value Decomposition (SVD) analysis
tiff(filename = "./CHAMP_SVDimages/ChampSVDPlot.tiff", width = 2000, height = 1500, units = "px", res = 300)
champ.SVD(beta = myNorm, pd = pd)
dev.off()

S
# Correct batch effects using ComBat (a statistical method to remove batch-related variation)
# Convert beta-values to M-values for ComBat
M <- logit2(myNorm)
myCombat <- ComBat(dat = as.matrix(M), batch = pheno_filtered$Sentrix_ID, mod = NULL)
myCombat <- ilogit2(myCombat)  # Convert back to beta-values

# Run SVD again after batch correction
champ.SVD(beta = as.data.frame(myCombat), pd = pheno_filtered, resultsDir = "./CHAMP_SVDimages/batch_corrected/")

# Convert to M-values again and correct for 'Sentrix_Position'
M <- logit2(myCombat)
myCombat <- ComBat(dat = as.matrix(M), batch = pheno_filtered$Sentrix_Position, mod = NULL)
myCombat <- ilogit2(myCombat)

# Save batch-corrected data and perform another SVD
fwrite(myCombat, "GSE49149_normalized_batch_corrected_beta.txt", quote = FALSE, row.names = TRUE, col.names = TRUE, sep = '\t')
champ.SVD(beta = as.data.frame(myCombat), pd = pheno_filtered, resultsDir = "./CHAMP_SVDimages/batch_position_corrected/")

# Final quality control after batch correction
champ.QC(beta = myCombat, pheno = pheno_filtered$Age, dendrogram = FALSE)
dev.off()  # Close plot device
