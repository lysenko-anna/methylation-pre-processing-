#Get the GEO dataset
library(tidyverse)
library(dplyr)
library(data.table)
directory = "D:/Lab/Atlas/Pancreas/GSE49149"

#Check sex
setwd(directory)
# Read only the selected columns
pheno <- read_csv("D:/Lab/Atlas/Pancreas/GSE49149/Final/GSE49149_pheno_corrected.csv")
print(head(pheno))

##### Modifications made to phenotypes #####
library(dplyr)

# Separate the 'Title' column into 'Sentrix_ID' and 'Sentrix_Position'
pheno <- pheno %>%
  mutate(Sentrix_ID = sub("_.*", "", Title),  # Extract the part before the underscore
         Sentrix_Position = Position) %>%
  select(-Position)  # Remove the original Position column

# Check the result
print(head(pheno))

library(dplyr)

# Update the 'Cellularity' column based on the condition
pheno <- pheno %>%
  mutate(Cellularity = ifelse(Cellularity == Age, "adjacent", Cellularity))

# Check the result
print(head(pheno))

write.csv(pheno, 
          "GSE49149_phenotypes.csv")
######

##### Unzip downloaded fron GEO files #####
# Define the file path and destination directory
file_path <- "D:/Lab/Atlas/Pancreas/GSE49149/GSE49149_RAW.tar"
dest_dir <- "D:/Lab/Atlas/Pancreas/GSE49149/Idats"

# Extract the contents of the .tar file
untar(file_path, exdir = dest_dir)

idat_path <- ("D:/Lab/Atlas/Pancreas/GSE49149/Idats")
#memory.limit(size=1000000000)
# List all .idat.gz files
gz_files <- list.files(idat_path, pattern = "\\.idat\\.gz$", full.names = TRUE)

# Decompress each .idat.gz file
for (file in gz_files) {
  gunzip(file, remove = TRUE)  # Set remove = TRUE if you want to delete the .gz file after decompression
}

###### Check sex #####
#Check Sex
library(minfi)
library(R.utils)

RGSet <- read.metharray.exp("D:/Lab/Atlas/Pancreas/GSE49149/Idats", force=T)
pData(RGSet)  # Shows phenotype or sample data
annotation_data <- getAnnotation(RGSet)
head(annotation_data)



library(IlluminaHumanMethylation450kmanifest)
MSet <- preprocessRaw(RGSet)
RSet <- ratioConvert(MSet, what = "both", keepCN = TRUE)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
GRset <- mapToGenome(RSet)
predictedSex <- getSex(GRset, cutoff = -2)$predictedSex


# Standardize predictedSex values to match 'Male' and 'Female'
test <- pheno %>%
  mutate(predictedSex = case_when(
    predictedSex == "M" ~ "Male",
    predictedSex == "F" ~ "Female",
    TRUE ~ predictedSex  # Keep other values unchanged if any
  ))

# Now compare Sex and predictedSex columns to check for mismatches
mismatched_samples <- test$Sex != test$predictedSex

# Get the mismatched samples
mismatched_samples_list <- test[mismatched_samples, ]

# Print the mismatched samples
print(mismatched_samples_list)

# Optionally, count the number of mismatches
num_mismatches <- sum(mismatched_samples)
cat(num_mismatches, "mismatched samples\n") #GSM1435245


#Filter probes
library(ChAMP)
library(R.utils)

# Define the file path
gz_file <- "D:/Lab/Atlas/Pancreas/GSE49149/Idats/GPL13534_HumanMethylation450_15017482_v.1.1.csv.gz"
output_file <- sub("\\.gz$", "", gz_file)  # Remove ".gz" extension to get the output file name

# Unzip the file
gunzip(gz_file, destname = output_file)

# Check if the file is unzipped
print(paste("File unzipped to:", output_file))


library(readr)

# Load phenotype data
pheno <- read_csv("D:/Lab/Atlas/Pancreas/GSE49149/GSE49149_phenotypes.csv")

# Extract sample IDs from the 'Sample' column (adjust the column name if needed)
phenotype_samples <- pheno$Sample

# Set the path to your IDAT files directory
idat_directory <- "D:/Lab/Atlas/Pancreas/GSE49149/Idats"

# List all IDAT files in the directory (assuming the file extension is .idat)
idat_files <- list.files(idat_directory, pattern = "\\.idat$", full.names = TRUE)

# Extract sample IDs from both green and red IDAT filenames (either Grn or Red)
idat_samples <- gsub("_.*_R01C01_(Grn|Red)\\.idat$", "", basename(idat_files))

# Select IDAT files that match the samples with phenotype data
idat_files_selected <- idat_files[idat_samples %in% phenotype_samples]

# View the selected IDAT files
print(idat_files_selected)

"D:/Lab/Atlas/Pancreas/GSE49149/Idats/GSM1194354_6042324127_R01C01_Grn.idat.gz"



myLoad <- champ.load(directory = "D:/Lab/Atlas/Pancreas/GSE49149/Idats",
                     filterXY=TRUE, #Remove sex chr
                     filterSNPs=TRUE,
                     arraytype="450K")
colnames(myLoad$beta) <- pheno$Sample

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

#Filter samples with missmatched Sex
# Specify the known mismatched samples directly
mismatched_samples <- c("GSM1435245")  # Replace with your known mismatched samples

# Filter the phenotype table to exclude mismatched samples
pheno_filtered <- pheno[!pheno$Sample %in% mismatched_samples, ]

# Filter the beta values table to exclude mismatched samples
# Ensure the column names in beta_values match Sample_ID in pheno
beta_values_filtered <- subbeta[, !colnames(subbeta) %in% mismatched_samples]

# Check the dimensions of the filtered tables
dim(pheno_filtered)
dim(beta_values_filtered)

# Assign the filtered data back to the myLoad object
myLoad$beta <- beta_values_filtered
myLoad$pd <- pheno_filtered


#Produce quality control graphs to look at the data
# Save the plot as a TIFF file
# Run champ.QC
champ.QC(beta = myLoad$beta,
         pheno = myLoad$pd$Sex,
         dendrogram = FALSE)

# Close the graphical device after saving the plot
dev.off()


#Normalization of Type I and Type II probes
myNorm <- champ.norm(beta = myLoad$beta,
                     arraytype = "450K",
                     resultsDir = "./CHAMP_QCimages/")  # Automatically saves in the specified folder



library(data.table)


# Use fwrite to save the normalized data
fwrite(myNorm, 
       file = "GSE49149_norm.txt", 
       quote = FALSE, 
       row.names = TRUE, 
       col.names = TRUE, 
       sep = "\t")


# Read the file saved with fwrite
myNorm2 <- fread("D:/Lab/Atlas/Pancreas/GSE49149/GSE49149_norm.txt", sep = "\t")

# Check the data
head(myNorm)

library(ChAMP)
library(dplyr)

# Check the sample names in both datasets
colnames(myNorm)  # Sample names in beta data (column names of myNorm)
head(pheno_filtered$Sample)  # Sample names in phenotype data (Sample column in pheno_filtered)

# Ensure that pheno_filtered only contains the samples present in myNorm and in the same order
pheno_filtered <- pheno_filtered %>%
  filter(Sample %in% colnames(myNorm)) %>%
  arrange(match(Sample, colnames(myNorm)))  # Reorder pheno_filtered to match the order in myNorm

# Convert pheno_filtered to a data.frame if it's a tibble
pheno_filtered <- as.data.frame(pheno_filtered)

# Now set rownames
rownames(pheno_filtered) <- pheno_filtered$Sample

# Proceed with the SVD analysis
champ.SVD(beta = as.data.frame(myNorm),
          pd = pheno_filtered,
          resultsDir = "./CHAMP_SVDimages/")


# Save the plot to a file
tiff("SVD_plot.tiff", width = 800, height = 600)
champ.SVD(beta = as.data.frame(myNorm),
          pd = pheno_filtered,
          resultsDir = "./CHAMP_SVDimages/")
dev.off()  # Close the device to save the file



library(sva)

#Run ComBat to correct batch effects
M <- logit2(myNorm)

myCombat=ComBat(dat=as.matrix(M), #it outputs an M-value matrix adjusted for batch
                batch=pheno_filtered$Sentrix_ID,
                mod=NULL)

#Convert back to beta-values after batch correction to run SVD again
myCombat=ilogit2(myCombat)

#Run SVD again
champ.SVD(beta=as.data.frame(myCombat),
          pd=pheno_filtered,
          resultsDir="./CHAMP_SVDimages/batch_corrected/")

M <- logit2(myCombat)

myCombat=ComBat(dat=as.matrix(M), #it outputs an M-value matrix adjusted for batch
                batch=pheno_filtered$Sentrix_Position,
                mod=NULL)

myCombat=ilogit2(myCombat)


fwrite(myCombat,
            file="GSE49149_normalized_batch_corrected_beta.txt",
            quote=FALSE,
            row.names=TRUE,
            col.names=TRUE,
            sep='\t')

champ.SVD(beta=as.data.frame(myCombat),
          pd=as.data.frame(pheno_filtered),
          resultsDir="./CHAMP_SVDimages/batch_position_corrected/")

#Produce quality control graphs to look at the data
champ.QC(beta = myCombat,
         pheno = pheno_filtered$Age,
         dendrogram = FALSE)


##### Remove the cancerous cells ####
pheno <- read_csv("D:/Lab/Atlas/Pancreas/GSE49149/Final/GSE49149_pheno_corrected.csv")
Data <- fread("D:/Lab/Atlas/Pancreas/GSE49149/Final/GSE49149_normalized_batch_corrected_beta.txt")
Data = as.data.frame(Data)
rownames(Data) = Data$V1
Data$V1 <- NULL


tumor_rows <- which(pheno$Cell_type == "Tumor")
tumor_samples <- pheno$Sample[tumor_rows]
pheno <- pheno[-tumor_rows, ]
Data <- Data[, !(colnames(Data) %in% tumor_samples)]
pheno <- pheno[, !(colnames(pheno) %in% c("Tissue", "Cellularity", "Cell_type"))]
#Produce new Champ images 

library(ChAMP)
champ.QC(beta = as.matrix(Data),
         pheno = pheno$Sex,
         dendrogram = FALSE,
         resultsDir = "./CHAMP_QCimages/Cancer-corrected")
champ.SVD(beta = Data,
          pd = as.data.frame(pheno),
          resultsDir="./CHAMP_SVDimages/Cancer-corrected")


fwrite(pheno, "D:/Lab/Atlas/Pancreas/GSE49149/Final/GSE49149_pheno_corrected.csv")
fwrite(Data, "D:/Lab/Atlas/Pancreas/GSE49149/Final/GSE49149_normalized_batch_corrected_beta.txt", 
       row.names = T,
       sep = "\t",
       quote = FALSE)

##### Check new pheno #####
summarize_age_sex <- function(pheno) {
  # Ensure Age is numeric
  pheno$Age <- as.numeric(pheno$Age)
  
  # Calculate Age metrics
  age_mean <- mean(pheno$Age, na.rm = TRUE)       # Mean Age
  age_sd <- sd(pheno$Age, na.rm = TRUE)           # Standard Deviation
  age_min <- min(pheno$Age, na.rm = TRUE)         # Minimum Age
  age_max <- max(pheno$Age, na.rm = TRUE)         # Maximum Age
  
  # Calculate % Male
  male_count <- sum(pheno$Sex == "Male", na.rm = TRUE)
  total_count <- sum(!is.na(pheno$Sex))           # Total valid Sex entries
  percent_male <- (male_count / total_count) * 100
  
  # Print formatted summary
  cat(sprintf("Age (mean ± SD): %.2f ± %.2f\n", age_mean, age_sd))
  cat(sprintf("Age range: %d - %d\n", age_min, age_max))
  cat(sprintf("%% Male: %.2f%%\n", percent_male))
}

# Usage Example
summarize_age_sex(pheno)



