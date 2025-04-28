####### Analysis #####

library(tidyverse)
library(readxl)
##### Load data ####
Data = read.table("D:/Lab/Cardiac Muscle/Studies/GSE197670/GSE197670_normalized_batch_corrected_beta_EPIC.txt")
Data = as.data.frame(beta_values_filtered)

pheno = read.csv("D:/Lab/Cardiac Muscle/Studies/GSE197670/GSE197670_pheno_EPIC.csv")

Entropy = read.table("D:/Lab/Cardiac Muscle/Studies/GSE197670/GSE197670_entropy.txt")

design = read_csv("D:/Lab/Cardiac Muscle/Studies/GSE197670/Design/GSE197670_design_matrix.csv")
design = as.matrix(design)

fit_B = readRDS("D:/Lab/Cardiac Muscle/Studies/GSE197670/Fits/GSE197670_fit_B_results.rds")
fit_M = readRDS("D:/Lab/Cardiac Muscle/Studies/GSE197670/Fits/GSE197670_fit_M_results.rds")
M = readRDS("D:/Lab/Cardiac Muscle/Studies/GSE197670/M.rds")

ResidualsMatrix_B = read.table("D:/Lab/Cardiac Muscle/Studies/GSE197670/Residuals/GSE197670_B_res.txt")
ResidualsMatrix_M = read.table("D:/Lab/Cardiac Muscle/Studies/GSE197670/Residuals/GSE197670_M_res.txt")


##### Limma #####
library(limma)
library(minfi)
library(dplyr)
library(broom) # for rownames_to_column
library(readr) # for write_csv

B = Data
M = logit2(B)
View(head(M))

# Ensure the base directory exists
base_dir <- "D:/Lab/Cardiac Muscle/Studies/GSE197670_EPIC"
if (!dir.exists(base_dir)) {
  dir.create(base_dir, recursive = TRUE)
}

# Save M as an .rds file
saveRDS(M, file = file.path(base_dir, "M.rds"))

# Create a model matrix for linear regression
design <- model.matrix(~ Age*Sex + HF.Etiology, data = pheno)

View(design)


# Ensure the Design directory exists
design_dir <- file.path(base_dir, "Design")
if (!dir.exists(design_dir)) {
  dir.create(design_dir, recursive = TRUE)
}

# Save the design matrix to a CSV file in the Design folder
write.csv(design, file = file.path(design_dir, "GSE197670_design_matrix.csv"), row.names = FALSE)


# Ensure the Fits directory exists
fits_dir <- file.path(base_dir, "Fits")
if (!dir.exists(fits_dir)) {
  dir.create(fits_dir, recursive = TRUE)
}

# Fit linear models
fit_B <- lmFit(B, design)
fit_B <- eBayes(fit_B)
saveRDS(fit_B, file = file.path(fits_dir, "GSE197670_fit_B_results.rds"))

fit_M <- lmFit(M, design)
fit_M <- eBayes(fit_M)
saveRDS(fit_M, file = file.path(fits_dir, "GSE197670_fit_M_results.rds"))


##### TopTables #####
library(tibble)

# Extracting coefficient names excluding '(Intercept)'
coefficients <- colnames(design)[!colnames(design) %in% "(Intercept)"]

# Function to generate and save TopTable for a specified coefficient
generate_top_table <- function(fit, coef_name, data_type, study_name) {
  # Clean the coefficient name for the folder and file names
  safe_coef_name <- gsub("[:/]", "_", coef_name)
  
  # Determine the folder name based on the cleaned coefficient name
  folder_name <- safe_coef_name
  
  # Create folder if it doesn't exist
  dir.create(file.path(base_dir, "TopTables", folder_name), showWarnings = FALSE, recursive = TRUE)
  
  # Generate TopTable for the specified coefficient
  if (data_type == "B") {
    CpGs <- topTable(fit, coef = coef_name, adjust = "BH", number = nrow(B), p.value = 1)
  } else if (data_type == "M") {
    CpGs <- topTable(fit, coef = coef_name, adjust = "BH", number = nrow(M), p.value = 1)
  }
  
  # Add standard error to the CpGs dataframe
  SE <- fit$stdev.unscaled * fit$sigma
  CpGs$SE <- SE[rownames(CpGs), coef_name]
  
  CpGs_fs <- rownames_to_column(CpGs)
  
  # Save the TopTable to the corresponding folder with study name in the filename
  write_csv(CpGs_fs, file.path(base_dir, "TopTables", folder_name, paste0("GSE197670_TopTable_", data_type, "_", safe_coef_name, ".csv")))
}

# List of study names
study_name <- "GSE197670"
# Loop through the coefficients and generate TopTables for both B and M
for (coef in coefficients) {
  generate_top_table(fit_B, coef, "B", study_name)
  generate_top_table(fit_M, coef, "M", study_name)
}


##### Residuals #####
# Define the directory name for residuals
residuals_dir <- file.path(base_dir, "Residuals")

# Create the Residuals folder if it doesn't exist
if (!dir.exists(residuals_dir)) {
  dir.create(residuals_dir, recursive = TRUE)
}

# Calculate residuals for both models
ResidualsMatrix_M <- residuals(fit_M, M)
ResidualsMatrix_B <- residuals(fit_B, B)

# Load the data.table package
library(data.table)

# Convert the matrix to a data.table explicitly and then write it
write.table(signif(ResidualsMatrix_M, digits = 4), 
            file = file.path(residuals_dir, "GSE197670_M_res.txt"))

# Convert the matrix to a data.table explicitly and then write it
write.table(signif(ResidualsMatrix_B, digits = 4), 
            file = file.path(residuals_dir, "GSE197670_B_res.txt"))



##### Entropy calculations #####
M_adj <- ResidualsMatrix_M + rowMeans(M)
B_adj <- ilogit2(M_adj)
plotDensities(B_adj, legend = F)

# Calculate entropy 
calculate_entropy <- function(B) {
  return(sum(B * log2(B) + (1 - B) * log2(1 - B)) / (length(B) * log2(1 / 2)))
}

# Beta matrix # 
Entropy <- apply(B_adj, 2, calculate_entropy)
Entropy <- as.data.frame(Entropy)
Entropy$ID <- rownames(Entropy)

# Create the Entropy folder if it doesn't exist
entropy_dir <- file.path(base_dir, "Entropy")
if (!dir.exists(entropy_dir)) {
  dir.create(entropy_dir, recursive = TRUE)
}

# Save the data frame to a CSV file
write.table(Entropy, file = file.path(entropy_dir, "GSE197670_entropy.txt"), sep = ",", row.names = FALSE, col.names = TRUE)

# Merge the entropy data with the pheno data
pheno <- left_join(pheno, Entropy, by = c("ID"))

# Save the updated pheno data frame
write.csv(pheno, file = file.path(base_dir, "GSE197670_pheno_with_entropy_EPIC.csv"), row.names = FALSE)
library(ggplot2)

# Save entropy plot
tiff(file.path(entropy_dir, "GSE197670_Entropy_EPIC.tiff"),
     width = 5,
     height = 3,
     units = 'in',
     res = 200)
ggplot(data = pheno, aes(x = Age, y = Entropy)) +
  geom_jitter(colour = "deepskyblue4") +
  geom_smooth(method = "lm", colour = "black") +
  labs(title = "GSE197670_EPIC Entropy",  # Updated title here
       y = "Entropy") +  # Updated to "y" from "ylab"
  theme_minimal()
dev.off()


fit <- lm(Entropy ~ Age, pheno)
summary(fit)

pval <- summary(fit)$coefficients[, 4][2]
pval 
#0.9638352 

effect <- summary(fit)$coefficients[, 1][2]
effect
#1.75039e-06
stderr <- summary(fit)$coefficients[, 2][2]
stderr
#3.792053e-05



###### Create tables for METAL analysis ######


##### VMP analysis #####
# Load necessary libraries
library(limma)
library(tibble)
library(dplyr)
library(ggplot2)
library(minfi)

design2 <- model.matrix(~ Age, pheno)
resid = ResidualsMatrix_M

shapirotest <- apply(as.matrix(resid), 1, shapiro.test)
pvals <- sapply(shapirotest, `[[`, 2)
CpGs_to_keep <- names(pvals[pvals > 1e-5])
resid <- resid[CpGs_to_keep, ]

resid_B <- ilogit2(resid)

sigma2 <- rowSums(resid^2) / ncol(resid)
w <- as.matrix((resid^2) - sigma2)
fit1 <- lmFit(w, design2)
fit2 <- eBayes(fit1)

sigma2_B <- rowSums(resid_B^2) / ncol(resid_B)
w_B <- as.matrix((resid_B^2) - sigma2_B)
fit1_B <- lmFit(w_B, design2)
fit2_B <- eBayes(fit1_B)


coef <- "Age"
results <- topTable(fit2, coef = coef, number = Inf, p.value = 1)
results_B <- topTable(fit2_B, coef = coef, number = Inf, p.value = 1)
results$logFC <- results_B[rownames(results), "logFC"]
SE <- fit2_B$sigma * fit2_B$stdev.unscaled
results$SE <- SE[rownames(results), coef]

fitted_sqrd <- rowSums(fitted(fit2)^2)
bp <- ncol(resid) * fitted_sqrd / rowSums(w^2)
df <- fit2$rank - 1
chisq_pval <- sapply(bp, pchisq, df = df, lower.tail = FALSE)

results_BP <- tibble(
  CPG = rownames(results),
  ALLELE1 = rep(1, nrow(results)),
  ALLELE2 = rep(2, nrow(results)),
  TESTSTAT = bp,
  N = ncol(w),
  COEF = results$logFC,
  SE = results$SE,
  PVALUE = chisq_pval
)

tiff(file.path("D:/Lab/Cardiac Muscle/All TopTables/Plot", paste0("GSE197670_EPIC_pvalhist_VMPs.tiff")),
     width = 5, height = 3, units = 'in', res = 200)
ggplot(results_BP, aes(x = PVALUE)) +
  geom_histogram(color = "darkblue", fill = "lightblue", bins = 20) +
  labs(title = paste("P-value Distribution for GSE197670_EPIC"),
       x = "p-value") +
  theme_minimal()
dev.off()

# Save results table
write.table(results_BP,
            file = file.path("D:/Lab/Cardiac Muscle/All TopTables/VMP_METAL_TopTables", paste0("GSE197670_EPIC.tbl")),
            quote = FALSE, row.names = FALSE, col.names = TRUE, sep = "\t")


##### Entropy calculations #####
setwd("D:/Lab/Cardiac Muscle/Studies/GSE197670_EPIC")

library(metafor)
library(tidyverse)
library(minfi)
library(limma)

#repeat for ALL age-related CpGS vs non-age-related CpGs

B <- data.table::fread("D:/Lab/Cardiac Muscle/Studies/GSE197670_EPIC/GSE197670_normalized_batch_corrected_beta_EPIC.txt")
B <- as.data.frame(B)
rownames(B) <- B$V1
B <- B %>% select(-V1)

pheno <- read.csv("D:/Lab/Cardiac Muscle/Studies/GSE197670_EPIC/GSE197670_pheno_with_entropy_EPIC.csv")

DMPs <- read.delim("D:/Lab/Cardiac Muscle/DMP_METAL_files/DMPs_005.txt")
#VMPs <- no, since none showed significance 
All_cps <- DMPs$CpG

B <- as.data.frame(B)
B_all <- B[rownames(B) %in% All_cps,]
M_all <- logit2(B_all)

#run limma and adjust for confounders EXCLUDING age
design = model.matrix(~Sex + HF.Etiology,
                      pheno)

fit1 <- lmFit(M_all,
              design)
fit2 <- eBayes(fit1)

resid_M <- residuals(fit2, M_all)

#add the mean beta value for each cpg to the residuals to get a beta matrix 

M_adj <- resid_M + rowMeans(M_all)
B_adj <- ilogit2(M_adj)
plotDensities(B_adj, legend = F)

#calcuate entropy 
calculate_entropy <- function(B){
  return(sum(B*log2(B) + (1-B)*log2(1-B))/(length(B)*log2(1/2)))
}

Entropy_all <- apply(B_adj, 2, calculate_entropy)
Entropy_all <- as.data.frame(Entropy_all)
Entropy_all$ID <- rownames(Entropy_all)
pheno$ID <- as.character(pheno$ID)
pheno <- left_join(pheno, Entropy_all, by = "ID")

library(ggplot2)

# Create the directory if it doesn't exist
dir.create("D:/Lab/Cardiac Muscle/Entropy figures/All-age related CPGs", recursive = TRUE, showWarnings = FALSE)

tryCatch({
  # Save the first plot to the correct folder
  tiff(filename = "D:/Lab/Cardiac Muscle/Entropy figures/All-age related CPGs/GSE197670_EPIC_Entropy_all_age_related_CpGs.tiff", 
       width = 10, height = 7, units = "in", res = 300)
  
  plot <- ggplot(data = pheno, aes(x = Age, y = Entropy_all)) +
    geom_jitter(colour = "deepskyblue4") +
    geom_smooth(method = "lm", colour = "black") +
    labs(title = "GSE197670_EPIC Entropy all age-related CpGs",
         y = "Entropy",  
         x = "Age") +    
    theme_minimal()
  
  print(plot)
}, finally = {
  if (!is.null(dev.list())) dev.off()
})



#calculate entropy on cpgs that do not change with age 

B_none <- B[!rownames(B) %in% All_cps,]
M_none <- logit2(B_none)


fit1 <- lmFit(M_none,
              design)
fit2 <- eBayes(fit1)

resid_M <- residuals(fit2, M_none)

#add the mean beta value for each cpg to the residuals to get a beta matrix 

M_adj <- resid_M + rowMeans(M_none)
B_adj <- ilogit2(M_adj)
plotDensities(B_adj, legend = F)

#calcuate entropy 
calculate_entropy <- function(B){
  return(sum(B*log2(B) + (1-B)*log2(1-B))/(length(B)*log2(1/2)))
}

Entropy_none <- apply(B_adj, 2, calculate_entropy)
Entropy_none <- as.data.frame(Entropy_none)
Entropy_none$ID <- rownames(Entropy_none)
pheno <- left_join(pheno, Entropy_none, by = "ID")

library(ggplot2)

# Create the directory if it doesn't exist
dir.create("D:/Lab/Cardiac Muscle/Entropy figures/Non-age-related CPGs", recursive = TRUE, showWarnings = FALSE)

tryCatch({
  tiff(filename = "D:/Lab/Cardiac Muscle/Entropy figures/Non-age-related CPGs/GSE197670_EPIC_Entropy_non_age_related_CpGs.tiff", 
       width = 10, height = 7, units = "in", res = 300)
  
  plot <- ggplot(data = pheno, aes(x = Age, y = Entropy_none)) +
    geom_jitter(colour = "deepskyblue4") +
    geom_smooth(method = "lm", colour = "black") +
    labs(title = "GSE197670_EPIC Entropy all non-age-related CpGs",
         y = "Entropy", 
         x = "Age") +
    theme_minimal()
  
  print(plot)
}, finally = {
  if (!is.null(dev.list())) dev.off()
})


#summary statistics all 
fit <- lm(Entropy_all ~ Age, pheno)
summary(fit)
pval = summary(fit)$coefficients[,4][2]
pval
#p-value: 0.01175214       
effect = summary(fit)$coefficients[,1][2]
effect
#effect: 0.0009970942      
stderr = summary(fit)$coefficients[,2][2]
stderr
#stderr: 0.0003643249

#summary statistics non-age-related cpgs
fit <- lm(Entropy_none ~ Age, pheno)
summary(fit)
pval = summary(fit)$coefficients[,4][2]
pval
#p-value: 0.0005125831       
effect = summary(fit)$coefficients[,1][2]
effect
#effect: -5.696658e-05     
stderr = summary(fit)$coefficients[,2][2]
stderr
#stderr: 1.635944e-05 



write.csv(pheno,
          "D:/Lab/Cardiac Muscle/Studies/GSE197670_EPIC/GSE197670_EPIC_pheno_final.csv")


# Load necessary libraries
library(readxl)
library(openxlsx)

# Read the existing Excel file
df <- read_excel("D:/Lab/Cardiac Muscle/Thesis_Tables.xlsx")

# Define the summary statistics for all-age related CpGs
fit_all <- lm(Entropy_all ~ Age, pheno)
all_pval <- summary(fit_all)$coefficients[,4][2]
all_effect <- summary(fit_all)$coefficients[,1][2]
all_stderr <- summary(fit_all)$coefficients[,2][2]

# Define the summary statistics for non-age related CpGs
fit_none <- lm(Entropy_none ~ Age, pheno)
non_pval <- summary(fit_none)$coefficients[,4][2]
non_effect <- summary(fit_none)$coefficients[,1][2]
non_stderr <- summary(fit_none)$coefficients[,2][2]

# Find the row corresponding to "GSE197670_EPIC"
row_index <- which(df$`Dataset ID` == "GSE197670_450K")

# Add the new columns with the summary statistics
df[row_index, "all_pval"] <- all_pval
df[row_index, "all_effect"] <- all_effect
df[row_index, "all_stderr"] <- all_stderr
df[row_index, "non_pval"] <- non_pval
df[row_index, "non_effect"] <- non_effect
df[row_index, "non_stderr"] <- non_stderr

write.xlsx(df, file = "D:/Lab/Cardiac Muscle/Thesis_Tables.xlsx", rowNames = FALSE)




#####
pheno = read.csv("D:/Lab/Cardiac Muscle/Studies/GSE197670_EPIC/GSE197670_EPIC_pheno_final.csv")
pheno_entropy <- pheno %>% select(ID, Age, Entropy_DMPsonly, Entropy_none, Entropy_full)

pheno_entropy <- pheno_entropy %>% pivot_longer(!c(Age,ID),
                                                names_to = "Entropy type",
                                                values_to = "Entropy value")

names(pheno_entropy)[names(pheno_entropy) == "Entropy type"] <- "Condition"

library(viridis)
setwd("D:/Lab/Cardiac Muscle/Studies/GSE197670_EPIC")

tiff('GSE197670_EPIC + alpha Entropy DMPs, all, none.tiff',
     width = 10,
     height = 5,
     units = 'in',
     res = 400)

# Rename the column to "Entropy"
colnames(pheno_entropy)[colnames(pheno_entropy) == "Entropy value"] <- "Entropy"

p <- ggplot(pheno_entropy, aes(x = Age, y = Entropy)) + 
  geom_jitter(aes(colour = Condition, alpha = 0.8)) +  # Removed alpha from aes()
  geom_smooth(aes(colour = Condition), method = "lm") + 
  scale_color_manual(values = c("purple", "navy", "coral")) +
  scale_fill_viridis(discrete = TRUE) +
  theme_minimal() +
  ggtitle("GSE197670_EPIC Entropy breakdown")

p
dev.off()
