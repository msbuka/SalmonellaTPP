######################
### Load libraries ###
######################
library(tidyverse)
library(vsn)
library(limma)
library(fdrtool)

##################
#### Load data ###
##################

### load sample mapping file (Adjust file path)
sample_mapping <- read.csv("Sample_mapping.csv")

### load results files (Adjust file path)
results_files <- c("S1975_protein.tsv", 
                   "S1976_protein.tsv", 
                   "S1977_protein.tsv", 
                   "S1978_protein.tsv", 
                   "S1979_protein.tsv", 
                   "S2032_protein.tsv",  
                   "S2033_protein.tsv", 
                   "S2034_protein.tsv", 
                   "S2035_protein.tsv", 
                   "S2036_protein.tsv", 
                   "S2037_protein.tsv",  
                   "S2038_protein.tsv", 
                   "S2039_protein.tsv", 
                   "S2040_protein.tsv", 
                   "S2041_protein.tsv", 
                   "S2074_protein.tsv", 
                   "S2075_protein.tsv", 
                   "S2076_protein.tsv", 
                   "S2077_protein.tsv", 
                   "S2078_protein.tsv", 
                   "S2079_protein.tsv",
                   "S2080_protein.tsv", 
                   "S2081_protein.tsv", 
                   "S2082_protein.tsv", 
                   "S2083_protein.tsv")

### combine data files into one data frame and filter
prot_melt <- bind_rows(lapply(results_files, function(MS_file){ prot_tab <- read_delim(MS_file, delim = "\t") %>%  
    filter(!grepl("contam_sp", `Protein ID`), ### remove contaminants
           `Unique Peptides` > 1, ### remove all proteins with less than 2 unique peptides
            grepl('Mus musculus',Organism)) ### keep only proteins assigned to host (Mus musculus)
 
  
##################################################
### Variance stabilization normalization (vsn) ### 
##################################################

### prepare for VSN 
### create a matrix of the intensity columns from the prot_tab data frame
signal_sum_mat <- as.matrix(prot_tab %>% dplyr::select(matches("channel_"))) 
###  filter out any infinite or NaN values
signal_sum_mat[is.infinite(log2(signal_sum_mat)) | is.nan(log2(signal_sum_mat))] <- NA 

### The following code performs normalization of signal intensity data using the VSN method.
### The first eight columns of the signal_sum_mat matrix contain the signal intensities of the temperatures 37°, 44, 49.8, 55.5, 62°C
### vsn2 is a function that fits a VSN model to the input matrix, and the resulting model is stored in vsn_fit_2.
  vsn_fit_1 <- vsn2(signal_sum_mat[,1:8])
  
### This line of code applies the predict function to the vsn_fit_1 model and the first eight columns of signal_sum_mat, 
### which generates a VSN-normalized matrix of the same dimensions as the input matrix. 
### The resulting matrix is then converted to a data frame and stored in vsn_norm_mat_1.
  vsn_norm_mat_1 <- as.data.frame(predict(vsn_fit_1, signal_sum_mat[,1:8]))
  
### The first column of vsn_norm_mat_1 contains the normalization offset for each sample (Here timepoint 0), 
### which is subtracted from all the other columns to obtain the final VSN-normalized data. 
### This line of code subtracts the first column from all the other columns in vsn_norm_mat_1. 
  vsn_norm_mat_1 <- vsn_norm_mat_1 - vsn_norm_mat_1[,1]
  
### The next eight columns of signal_sum_mat contain the remaining signal intensities of the temperatures 40.4, 46.9, 52.9, 58.6, 66.3
  vsn_fit_2 <- vsn2(signal_sum_mat[,9:16])
  vsn_norm_mat_2 <- as.data.frame(predict(vsn_fit_2, signal_sum_mat[,9:16]))
  vsn_norm_mat_2 <- vsn_norm_mat_2 - vsn_norm_mat_2[,1]
  
### The VSN-normalized data for both sets of temperature are combined into a single data frame
  vsn_norm_mat <- cbind(vsn_norm_mat_1,vsn_norm_mat_2)

### This line of code adds a column to vsn_norm_mat that contains protein IDs, which are obtained from the prot_tab data frame
  vsn_norm_mat$protein_id <- prot_tab$`Protein ID`
 
### Merge the normalized data with the remaining data in prot_tab
### first prepare the normalized data for merging
  ### reshape the data frame so that the values in each column are stacked on top of each other, creating a long format.
  vsn_norm_mat %>% ### input data set 
    tbl_df %>%  ### converts the input data set into a tibble data frame
    gather(key,value,-protein_id) %>% ### reshape
                                      ### key = name of the new column that will contain the original column names, 
                                      ### value = the name of the new column that will contain the original values. 
                                      ### The -protein_id argument specifies that the protein_id column should not be included in the new columns.

  ### create a new column called MS_sample_number, which is derived from the MS_file column. Here we extract the experiment number which,in this case is made up by a string of an S and 4 digits. 
    mutate(MS_sample_number = sub(".*/(S\\d{4})_.*\\.tsv", "\\1", MS_file)) %>% ##ä
  ### rename the value column to log2fc    
    dplyr::rename(log2fc = value) %>%
    
### Now combine the current data frame with prot_tab using a left join (keeping all rows from vsn_norm_mat and matching rows from prot_tab) 
### The gene_name, protein_id, and description columns are joined with the current data frame based on the matching values in the protein_id column.
    left_join(prot_tab %>% dplyr::select(gene_name = Gene,protein_id = `Protein ID`,description = `Protein Description`))
  
}))

### The output of this code is a data frame named "prot_melt"


###############################################################
### Abundance scores, stability scores and stability weight ###
###############################################################

### 1) Prepare the data frame

prot_merged <- prot_melt %>% ### Create a new data frame called prot_merged on the basis of prot_melt
  left_join(sample_mapping %>% ### Join the sample_mapping to prot_melt using a left join based on the specified key
              mutate(key = paste0('channel_',TMT_label)) %>% ### create a new column called key containing the string 'channel_' and the values from the TMT_label column. 
              dplyr::select(Time_point,Replicate,MS_sample_number,key,Temperature,Sample), ###  select specific columns from sample_mapping
            by = c('key' = 'key','MS_sample_number' = 'MS_sample_number')) %>% ### merge by "key"
  dplyr::select(-key) %>% ### Remove the key column 
 
  ### make new columns in which we count how often a replicate is there.
  ### For abundance we use only the temperature 37C - so there should be 3 replicates for each timepoint at 37C. 
   group_by(protein_id, Time_point, Sample) %>% ### Group the data frame by the protein_id, Time_point, and Sample 
  mutate(rep_abun = sum(Temperature == 37), ### Create two new columns in the data frame called rep_abun and rep_stab. 
                                            ### The rep_abun column is calculated by summing the number of times Temperature equals 37 within each group
         rep_stab = n()) %>%                ### The rep_stab column is calculated by counting the number of rows within each group.
  ungroup() %>%  ### Remove the grouping 
  filter(Time_point != 0) ### Filter the data frame to remove any rows where the Time_point column is equal to 0.

#### 2) generate a table of "abundance scores" 
abun_scores_raw <- prot_merged %>%  
  filter(rep_abun >= 2, ### only include rows where the rep_abun column is greater than or equal to 2...
         Temperature %in% sort(unique(sample_mapping$Temperature))[1:2]) %>% ### ... and where the Temperature column is one of the first two unique values found in the Temperature column of the sample_mapping (37 and 40.4°C)
  group_by(protein_id,Time_point,Replicate,Sample) %>%  ### group the filtered data frame by four columns: protein_id, Time_point, Replicate, and Sample. 
  mutate(abun_score = mean(log2fc, na.rm = TRUE)) %>%  ### create a column called abun_score, which contains the mean value of the log2fc column
  ungroup() %>%  ### remove the grouping 
  dplyr::select(protein_id,Time_point,Replicate,Sample,abun_score) %>% ### Select only the five columns that we want to include in the final output.
  distinct() ### Removes any duplicate rows, so that we only have one row per unique combination of protein_id, Time_point, Replicate, and Sample.

### 3) Bring the table of "abundance scores" into the right format
abun_scores <- abun_scores_raw %>% 
  #### make a new column called "condition" that contains a string composed of "t" followed by the value of Time_point, followed by the text "_sample", followed by the value of Sample, followed by the text "_rep", followed by the value of Replicate. 
  mutate(condition = paste0('t',Time_point,'_sample',Sample,'_rep',Replicate)) %>% 
  dplyr::select(protein_id,condition,abun_score) %>% ### This selects only the three columns that we want to include in the final output
  spread(condition,abun_score) ### Transform the data frame from a "long" format to a "wide" format -> Specifically, it creates a new column for each unique value of condition, and fills in the values of abun_score in the corresponding cells.

abun_scores <- as.data.frame(abun_scores) ### convert the data type of the abun_scores to a data frame if it wasn't already. 
  ### Set the protein_id to row names and remove the column protein_id
rownames(abun_scores) <- abun_scores$protein_id 
abun_scores$protein_id <- NULL

#### 4) generate a table of "stability scores" 

stab_scores_raw <- prot_merged %>% 
  filter(rep_abun >= 2, ### only include rows where the rep_abun column is greater than or equal to 2...
         rep_stab >= 10) %>% ### and where the rep_stab column is greater than or equal to 10...
  left_join(abun_scores_raw) %>%  ### Join the abun_scores_raw data to stab_scores_raw 
  group_by(protein_id,Time_point,Replicate,Sample) %>% ### group the filtered data frame by four columns: protein_id, Time_point, Replicate, and Sample. 
  mutate(stab_score = sum(log2fc-abun_score, na.rm = TRUE), ### Add two new columns: stab_score calculates the sum of the difference between log2fc and abun_score, with any missing values removed (na.rm = TRUE)
         stab_score_weight = n()) %>%                       ### stab_score_weight calculates the number of rows in each group (n()).
  ungroup() %>%  ###  remove the grouping
  dplyr::select(protein_id,Time_point,Replicate,Sample,stab_score,stab_score_weight) %>%  ### Select only the columns protein_id, Time_point, Replicate, Sample, stab_score, and stab_score_weight.
  distinct() %>%  ### removes any duplicate rows from the data set.
  filter(stab_score_weight > 4) ### Remove rows where stab_score_weight is less than or equal to 4.

### 5) Bring the table of "stability scores" into the right format

stab_scores <- stab_scores_raw %>% 
  mutate(condition = paste0('t',Time_point,'_sample',Sample,'_rep',Replicate)) %>% ### new column called condition = string composed of "t" followed by the value of Time_point, followed by the text "_sample", followed by the value of Sample, followed by the text "_rep", followed by the value of Replicate. 
  dplyr::select(protein_id,condition,stab_score) %>% ### Select only the columns protein_id, condition, and stab_score.
  spread(condition,stab_score) ### Transform the data frame from a "long" format to a "wide" format

stab_scores <- as.data.frame(stab_scores) ### converts the resulting object to a data frame.
  ### Set the protein_id to row names and remove the column protein_id
rownames(stab_scores) <- stab_scores$protein_id
stab_scores$protein_id <- NULL

#### 6) generate a table of "stability weights" from stab_scores_raw and bring it in the right format

stab_weights <- stab_scores_raw %>% 
  mutate(condition = paste0('t',Time_point,'_sample',Sample,'_rep',Replicate)) %>% ## new column called condition = string composed of "t" followed by the value of Time_point, followed by the text "_sample", followed by the value of Sample, followed by the text "_rep", followed by the value of Replicate. 
  dplyr::select(protein_id,condition,stab_score_weight) %>% ### Select only the columns protein_id, condition, and stab_score_weight.
  spread(condition,stab_score_weight) ### Transform the data frame from a "long" format to a "wide" format 

stab_weights <- as.data.frame(stab_weights)### converts the resulting object to a data frame.
  ### Set the protein_id to row names and remove the column protein_id
rownames(stab_weights) <- stab_weights$protein_id
stab_weights$protein_id <- NULL


#######################################
### Statistical analysis with limma ###
#######################################

### For each time point and sample combination, the corresponding subset of abundance (abun_scores_temp), stability (stab_scores_temp), and stability weights (stab_weights_temp) data is extracted from the full dataset. 
limma_results <- bind_rows(lapply(paste0(rep(paste0('t',unique(sample_mapping$Time_point)[-1]),each = length(unique(sample_mapping$Sample))),'_sample',unique(sample_mapping$Sample)), function(time_point) {
  abun_scores_temp <- abun_scores[,grep(time_point,colnames(abun_scores))]
  abun_scores_temp <- abun_scores_temp[rowSums(!is.na(abun_scores_temp)) > 1,]
  stab_scores_temp <- stab_scores[,grep(time_point,colnames(stab_scores))]
  stab_weights_temp <- stab_weights[,grep(time_point,colnames(stab_weights))]
  stab_scores_temp <- stab_scores_temp[rowSums(!is.na(stab_scores_temp)) > 1,]
  stab_weights_temp <- stab_weights_temp[rownames(stab_scores_temp),]

  ### Rows with missing values are removed, and metadata information is generated using the remaining data. 
  
  metadata <- data.frame(ID = colnames(abun_scores_temp),
                         time_point = gsub("_.+","",colnames(abun_scores_temp)),
                         sample = gsub("_.+","",gsub(".+_sample","",colnames(abun_scores_temp))),
                         rep = gsub(".+_","",colnames(abun_scores_temp)))
  rownames(metadata) <- metadata$ID
  
  abundance.dataE <- ExpressionSet(assayData = as.matrix(abun_scores_temp),
                                   phenoData = AnnotatedDataFrame(metadata))
  stability.dataE <- ExpressionSet(assayData = as.matrix(stab_scores_temp),
                                   phenoData = AnnotatedDataFrame(metadata))
  
### linear models are fit using lmFit with eBayes empirical Bayes adjustment.  
  abun_comparison <- eBayes(lmFit(abundance.dataE))
  stab_comparision <- eBayes(lmFit(stability.dataE,weights = stab_weights_temp^2))

### For each time point and sample combination, two tables of results (abun_res and stab_res) are created using topTable, which includes the log-fold change (logFC), p-value (P.Value), adjusted p-value (adj.P.Val), and protein ID for each protein.    
  abun_res <- limma::topTable(abun_comparison, sort.by = "t",  
                              coef = 1, number = Inf)
  stab_res <- limma::topTable(stab_comparision, sort.by = "t",  
                              coef = 1, number = Inf)
  
  abun_res$protein_id <- rownames(abun_res)
  abun_res$comparison <- 'abundance'
  abun_res$time_point <- unique(metadata$time_point)
  abun_res$sample <- unique(metadata$sample)
  stab_res$protein_id <- rownames(stab_res)
  stab_res$comparison <- 'stability'
  stab_res$time_point <- unique(metadata$time_point)
  stab_res$sample <- unique(metadata$sample)

### The fdrtool function is then used to calculate the false discovery rate (FDR) for each comparison, and the results are appended to the respective tables. 
  abun_res$FDR_tool <- fdrtool(abun_res$t,verbose = FALSE, plot = FALSE)$qval
  stab_res$FDR_tool <- fdrtool(stab_res$t,verbose = FALSE, plot = FALSE)$qval

### The two tables are combined using rbind and formatted
  rbind(abun_res,stab_res) %>% 
    tbl_df() %>% 
    dplyr::select(protein_id,time_point,sample,comparison,logFC,P.Value,adj.P.Val,FDR_tool) %>% 
    dplyr::rename(p_value = P.Value,
                  FDR_limma = adj.P.Val)
}))

### scaled_limma_results is then created by grouping the results by comparison (abundance or stability) and calculating the z-score of the log-fold change
### Proteins are considered significant if the absolute value of their z-score is greater than the 97.5th percentile of a standard normal distribution (qnorm(0.975)) and the FDR is less than 0.05. 
scaled_limma_results <- limma_results %>% 
  group_by(comparison) %>% 
  mutate(z_score = scale(logFC)) %>% 
  ungroup() %>% 
  mutate(hit = abs(z_score) > qnorm(0.975) & FDR_tool < 0.05)

########################
### Save data to csv ###
########################

### combine all data into one data frame and save as file (Adjust file path)
write.csv(scaled_limma_results %>% 
          left_join(prot_merged %>% 
          dplyr::select(gene_name,protein_id,description) %>% distinct()) %>% 
          dplyr::select(gene_name,protein_id,description,time_point,sample,comparison,logFC,z_score,p_value,FDR_limma,FDR_tool,hit),
         'results_TPP_Salmonella_infection_timecourse.csv', row.names = FALSE)