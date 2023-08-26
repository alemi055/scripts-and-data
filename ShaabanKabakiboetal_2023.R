
#################################################################################################
#################################################################################################
#################################################################################################

# Audrée Lemieux
# Updated on August 25, 2023

# This script was adapted from Lemieux, A., ..., Dubé, M. & Kaufmann, D. E. Enhanced detection of
# antigen-specific T cells by a multiplexed AIM assay. In preparation. (2023)

#################################################################################################

# Load libraries
library(flowCore)
library(FlowSOM)
library(M3C)
library(pheatmap)
library(ggplot2)
library(reshape2)
library(ggfortify)
library(grid)
library(ggplotify)
library(ggpubr)
library(data.table)
library(Biobase)
library(Rphenograph)
library(gridExtra)
# library(directlabels)
library(ggrepel)
library(scales)


### To modify ###
data_dir <- "/Unsupervised Analysis/data"
today <- "230825"
params <- c("IL2", "IFNg", "TNFa", "CD107a", "IL17a") # Channels of interest

#################################################################################################

                                          #####################
                                          # Declare functions #
                                          #####################

# Function 1 out 2
check_colnames <- function(List){
  # (list) -> bool
  #
  # Input:
  #   - List: list of the datat of each donor
  #
  # Returns a list:
  #   [[1]] True if the column names are the same for all files, or False if not
  #   [[2]] Positions of files/donors whose column names are different
  
  colnames_og <- colnames(List[[1]])
  flag <- TRUE
  n <- length(List)
  i <- 2
  pos <- NULL
  
  while (flag){
    if (i <= n){
      tmp <- colnames(List[[i]])
      if (identical(colnames_og, tmp)){ # If the column names are identical, keep going
        flag <- TRUE
        i <- i + 1
      }else{ # If not, flag the position
        flag <- TRUE
        pos <- c(pos, i)
        i <- i + 1
      }
      
    }else{
      flag <- FALSE
    }
  }
  
  if (length(pos) > 0){
    final_list <- list(FALSE, pos)
  }else{
    final_list <- list(TRUE, NULL)
  }
  return(final_list)
}


# Function 2 out of 2
flowSet_excl_1 <- function(flowSet_AIM){
  # (flowSet) -> flowSet
  #
  # Input:
  #   - flowSet_AIM: flowSet of the AIM data, containing all monkeys
  #
  # Returns a flowSet excluding monkeys who only have 1 gated event
  
  num_gated_events <- fsApply(flowSet_AIM, function(ff) { # Get total number of gated events
    nrow(ff)
  })
  
  pos <- which(num_gated_events[,1] == 1) # Get pos of monkeys who only have 1 gated event
  if (length(pos) > 0){
    new_flowSet <- flowSet_AIM[-pos]
  }else{
    new_flowSet <- flowSet_AIM
  }
  return(new_flowSet)
}

#################################################################################################


                                      ###################
                                      # 1 Load the data #
                                      ###################

PrimaryDirectory <- getwd()
current_dir <- getwd()

# Load CSV files
FileNames <- list.files(paste0(data_dir, "/"), pattern = ".csv")
DataList=list()

for (File in FileNames) {
  # tempdata <- fread(File, check.names = FALSE) # fread was not working for me
  tempdata <- read.csv(paste0(data_dir, "/", File))
  File <- gsub(".csv", "", File)
  DataList[[File]] <- tempdata
}

rm(tempdata, File) # Clean space
AllSampleNames <- names(DataList) # Get the names of all samples

# Clean the names
AllSampleNames <- gsub("  ", " ", AllSampleNames)
AllSampleNames <- gsub(" ", "_", AllSampleNames)
AllSampleNames <- gsub(".fcs", "", AllSampleNames)
AllSampleNames <- gsub("export_2022-04-05_", "", AllSampleNames)
AllSampleNames <- gsub("export_2022-04-12_", "", AllSampleNames)
AllSampleNames <- gsub("export_2022-04-27_", "", AllSampleNames)
AllSampleNames <- gsub("export_2022-05-11_", "", AllSampleNames)
AllSampleNames <- gsub("export_2022-08-25_", "", AllSampleNames)
AllSampleNames <- gsub("export_2022-11-03_", "", AllSampleNames)
AllSampleNames <- gsub("export_2022-11-07_", "", AllSampleNames)
AllSampleNames <- gsub("export_2022-11-18_", "", AllSampleNames)
AllSampleNames <- gsub("export_2022-11-22_", "", AllSampleNames)
names(DataList) <- AllSampleNames

# Make sure that column names are the same for all files
# In our case, this is caused by missing colours (e.g., "Comp.APC.A", "Comp.APC.R700.A", etc.) in the names
# Sample #39 also has two random parameters "Event"
tmp <- check_colnames(DataList)
if (!tmp[[1]]){ # If not, fix the column names
  print(paste0("-- Warning, the column names of donors \'", paste(AllSampleNames[tmp[[2]]], collapse = "\', \'"), "\' are different --"))

  # Remove parameters with "Event"
  for (pos in tmp[[2]]){
    events_pos <- grep("Event", colnames(DataList[[pos]]))
    if (length(events_pos) > 0){
      DataList[[pos]] <- DataList[[pos]][,-events_pos]
    }
  }

  # Since column names are the same (only missing targets), use colnames of 1st data
  colnames_og <- colnames(DataList[[1]])
  for (pos in tmp[[2]]){
    colnames(DataList[[pos]]) <- colnames_og
  }

  # Double-check again
  cytokine_cluster <- check_colnames(DataList)
  if (cytokine_cluster[[1]]){
    print(paste0("-- Yes, the column names are now all the same. --"))
  }else{
    print(paste0("-- No, the column names are not all the same. Go check back. --"))
  }
}

rm(tmp, cytokine_cluster, events_pos) # Clean space

#################################################################################################

                        ##################################################
                        # 2 Do the downsampling before converting to FCS #
                        ##################################################

# Print total number of gated events
k <- sapply(DataList, function(sample) {
  # print(nrow(ff))
  nrow(sample)
})
k <- as.data.frame(k)

# Get histogram
ggplot(k, aes(k)) +
  geom_histogram(col = "black", fill = "grey80") +
  theme_bw() +
  labs(x = "Number of gated events", y = "Frequency") +
  geom_vline(aes(xintercept = median(k)), col = "red", linetype = "dashed", linewidth = 1) +
  annotate(geom="text", label = paste0("median = ", median(k$k)), x = 2500, y = 100, color = "red")

# Downsample to 1000 before conversion to FCS (too many cells otherwise)
new_DataList <- DataList
sampling.ceiling <- 1000

for (i in 1:length(DataList)){
  idx <- sample.int(nrow(DataList[[i]]), min(sampling.ceiling, nrow(DataList[[i]])))
  new_DataList[[i]] <- DataList[[i]][idx,]
}

# Print total number of gated events after downsampling
k2 <- sapply(new_DataList, function(sample) {
  # print(nrow(ff))
  nrow(sample)
})
k2 <- as.data.frame(k2)

#################################################################################################

                          ###########################################
                          # 3 Get the proportion of cytokine+ cells #
                          ###########################################

# Before converting to FCS, for each cell, identify if cell is pos. or neg. for each cytokine
# FlowSOM transforms the data, which means that we need to do this beforehand
YN_cyt <- list()
cutoff_values <- as.data.frame(cbind(798, 1002, 517, 3500, 649))
cutoff_values <- as.data.frame(rbind(cutoff_values, cbind(818, 843, 510, 3500, 750)))
colnames(cutoff_values) <- params
row.names(cutoff_values) <- c("CD4", "CD8")

for (i in 1:length(new_DataList)){
  tmp <- new_DataList[[i]]
  df_cyt <- as.data.frame(cbind(rep(1:nrow(tmp)), rep(1:nrow(tmp)), rep(1:nrow(tmp)), rep(1:nrow(tmp)), rep(1:nrow(tmp))))
  colnames(df_cyt) <- params
  
  for (j in 1:length(params)){
    marker <- params[j]
    pos <- grep(marker, colnames(tmp))
    
    # Is it CD4 or CD8?
    if (grepl("CD4", names(new_DataList)[i])){
      Y_pos <- which(tmp[,pos] >= as.numeric(cutoff_values[1,j]))
      N_pos <- setdiff(1:nrow(tmp), Y_pos)
    }else{
      Y_pos <- which(tmp[,pos] >= as.numeric(cutoff_values[2,j]))
      N_pos <- setdiff(1:nrow(tmp), Y_pos)
    }
    
    # Write "Y" if the cell is positive for the AIM, or "N" if negative
    df_cyt[Y_pos, j] <- "Y"
    df_cyt[N_pos, j] <- "N"
  }
  YN_cyt[[i]] <- df_cyt
}
names(YN_cyt) <- names(new_DataList)
rm(df_cyt, tmp, i, j, marker, N_pos, pos, Y_pos) # Clean space

# Get the proportion of cytokine+ cells
cytokine_prop <- NULL

for (i in 1:length(YN_cyt)){
  tmp <- YN_cyt[[i]]
  tmp_df <- NULL
  
  for (j in 1:5){
    tmp_table <- as.data.frame(table(tmp[,j]))
    tmp_df <- cbind(tmp_df, as.numeric(tmp_table[1,2]), as.numeric(tmp_table[2,2]))
    colnames(tmp_df)[(ncol(tmp_df)-1):ncol(tmp_df)] <- c(paste0("N_", params[j]), paste0("Y_", params[j]))
  }
  
  # Add to main df
  cytokine_prop <- rbind(cytokine_prop, tmp_df)
}
cytokine_prop <- as.data.frame(cytokine_prop)
row.names(cytokine_prop) <- names(YN_cyt)
write.csv(cytokine_prop, file = paste0(today, "_cytokine+_prop.csv"), quote = F) # Save as CSV
rm(tmp, tmp_df, tmp_table, i, j) # Clean space

save.image(paste0(PrimaryDirectory, "/", today, "_unsupervised_ICS.RData")) # Save image

#################################################################################################

                                  ##############################
                                  # 4 Convert CSV to FCS files #
                                  ##############################

# Convert CSV to FCS
outDir <- 'new_FCS'
newdir <- paste(current_dir, outDir, sep = "/")

dir.create(newdir, showWarnings = FALSE)
setwd(newdir)
current_dir <- getwd()

# Iterate on samples and save as FCS
for (i in 1:length(AllSampleNames)){
  data_subset <- rbindlist(as.list(new_DataList[i]))
  cols <- colnames(data_subset)
  # dim(data_subset)
  a <- names(new_DataList)[i]

  metadata <- data.frame(name = colnames(data_subset), desc = paste('column', colnames(data_subset), "from dataset"))

  # Create FCS file metadata - ranges, min, and max settings
  metadata$minRange <- apply(data_subset, 2, min)
  metadata$maxRange <- apply(data_subset, 2, max)

  data_subset.ff <- new("flowFrame", exprs = as.matrix(data_subset), parameters = AnnotatedDataFrame(metadata)) # In order to create a flow frame, data needs to be read as matrix by exprs
  # print(colnames(data_subset.ff))
  write.FCS(data_subset.ff, paste0(a, ".fcs"))
}

rm(data_subset, data_subset.ff, metadata, a, cols, i, outDir, FileNames, colnames_og, pos) # Clean space
save.image(paste0(PrimaryDirectory, "/", today, "_unsupervised_ICS.RData"))

#################################################################################################

                                  ###############################
                                  # 5 Prepare data for analysis #
                                  ###############################

# Get FCS files
fcs_files <- list.files(paste0(newdir, "/"), pattern = ".fcs")

# Read all AIM files as a FlowSET
fs.aim <- read.flowSet(paste0(newdir, "/", fcs_files), ignore.text.offset = T, truncate_max_range = FALSE)

# Find which IDs are associated to the parameters
colnames <- sapply(strsplit(colnames(new_DataList[[1]]), "\\.\\.\\.\\."), function(x){ # Get only markers' names
  x[2]
})

ids <- c()
for (i in params){
  ids <- c(ids, grep(paste0("^", i, "$"), colnames))
}

# Use FlowSOM to integrate the FCS files, apply compensation, scaling, and transformation
fs <- FlowSOM(fs.aim, compensate = FALSE, transform = TRUE, toTransform = ids, scale = TRUE, colsToUse = ids, nClus = 10)

# Print channels corresponding to ids, to validate markers
print(paste(c('Channels = ', paste(colnames(fs[["data"]])[ids],collapse=', ')),collapse=''))

# Clean names
names(fs$metaData) <- fs.aim@phenoData@data[["name"]] # Assign donor's name to metaData

rm(current_dir) # Clean space
save.image(paste0(PrimaryDirectory, "/", today, "_unsupervised_ICS.RData"))

#################################################################################################

                                          ##################
                                          # 6 Run analysis #
                                          ##################

# Extract the transformed MFI values from the flow set
df <- data.frame(fs$data)

# Perform UMAP only on the channels of interest
um <- umap(t(df[,ids]))
save.image(paste0(PrimaryDirectory, "/", today, "_unsupervised_ICS.RData")) # Save image

# Parse sample-wise metadata info
out <- c()
for (j in names(fs$metaData)){
  idx <- fs$metaData[j]
  r <- (idx[[1]][2] - idx[[1]][1]) + 1 # Get number of gated events
  out <- c(out, rep(j, r))
}
out <- gsub(".fcs", "", out)

# Create output data frame with UMAP coordinates
dat <- um$data[,1:2]
colnames(dat) <- c('UMAP_1','UMAP_2')

# # Cluster with Phenograph
# pheno_cluster150 <- Rphenograph(df[,ids], k = 150)
# save.image(paste0(PrimaryDirectory, "/", today, "_unsupervised_ICS.RData")) # Save image
pheno_cluster400 <- Rphenograph(df[,ids], k = 400)
save.image(paste0(PrimaryDirectory, "/", today, "_unsupervised_ICS.RData")) # Save image

pheno_cluster <- pheno_cluster400 # Only keep pheno_cluster whose k = 400
dat$pheno_cluster <- paste0('C', pheno_cluster[[2]]$membership) # Rename each cluster by "C"

rm(fcs_files, j, r, newdir, pheno_cluster150, pheno_cluster400) # Clean space
save.image(paste0(PrimaryDirectory, "/", today, "_unsupervised_ICS.RData"))

# Add sample annotation based on the out variable
dat$everything <- out # For each dot, get the sample
dat$donor <- sapply(dat$everything, function(x){strsplit(x, "_")[[1]][1]}) # For each dot, get the donor
dat$type <- sapply(dat$everything, function(x){strsplit(x, "_")[[1]][2]}) # For each dot, get the type
dat$ratio <- sapply(dat$everything, function(x){strsplit(x, "_")[[1]][4]}) # For each dot, get the ratio
dat$T_cell <- sapply(dat$everything, function(x){strsplit(x, "_")[[1]][5]}) # For each dot, get the type of T cell (CD4 or CD8)

# Add transformed MFI
a <- data.frame(fs$data)[,ids]

# Add transformed MFI to data frame
dat.num <- data.frame(cbind(dat, a))

# Add information on cytokine+ cells 
dat.num <- cbind(dat.num, 1, 1, 1, 1, 1)
colnames(dat.num)[(ncol(dat.num)-4):ncol(dat.num)] <- paste0("YN_", colnames(YN_cyt[[1]]))

for (i in 1:length(YN_cyt)){
  name <- names(YN_cyt)[i]
  pos <- which(dat.num$everything == name)
  dat.num[pos,(ncol(dat.num)-4):ncol(dat.num)] <- YN_cyt[[i]]
}

# Save as CSV: cytokine+ cells per donor in each cluster
dat.num$everything_cluster <- paste0(dat.num$everything, ".", dat.num$pheno_cluster)
count <- 1
cytokine_cluster_donor <- NULL

for (j in 15:19){
  tmp <- as.data.frame.matrix(table(dat.num$everything_cluster, dat.num[,j]))
    
  # Check: is there only one column? If so, that means that the other column is empty.
  if (ncol(tmp) == 1){
    if (colnames(tmp) ==  "N"){
      tmp$Y <- rep(0, nrow(tmp))
    }else{
      tmp$N <- rep(0, nrow(tmp))
    }
  }
  colnames(tmp) <- paste0(colnames(tmp), "_", params[count])
    
  if (count == 1){
    cytokine_cluster_donor <- tmp
    count <- count + 1
  }else{
    cytokine_cluster_donor <- cbind(cytokine_cluster_donor, tmp)
    count <- count + 1
  }
}
write.csv(cytokine_cluster_donor, paste0(today, "_cytokine_cluster_donor.csv"), quote = F)

# Save as CSV (cytokine+ cells in each cluster)
count <- 1
cytokine_cluster <- NULL

for (j in 15:19){
  tmp <- as.data.frame.matrix(table(dat.num$pheno_cluster, dat.num[,j]))
  colnames(tmp) <- paste0(colnames(tmp), "_", params[count])
  
  if (count == 1){
    cytokine_cluster <- tmp
    count <- count + 1
  }else{
    cytokine_cluster <- cbind(cytokine_cluster, tmp)
    count <- count + 1
  }
}
write.csv(cytokine_cluster, paste0(today, "_cytokine_cluster.csv"), quote = F)

rm(a, dat, out, i, idx, name, pos, tmp, data_tmp, j, count, i, tmp) # Clean space
save.image(paste0(PrimaryDirectory, "/", today, "_unsupervised_ICS.RData")) # Save image

#################################################################################################

                                          ###############
                                          # 7 Plot data #
                                          ###############

# Put clusters in alphabetical order
dat.num$pheno_cluster <- factor(dat.num$pheno_cluster, levels = unique(dat.num$pheno_cluster))

# Find median coordinates of each cluster
med_coord <- NULL
clusters <- as.character(unique(dat.num$pheno_cluster))
for (i in clusters){
  tmp <- dat.num[dat.num$pheno_cluster == i,]
  UMAP_1 <- median(tmp$UMAP_1)
  UMAP_2 <- median(tmp$UMAP_2)
  med_coord <- rbind(med_coord, cbind(UMAP_1, UMAP_2, pheno_cluster = i, matrix("", nrow = 1, ncol = ncol(dat.num) - 3)))
}
med_coord <- as.data.frame(med_coord)
colnames(med_coord) <- colnames(dat.num)
class(med_coord$UMAP_1) <- class(med_coord$UMAP_2) <- "numeric"
rm(UMAP_1, UMAP_2) # Clean space

# Plot UMAP (total cells)
umap_total <- ggplot(dat.num, aes(UMAP_1, UMAP_2, col = pheno_cluster)) +
  geom_point(size = 0.5) +
  theme_bw() +
  labs(colour = "Cluster", title = expression(Total~CD4^{"+"}~and~CD8^{"+"}~T~cells)) +
  theme(plot.title = element_text(hjust = 0.5)) + # Center title
  guides(colour = guide_legend(override.aes = list(size = 3))) + # Increase size of dots in legend
  geom_label_repel(data = med_coord, aes(label = pheno_cluster), box.padding = 2, show.legend = F) + # Add labels on map
  scale_y_continuous(breaks = pretty_breaks(), limits = c(min(dat.num$UMAP_2), max(dat.num$UMAP_2))) + # To remove decimals
  scale_x_continuous(limits = c(min(dat.num$UMAP_1), max(dat.num$UMAP_1)))
umap_total

##########################################################

# Figure3A
fig3a_data <- subset(dat.num, type == "UD" & ratio != "1-1-3")
fig3a_data$T_cell_facet <- factor(fig3a_data$T_cell, labels = c("bold(CD4^{'+'})", "bold(CD8^{'+'})"))
fig3a_data$ratio_facet <- factor(fig3a_data$ratio, labels = c("bold('1:0')", "bold('1:4')"))
med_coord$T_cell_facet <- "bold(CD4^{'+'})"
med_coord$ratio_facet <- "bold('1:0')"

fig3a <- ggplot(fig3a_data, aes(UMAP_1, UMAP_2, col = pheno_cluster)) +
  geom_point(size = 0.35) +
  scale_color_discrete(drop = F) + # Keep same colour scheme was total cells graph
  facet_grid(T_cell_facet ~ ratio_facet, switch = "y", labeller = label_parsed) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = "none", # Remove legend
        strip.text = element_text(size = 12),
        strip.background = element_blank(),
        strip.placement = "outside") +
  scale_y_continuous(breaks = pretty_breaks(), limits = c(min(dat.num$UMAP_2), max(dat.num$UMAP_2))) + # To remove decimals
  scale_x_continuous(limits = c(min(dat.num$UMAP_1), max(dat.num$UMAP_1))) +
  geom_label_repel(data = subset(med_coord, T_cell_facet == "bold(CD4^{'+'})" & ratio_facet == "bold('1:0')"), aes(label = pheno_cluster), show.legend = F, size = 1.25, position = position_jitter()) # Add labels on map
fig3a

ggsave(file = "230825_Figure3A.svg", plot = fig3a, height = 9, width = 12, units = "cm")

##########################################################

# Figure 6A
fig6a_data <- subset(dat.num, type != "ART" & ratio == "1-4-0")
fig6a_data$T_cell_facet <- factor(fig6a_data$T_cell, labels = c("bold(CD4^{'+'})", "bold(CD8^{'+'})"))
fig6a_data$type_facet <- factor(fig6a_data$type, labels = c("bold('UD')", "bold('UNT')"))
med_coord$T_cell_facet <- "bold(CD4^{'+'})"
med_coord$type_facet <- "bold('UD')"

fig6a <- ggplot(fig6a_data, aes(UMAP_1, UMAP_2, col = pheno_cluster)) +
  geom_point(size = 0.35) +
  scale_color_discrete(drop = F) + # Keep same colour scheme was total cells graph
  facet_grid(T_cell_facet ~ type_facet, switch = "y", labeller = label_parsed) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = "none", # Remove legend
        strip.text = element_text(size = 12),
        strip.background = element_blank(),
        strip.placement = "outside") +
  scale_y_continuous(breaks = pretty_breaks(), limits = c(min(dat.num$UMAP_2), max(dat.num$UMAP_2))) + # To remove decimals
  scale_x_continuous(limits = c(min(dat.num$UMAP_1), max(dat.num$UMAP_1))) +
  geom_label_repel(data = subset(med_coord, T_cell_facet == "bold(CD4^{'+'})" & type_facet == "bold('UD')"), aes(label = pheno_cluster), show.legend = F, size = 1.25, position = position_jitter()) # Add labels on map
fig6a

ggsave(file = "230825_Figure6A.svg", plot = fig6a, width = 12, height = 9, units = "cm")

##########################################################

# Figure 7A
fig7a_data <- subset(dat.num, type == "ART" & ratio == "1-4-0")
fig7a_data$T_cell_facet <- factor(fig7a_data$T_cell, labels = c("bold(CD4^{'+'})", "bold(CD8^{'+'})"))
med_coord$T_cell_facet <- "bold(CD4^{'+'})"

fig7a <- ggplot(fig7a_data, aes(UMAP_1, UMAP_2, col = pheno_cluster)) +
  geom_point(size = 0.3) +
  scale_color_discrete(drop = F) + # Keep same colour scheme was total cells graph
  facet_grid(T_cell_facet ~ ., switch = "y", labeller = label_parsed) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = "none", # Remove legend
        strip.text = element_text(size = 12),
        strip.background = element_blank(),
        strip.placement = "outside") +
  scale_y_continuous(breaks = pretty_breaks(), limits = c(min(dat.num$UMAP_2), max(dat.num$UMAP_2))) + # To remove decimals
  scale_x_continuous(limits = c(min(dat.num$UMAP_1), max(dat.num$UMAP_1)))
  # geom_label_repel(data = subset(med_coord, T_cell_facet == "bold(CD4^{'+'})"), aes(label = pheno_cluster), show.legend = F, size = 1.75, position = position_jitter()) # Add labels on map
fig7a

ggsave(file = "230825_Figure7A.svg", plot = fig7a, height = 7, width = 6, units = "cm")

#################################################################################################

              ###############################################################
              # 8 "Manually" scale cytokine values at the single-cell level #
              ###############################################################

# Instead of using the MFI to create the UMAP, scale manually using the min and max FI of each cytokine
# Normalize between 0 and 1
# Formula for scaling is: zi = (xi – min(x)) / (max(x) – min(x))
# https://www.statology.org/normalize-data-between-0-and-1/

normalized_values <- dat.num

for (i in 9:13){
  max <- min <- NULL
  max <- max(dat.num[,i])
  min <- min(dat.num[,i])
  
  normalized_values[,i] <- sapply(dat.num[,i], function(x){
    (x-min)/(max-min)
  })
}

##########################################################

# Create custom palette for heatmap
custom_pal <- c("#000000", "#4C0C1A", "#F33248", "#F2A87E", "#E6D08A", "#F9D984", "#FCEEA8", "#FFFDD1", "#FFFDF0", "#FFFFFF")

# CD4
umap_cytokines_CD4_norm <- list()
for (i in 9:13){
  tmp <- colnames(normalized_values)[i]
  umap_cytokines_CD4_norm[[tmp]] <- ggplot(normalized_values[normalized_values$T_cell == "CD4",], aes_string("UMAP_1", "UMAP_2", col = tmp)) +
    geom_point(size = 0.2) +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5)) +
    labs(y = "", x = "", title = tmp, colour = "Value") +
    scale_color_gradientn(colours = rev(custom_pal), limits = c(min(normalized_values[,9:13]), max(normalized_values[,9:13])))
  # scale_color_gradientn(colours = rev(heat.colors(10)), limits = c(min(dat.num[,9:13]), max(dat.num[,9:13])))
}
fig_cytokines_CD4_norm <- do.call("ggarrange", c(umap_cytokines_CD4_norm, ncol = 5, nrow = 1, common.legend = T, legend = "right"))
annotate_figure(fig_cytokines_CD4_norm, bottom = "UMAP_1", left = "UMAP_2")
ggsave(file = "230825_FigureS3A_CD4.svg", plot = annotate_figure(fig_cytokines_CD4_norm, bottom = "UMAP_1", left = "UMAP_2"), width = 20, height = 5, units = "cm")

# CD8
umap_cytokines_CD8_norm <- list()
for (i in 9:13){
  tmp <- colnames(normalized_values)[i]
  umap_cytokines_CD8_norm[[tmp]] <- ggplot(normalized_values[normalized_values$T_cell == "CD8",], aes_string("UMAP_1", "UMAP_2", col = tmp)) +
    geom_point(size = 0.2) +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5)) +
    labs(y = "", x = "", title = tmp, colour = "Value") +
    scale_color_gradientn(colours = rev(custom_pal), limits = c(min(normalized_values[,9:13]), max(normalized_values[,9:13])))
  # scale_color_gradientn(colours = rev(heat.colors(10)), limits = c(min(dat.num[,9:13]), max(dat.num[,9:13])))
}
fig_cytokines_CD8_norm <- do.call("ggarrange", c(umap_cytokines_CD8_norm, ncol = 5, nrow = 1, common.legend = T, legend = "right"))
annotate_figure(fig_cytokines_CD8_norm, bottom = "UMAP_1", left = "UMAP_2")
ggsave(file = "230825_FigureS3A_CD8.svg", plot = annotate_figure(fig_cytokines_CD8_norm, bottom = "UMAP_1", left = "UMAP_2"), width = 20, height = 5, units = "cm")

##########################################################

# Organize data
sub2 <- normalized_values[,9:13] # Subset only phenotype channels
sub2$cluster <- normalized_values[,] # Add data
sub2 <- melt(sub2$cluster) # Reshape data w/ respect to variable cluster
# sub <- aggregate(sub$value, by = list(sub$cluster, sub$variable), mean)
sub2 <- aggregate(sub2$value, by = list(sub2$pheno_cluster, sub2$variable), mean) # Calculate the mean value for each cluster, for each variable.
d2 <- dcast(sub2, formula = Group.1 ~ Group.2) # Create df: rows = clusters; columns = variables
rownames(d2) <- d2[,1]
d2 <- d2[,-1]
d2 <- d2[,3:ncol(d2)]

ph_normalized <- as.ggplot(pheatmap(d2, border_col = 'white',
                                    cluster_col = F, cluster_row = F,
                                    color = rev(custom_pal), # Color palette can be changed
                                    angle_col = 90, # Turn the x-axis labels
                                    breaks = seq(0,1,0.1) # So that the legend stays from 0 to 1
))
ph_normalized

ggsave(file = "230825_Figure3B.svg", plot = ph_normalized, width = 6.30, height = 6.30)

##########################################################

dat.num$type.ratio <- paste0(dat.num$type, ".", dat.num$ratio)

# Save data as tables/CSV files
cluster_prop <- t(table(dat.num$pheno_cluster))
cluster_prop_T_cell <- t(table(dat.num$pheno_cluster, dat.num$T_cell))
cluster_prop_type_CD4 <- t(table(dat.num[dat.num$T_cell == "CD4",]$pheno_cluster, dat.num[dat.num$T_cell == "CD4",]$type))
cluster_prop_type_CD8 <- t(table(dat.num[dat.num$T_cell == "CD8",]$pheno_cluster, dat.num[dat.num$T_cell == "CD8",]$type))
cluster_prop_ratio_CD4 <- t(table(dat.num[dat.num$T_cell == "CD4",]$pheno_cluster, dat.num[dat.num$T_cell == "CD4",]$ratio))
cluster_prop_ratio_CD8 <- t(table(dat.num[dat.num$T_cell == "CD8",]$pheno_cluster, dat.num[dat.num$T_cell == "CD8",]$ratio))
cluster_prop_type.ratio_CD4 <- t(table(dat.num[dat.num$T_cell == "CD4",]$pheno_cluster, dat.num[dat.num$T_cell == "CD4",]$type.ratio))
cluster_prop_type.ratio_CD8 <- t(table(dat.num[dat.num$T_cell == "CD8",]$pheno_cluster, dat.num[dat.num$T_cell == "CD8",]$type.ratio))
cluster_donor_CD4 <- t(table(dat.num[dat.num$T_cell == "CD4",]$pheno_cluster, dat.num[dat.num$T_cell == "CD4",]$everything))
cluster_donor_CD8 <- t(table(dat.num[dat.num$T_cell == "CD8",]$pheno_cluster, dat.num[dat.num$T_cell == "CD8",]$everything))


# Save as CSV
CSV_filename <- paste0(PrimaryDirectory, "/", today, "_Tayma_ICS_CD4_CD8.csv")
write.table(cluster_prop, CSV_filename, sep = ",")
write.table(cluster_prop_T_cell, CSV_filename, sep = ",", append = T)
write.table(cluster_prop_type_CD4, CSV_filename, sep = ",", append = T)
write.table(cluster_prop_type_CD8, CSV_filename, sep = ",", append = T)
write.table(cluster_prop_ratio_CD4, CSV_filename, sep = ",", append = T)
write.table(cluster_prop_ratio_CD8, CSV_filename, sep = ",", append = T)
write.table(cluster_prop_type.ratio_CD4, CSV_filename, sep = ",", append = T)
write.table(cluster_prop_type.ratio_CD8, CSV_filename, sep = ",", append = T)
write.table(d2[,3:ncol(d2)], CSV_filename, sep = ",", append = T)
write.table(cluster_donor_CD4, CSV_filename, sep = ",", append = T)
write.table(cluster_donor_CD8, CSV_filename, sep = ",", append = T)

##########################################################

rm(i, name, order_params, tmp, count, colnames)
save.image(paste0(PrimaryDirectory, "/", today, "_unsupervised_ICS.RData")) # Save image

#################################################################################################
#################################################################################################


q(save="no")


#################################################################################################
#################################################################################################
#################################################################################################
