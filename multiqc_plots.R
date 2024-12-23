library("ggplot2")
library(stringr)
library(forcats)

#adjust to preferred project directory
project_path <- "C:/Users/yanni/OneDrive - Universitaet Bern/Master"

#set working directory to project directory
setwd(project_path)

#read in summary statistics data file for featureCounts
featCounts <- read.csv(file = paste0(project_path,"/featureCounts_assignment_plot.csv"), header = T)

#extract sample numbers and add SRR as a prefix to sample IDs
featCounts$Category <- paste0("SRR",str_extract_all(featCounts$Category, "\\d+"))

#create new dataframe to make it easier to plot
df <- data.frame(Category = rep(NA, 80))

#add sample ids to df each repeated 5 times
df$Category <- rep(featCounts$Category, each = 5)

#add types of assignment repeated 16 times 
df$type <- rep(colnames(featCounts[-1]), by = 16)

#extract all numerical values from the featCounts data frame
numerical_values <- apply(featCounts[,2:6], 1, function(row) {
  as.numeric(row[!is.na(row)])  # Convert to numeric and remove NA values
})

#create empty vector to be filled in the for loop
values <- c()

#extract each row from featCounts in the numerical values vector and add it to values
for (i in 1:16) {
  values <- c(values, numerical_values[,i])
}

#add values to the dataframe
df$values <- values

#create stacked percentage barplot for featCounts
featCounts_plot <- ggplot(df, mapping = aes(fill = fct_reorder(type, values), x = values, y = Category)) +
  geom_bar(position = "fill", stat = "identity") +
  scale_fill_viridis_d(option = "plasma", direction = -1) +
  theme_bw() +
  labs(x = "Percentages") +
  theme(axis.title.x = element_text(size = 16),
        axis.title.y = element_blank(),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        legend.text = element_text(size = 14),
        legend.title = element_blank()) +
  guides(fill = guide_legend(reverse = TRUE))

#initialize total vector to store total number of reads
total <- c()

for(i in 1:16){
  #calculate total number of reads for each sample and add to vector
  total <- c(total, sum(featCounts[i,-1]))
} 

#add total number of reads to df
df$total <- rep(total, each = 5)

#calculate percentages
df$percentages = df$values / df$total

#calculate mean percentage of aligned reads
mean_assigned = mean(df$percentages[df$type == "Assigned"]) * 100

#calculate standard error of the mean percentage of aligned reads
std_error = sd(df$percentages[df$type == "Assigned"]) / sqrt(length(df$percentages[df$type == "Assigned"])) * 100

#-----------------------------------------------------------------------------------------------------------------------------------

#read in summary statistics of mapping step
mapping <- read.csv(file = paste0(project_path,"/bowtie2_pe_plot.csv"), header = T)

#remove redundant last row
mapping <- mapping[-17,]

mapping_df <- data.frame(Category = rep(NA, 96))

mapping_df$Category <- rep(paste0("SRR",78219,str_extract_all(mapping$Category, "\\d+")), each = 6)

mapping_df$type <- rep(colnames(mapping[-1]), by = 16)

#extract all numerical values from the mapping data frame
numerical_values <- apply(mapping[,2:7], 1, function(row) {
  as.numeric(row[!is.na(row)])  # Convert to numeric and remove NA values
})

#create empty vector to be filled in the for loop
values <- c()

#extract each row from featCounts in the numerical values vector and add it to values
for (i in 1:16) {
  values <- c(values, numerical_values[,i])
}

mapping_df$values <- values

#create stacked percentage barplot for mapping
mapping_plot <- ggplot(mapping_df, mapping = aes(fill = fct_reorder(type, values), x = values, y = Category)) +
  geom_bar(position = "fill", stat = "identity") +
  scale_fill_viridis_d(option = "plasma", direction = -1) +
  theme_bw() +
  labs(x = "Percentages") +
  theme(axis.title.x = element_text(size = 16),
        axis.title.y = element_blank(),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        legend.text = element_text(size = 14),
        legend.title = element_blank()) +
  guides(fill = guide_legend(reverse = TRUE))

#add featureCount and mapping plot together
aplot::plot_list(mapping_plot, featCounts_plot, ncol = 1, tag_levels = "A")

#initialize total vector to store total number of reads
total <- c()

for(i in 1:16){
  #calculate total number of reads for each sample and add to vector
  total <- c(total, sum(mapping[i,-1]))
} 

#add total number of reads to mapping
mapping_df$total <- rep(total, each = 6)

#calculate percentages
mapping_df$percentages <- mapping_df$values / mapping_df$total

overall_alignment <- rep(1, 16) - mapping_df$percentages[mapping_df$type == "PE.neither.mate.aligned"]

#calculate mean percentage of overall aligned reads
mean_overall_assigned <- mean(overall_alignment) * 100

#calculate standard error of the mean percentage of overall aligned reads
std_error <- sd(overall_alignment) / sqrt(length(overall_alignment)) * 100

#calculate mean percentage of concordantly aligned reads
mean_assigned_concordantly <- mean(mapping_df$percentages[which(mapping_df$type == "PE.mapped.uniquely")]) * 100

#calculate standard error of the mean percentage of overall aligned reads
std_error_concordantly <- sd(mapping_df$percentages[which(mapping_df$type == "PE.mapped.uniquely")])/sqrt(length(mapping_df$percentages[which(mapping_df$type == "PE.mapped.uniquely")])) * 100

#-----------------------------------------------------------------------------------------------------------------------------------

#read in summary statistics of fastqc
fastqc <- read.csv(file = paste0(project_path,"/fastqc_per_base_sequence_quality_plot.csv"), header = T)

#adjust plot visualization to show to plots per row, make plots bigger with cex
par(mfrow = c(1,2), cex = 1.3)

#initialize x-axis
x = fastqc$Position..bp.
#initialize plot by plotting first sample
plot(x, y = fastqc$SRR7821918_1, type = 'l', col = "firebrick", xlim = c(1, max(fastqc$Position..bp.)), ylim = c(0, 50), xlab = "Position (bp)", ylab = "Phred Score", main = "Before cleaning")

#iterate column wise through fastqc dataframe, save column name, if its uncleaned mate 1 make a red line on the plot, if uncleaned mate 2 make blue line
#if not corresponding to uncleaned samples, do nothing
for(i in 3:length(colnames(fastqc))){
  col_name = colnames(fastqc)[i]
  if (grepl("_1$", col_name)) {
    lines(x, y = fastqc[,i], col = "firebrick")
  } else if (grepl("_2$", col_name))  {
    lines(x, y = fastqc[,i], col = "dodgerblue")
  } else{}
}
#add black lines corresponding to quality thresholds
abline(h = c(28, 20), col = "black", lty = 2)
#add figure legend
legend("bottomright", legend = c("Mate 1", "Mate 2"), col = c("firebrick", "dodgerblue"), lty = c(1,1), cex = 0.9)

#initialize second plot for cleaned reads
plot(x, y = fastqc$SRR7821918_1_cleaned, type = 'l', col = "firebrick", xlim = c(1, max(fastqc$Position..bp.)), ylim = c(0, 50), xlab = "Position (bp)", ylab = "Phred Score", main = "After cleaning")

#iterate column wise through fastqc dataframe, save column name, if its cleaned mate 1 make a red line on the plot, if cleaned mate 2 make blue line
#if not corresponding to cleaned samples, do nothing
for(i in 3:length(colnames(fastqc))){
  col_name = colnames(fastqc)[i]
if (grepl("_1_cleaned$", col_name)) {
  lines(x, y = fastqc[,i], col = "firebrick")
} else if(grepl("_2_cleaned$", col_name)) {
  lines(x, y = fastqc[,i], col = "dodgerblue")
} else{}
}

#add black lines corresponding to quality thresholds
abline(h = c(28, 20), col = "black", lty = 2)
#add figure legend
legend("bottomright", legend = c("Mate 1", "Mate 2"), col = c("firebrick", "dodgerblue"), lty = c(1,1), cex = 0.9)
