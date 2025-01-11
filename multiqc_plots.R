library("ggplot2")
library(stringr)
library(forcats)
library(gridExtra)
library(dplyr)

#adjust to preferred project directory
project_path <- "C:/Users/yanni/OneDrive - Universitaet Bern/Master"

#set working directory to project directory
setwd(project_path)

#read in summary statistics data file for featureCounts
featCounts <- read.csv(file = paste0(project_path,"/featureCounts_assignment_plot.csv"), header = T)

#extract sample numbers and add SRR as a prefix to sample IDs
featCounts$Category <- paste0("SRR",str_extract_all(featCounts$Category, "\\d+"))

#create new dataframe to make it easier to plot, each sample will have 5 entries --> 5*16=80
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
  geom_bar(position = "fill", stat = "identity") + #add barplot
  scale_fill_viridis_d(option = "plasma", direction = -1) + #set color palette
  theme_bw() + #adjust theme
  labs(x = "Percentages") + #rename x-axis
  theme(axis.title.x = element_text(size = 16), #adjust size of axes-titles and ticks, remove redundant y-axis-label and legend-title
        axis.title.y = element_blank(),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        legend.text = element_text(size = 13),
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
df$percentages <- df$values / df$total

#calculate mean percentage of aligned reads
mean_assigned <- mean(df$percentages[df$type == "Assigned"]) * 100

#calculate standard error of the mean percentage of aligned reads
std_error <- sd(df$percentages[df$type == "Assigned"]) / sqrt(length(df$percentages[df$type == "Assigned"])) * 100

#calculate mean percentage of ambiguous reads
mean_ambiguity <- mean(df$percentages[df$type == "Unassigned..Ambiguity"]) * 100

#calculate standard error of the mean percentage of ambiguous reads
std_error_ambiguity <- sd(df$percentages[df$type == "Unassigned..Ambiguity"]) / sqrt(length(df$percentages[df$type == "Unassigned..Ambiguity"])) * 100

#-----------------------------------------------------------------------------------------------------------------------------------

#read in summary statistics of mapping step
mapping <- read.csv(file = paste0(project_path,"/bowtie2_pe_plot.csv"), header = T)

#remove redundant last row
mapping <- mapping[-17,]

#create empty data frame to fit data, each sample will have 6 different mapping types --> 6*16=96
mapping_df <- data.frame(Category = rep(NA, 96))

#extract sample id
mapping_df$Category <- rep(paste0("SRR",78219,str_extract_all(mapping$Category, "\\d+")), each = 6)

#extract the mapping types
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
  geom_bar(position = "fill", stat = "identity") + #create barplot
  scale_fill_viridis_d(option = "plasma", direction = -1) + #set color palette
  theme_bw() + #adjust theme
  labs(x = "Percentages") + #adjust x-axis-label
  theme(axis.title.x = element_text(size = 16), #adjust size off axis titles and ticks, remove redundant legend title and y-axis-title
        axis.title.y = element_blank(),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        legend.text = element_text(size = 13),
        legend.title = element_blank()) +
  guides(fill = guide_legend(reverse = TRUE)) #adjust legend to show filled squares

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

#number of samples: 4*normal because two time points (before cleaning and after cleaning) and two mates.
number_of_samples <- 16*4

#create empty data frame to fit reformatted data
fastqc_df <- data.frame(sample_id = rep(NA,length(fastqc$Position..bp.)*number_of_samples)) 

#extract sample_id from original data frame, remove information to be stored in other columns
fastqc_df$sample_id <- rep(sub("_.*", "", colnames(fastqc)[-1]), each = length(fastqc$Position..bp.))

#extract position from original data frame
fastqc_df$position <- rep(fastqc$Position..bp.,number_of_samples)

#extract before cleaning, after cleaning information from the sample names from original dataframe
fastqc_df$time_point <- rep(sapply(colnames(fastqc)[-1], function(x){if(grepl(x, pattern = "_cleaned")){x <- "after"}else{x <- "before"}}), each = length(fastqc$Position..bp.))
#convert time information to factor
fastqc_df$time_point <- as.factor(fastqc_df$time_point)

#extract mate information from original data frame
fastqc_df$mate <- rep(sapply(colnames(fastqc)[-1], function(x){if(grepl(x, pattern = "_1")){x <- 1}else{x <- 2}}), each = length(fastqc$Position..bp.))
#convert mate column to factor
fastqc_df$mate <- as.factor(fastqc_df$mate)

#extract all numerical values from the fastqc data frame
numerical_values <- apply(fastqc[,2:(number_of_samples+1)], 1, function(column) {
  as.numeric(column[!is.na(column)])  # Convert to numeric and remove NA values
})

#create empty vector to be filled in the for loop
values <- numeric()

#extract column from the numerical values vector and add it to values
for (i in 1:number_of_samples) {
  values <- c(values, numerical_values[i,])
}

#add phred_score values to dataframe
fastqc_df$phred_score <- values

#order time_points, so that before appears to the left of after in the plot
fastqc_df <- fastqc_df %>%
  mutate(time_point = factor(time_point, levels = c("before", "after")))

#create ggplot
fastqc_plot <- ggplot(fastqc_df, mapping = aes(x = position, y = phred_score))+
  geom_point(aes(color = mate))+  #add data points
  facet_wrap(~ time_point)+   #separate by time, creates two plots one before- and one after cleaning of the reads
  scale_color_manual(values = c("#0072B2", "#D95F02"))+ #define nice colors for the data
  geom_smooth(se = F, aes(color = mate))+ #add trendlines of the data
  geom_abline(slope = 0, intercept = 28, linetype = 'dashed')+ #add quality thresholds from multiqc as black dashed lines
  geom_abline(slope = 0, intercept = 20, linetype = 'dashed')+
  theme_bw()+ #adjust theme of the plot
  theme(strip.text = element_text(size = 16),  #increase text size of text appearing on the plot
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 14))+
  scale_y_continuous(limits = c(0,45), breaks = seq(0,45, by = 10))+ #adjust axes and axes-ticks
  scale_x_continuous(breaks = seq(0,76, by = 10))+
  ylab('Phred Score')+ #rename axes
  xlab('Position (bp)')+
  labs(color = 'Mate') #rename legend title

#arrange all three created plots into one, fastqc plot should fill out first row
grid.arrange(fastqc_plot, mapping_plot, featCounts_plot, nrow = 2, 
             layout_matrix = rbind(c(1, 1), c(2, 3)))
