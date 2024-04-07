setwd("/Users/marti/Desktop/Desktop/OneDrive/University of york 
      third year/group project/core data analysis/R /omics data science /omic workshops ")
library(tidyverse)
hspc <- read.csv("data/surfaceome_hspc.csv")


# data distribution in expression levels per cells ------------------------

hspc |> 
  pivot_longer(cols = -ensembl_gene_id,
               names_to = "cell",
               values_to = "expression") |> 
  ggplot(aes(x = expression)) +
  geom_histogram(fill = "red", color = "black", binwidth = 0.5) +
  theme_minimal() +
  theme(panel.grid = element_blank(),  # Remove grid lines
        axis.line.x = element_line(color = "black", size = 1),  # Adjust x-axis line
        axis.line.y = element_line(color = "black", size = 1),
        axis.text = element_text(size = 30),  # Increase axis text size
        axis.title = element_text(size = 35),  # Increase axis title size
        axis.line = element_line(size = 2),  # Make axis lines thicker
        plot.title = element_text(size = 40, hjust = 0.5),  # Adjust plot title size and alignment
        plot.margin = margin(0.4, 0.4, 0.4, 0.4, "cm"))  # Adjust plot margins to move it downwards slightly

# summarising for each HSPC the min and max expression level  --------

summary(hspc)
summary(hspc[1:20])
# ensembl_gene_id       HSPC_001         HSPC_002         HSPC_003          HSPC_004         HSPC_006     
# Length:280         Min.   : 0.000   Min.   : 0.000   Min.   : 0.0000   Min.   : 0.000   Min.   : 0.000  
# Class :character   1st Qu.: 0.000   1st Qu.: 0.000   1st Qu.: 0.0000   1st Qu.: 0.000   1st Qu.: 0.000  
# Mode  :character   Median : 0.000   Median : 0.000   Median : 0.9929   Median : 0.000   Median : 1.276  
# Mean   : 2.143   Mean   : 1.673   Mean   : 2.5964   Mean   : 1.851   Mean   : 2.338  
# 3rd Qu.: 2.120   3rd Qu.: 2.239   3rd Qu.: 6.1559   3rd Qu.: 2.466   3rd Qu.: 3.536  
# Max.   :12.567   Max.   :11.976   Max.   :11.1138   Max.   :11.133   Max.   :10.014  
# HSPC_008         HSPC_009        HSPC_011         HSPC_012         HSPC_014         HSPC_015     
# Min.   : 0.000   Min.   :0.000   Min.   : 0.000   Min.   : 0.000   Min.   : 0.000   Min.   : 0.000  
# 1st Qu.: 0.000   1st Qu.:0.000   1st Qu.: 0.000   1st Qu.: 0.000   1st Qu.: 0.000   1st Qu.: 0.000  
# Median : 0.000   Median :0.000   Median : 0.000   Median : 1.750   Median : 0.000   Median : 0.000  
# Mean   : 2.375   Mean   :2.220   Mean   : 2.285   Mean   : 2.431   Mean   : 2.295   Mean   : 2.515  
# 3rd Qu.: 3.851   3rd Qu.:3.594   3rd Qu.: 3.193   3rd Qu.: 3.741   3rd Qu.: 3.150   3rd Qu.: 3.789  
# Max.   :11.574   Max.   :9.997   Max.   :11.260   Max.   :10.905   Max.   :11.051   Max.   :10.751  
# HSPC_016          HSPC_017         HSPC_018         HSPC_020         HSPC_021         HSPC_022     
# Min.   : 0.0000   Min.   : 0.000   Min.   : 0.000   Min.   : 0.000   Min.   : 0.000   Min.   : 0.000  
# 1st Qu.: 0.0000   1st Qu.: 0.000   1st Qu.: 0.000   1st Qu.: 0.000   1st Qu.: 0.000   1st Qu.: 0.000  
# Median : 0.9488   Median : 0.000   Median : 1.248   Median : 0.000   Median : 0.000   Median : 0.000  
# Mean   : 2.6115   Mean   : 2.146   Mean   : 2.710   Mean   : 2.509   Mean   : 2.170   Mean   : 2.287  
# 3rd Qu.: 5.9412   3rd Qu.: 2.357   3rd Qu.: 6.006   3rd Qu.: 4.470   3rd Qu.: 2.996   3rd Qu.: 3.351  
# Max.   :11.3082   Max.   :12.058   Max.   :11.894   Max.   :11.281   Max.   :10.709   Max.   :11.814  
# HSPC_023         HSPC_024     
# Min.   : 0.000   Min.   : 0.000  
# 1st Qu.: 0.000   1st Qu.: 0.000  
# Median : 0.000   Median : 0.000  
# Mean   : 2.314   Mean   : 2.195  
# 3rd Qu.: 2.749   3rd Qu.: 2.944  
# Max.   :12.113   Max.   :11.279  

# summarising expression levels (looking at epxression levels of cells) -----------------------------------------------

# # below is a pointrange plot Pointrange puts a dot 
# at the mean and a line between a minimum and a maximum 
# such as +/- one s.d. Not unlike a boxplot, but when you 
# need the boxes too be very narrow!

hspc_summary_samp <- hspc |>
  pivot_longer(cols = -ensembl_gene_id,
               names_to = "cell",
               values_to = "expr") |>
  group_by(cell) |>
  summarise(min = min(expr),
            lowerq = quantile(expr, 0.25),
            mean = mean(expr),
            median = median(expr),
            sd = sd(expr),
            upperq = quantile(expr, 0.75),
            max = max(expr),
            n_zero = sum(expr == 0))

hspc_summary_samp |> 
  ggplot(aes(x = cell, y = mean)) +
  geom_pointrange(aes(ymin = mean - sd, 
                      ymax = mean + sd ),
                  size = 0.1)

hspc_summary_samp |> 
  ggplot(aes(x = cell, y = mean)) +
  geom_pointrange(aes(ymin = mean - sd, ymax = mean + sd),
                  size = 0.1,  # Adjust the size here for the points
                  color = "red") +
  geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd),
                width = 0.1,  # Adjust the width of the lines here
                linewidth = 0.1,  # Adjust the size of the lines here
                color = "darkred") +
  theme_minimal() +
  theme( axis.line.x = element_line(color = "black", size = 1),  # Adjust x-axis line
         axis.line.y = element_line(color = "black", size = 1),
         axis.text = element_text(size = 30),  # Increase axis text size
         axis.title = element_text(size = 35),  # Increase axis title size
         axis.line = element_line(size = 2), axis.text.x = element_blank()) # Remove x-axis text


# You will need to use the Zoom button to pop the plot window out so you can
# make it as wide as possible
# 
# The things to notice are:
#   
# the average expression in cells is similar for all cells. 
# This is good to know - if some cells had much lower expression
# perhaps there is something wrong with them, or their sequencing, 
# and they should be excluded.
# # the distributions are roughly similar in width too
# # The default order of cell is alphabetical. It can be easier to 
# see these (non-) effects if we order the lines by the size of the mean.

hspc_summary_samp |> 
  ggplot(aes(x = reorder(cell, mean), y = mean)) +
  geom_pointrange(aes(ymin = mean - sd, 
                      ymax = mean + sd ),
                  size = 0.1)
# so this is a plot with the cells reordered instead of being alphabetical they 
# are now in order of smallest to biggest mean.

hspc_summary_samp |> 
  ggplot(aes(x = reorder(cell, mean), y = mean)) +
  geom_pointrange(aes(ymin = mean - sd, ymax = mean + sd),
                  size = 0.1,  # Adjust the size here for the points
                  color = "red") +
  geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd),
                width = 0.1,  # Adjust the width of the lines here
                linewidth = 0.1,  # Adjust the size of the lines here
                color = "darkred") +
  theme_minimal() +
  theme( axis.line.x = element_line(color = "black", size = 1),  # Adjust x-axis line
         axis.line.y = element_line(color = "black", size = 1),
         axis.text = element_text(size = 30),  # Increase axis text size
         axis.title = element_text(size = 35),  # Increase axis title size
         axis.line = element_line(size = 2), axis.text.x = element_blank()) +
  labs(x = "Cells reordered", y = "Mean Expression") 


# summarising expression levels (looking at each gene) --------------------------

# There are fewer genes in this dataset, but still more than you can understand
# without the overview provided by a plot. We will again pivot the data to 
# tidy and then summarise the expression for each gene.
# 
# 🎬 Summarise the expression for each genes:

hspc_summary_gene <- hspc |>
  pivot_longer(cols = -ensembl_gene_id,
               names_to = "cell",
               values_to = "expr") |>
  group_by(ensembl_gene_id) |>
  summarise(min = min(expr),
            lowerq = quantile(expr, 0.25),
            sd = sd(expr),
            mean = mean(expr),
            median = median(expr),
            upperq = quantile(expr, 0.75),
            max = max(expr),
            total = sum(expr),
            n_zero = sum(expr == 0))

# View the hspc_summary_gene dataframe. Remember these are normalised and 
# logged (base 2) so we should not see very large values.
# 
# Notice that we have:
#   
#   no genes with 0 in every cell
# very few genes (9) with no zeros at all
# quite a few genes with zero in many cells but this matters less than zeros 
# in the frog samples because we had just 6 samples and we have 701 cells.
# As we have a lot of genes, it is again helpful to plot the mean expression 
# with pointrange to get an overview. We do not need to log the values but 
# ordering the genes will help.
# 
# 🎬 Plot the logged mean counts for each gene in order of size using 
# geom_pointrange():

hspc_summary_gene |> 
  ggplot(aes(x = reorder(ensembl_gene_id, mean), y = mean)) +
  geom_pointrange(aes(ymin = mean - sd, 
                      ymax = mean + sd ),
                  size = 0.1)

hspc_summary_gene |> 
  ggplot(aes(x = reorder(ensembl_gene_id, mean), y = mean)) +
  geom_pointrange(aes(ymin = mean - sd, ymax = mean + sd, color = "red"),
                  size = 0.1) +
  geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd, color = "darkred"),
                width = 0.1) +
  theme_minimal() +
  theme(axis.line.x = element_line(color = "black", size = 1),  # Adjust x-axis line
        axis.line.y = element_line(color = "black", size = 1),
        axis.text = element_text(size = 30),  # Increase axis text size
        axis.title = element_text(size = 35),  # Increase axis title size
        axis.line = element_line(size = 2), 
        axis.text.x = element_blank()) +
  labs(x = "genes", y = "Mean Expression") +
  scale_color_identity()



# Note again that the variability between genes 
# (average expression between 0.02 and and 10.03) is far greater than 
# between cells (average expression from1.46 to 3.18) which is expected.

write_csv(hspc_summary_gene, "processed data/hspc_summary_gene.csv")


install.packages(BiocManager)
BiocManager::install("DESeq2")
BiocManager::install("scran")




