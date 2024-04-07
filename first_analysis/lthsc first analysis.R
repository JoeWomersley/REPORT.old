setwd("/Users/marti/Desktop/Desktop/OneDrive/University of york third year/group project/core data analysis/R /omics data science /omic workshops ")
library(tidyverse)
lthsc <- read.csv("data/surfaceome_lthsc.csv")

# data distribution in expression levels per cells ------------------------

lthsc |> 
  pivot_longer(cols = -ensembl_gene_id,
               names_to = "cell",
               values_to = "expression") |> 
  ggplot(aes(x = expression)) +
  geom_histogram(fill = "green", color = "black", binwidth = 0.5) +
  theme_minimal() +
  theme(panel.grid = element_blank(),  # Remove grid lines
        axis.line.x = element_line(color = "black", size = 1),  # Adjust x-axis line
        axis.line.y = element_line(color = "black", size = 1),
        axis.text = element_text(size = 30),  # Increase axis text size
        axis.title = element_text(size = 35),  # Increase axis title size
        axis.line = element_line(size = 2),  # Make axis lines thicker
        plot.title = element_text(size = 40, hjust = 0.5),  # Adjust plot title size and alignment
        plot.margin = margin(0.4, 0.4, 0.4, 0.4, "cm"))

# summarising for each ltht the min and max expression level  --------

summary(lthsc)

summary(lthsc[1:20])
# ensembl_gene_id      LT.HSC_001       LT.HSC_002       LT.HSC_003       LT.HSC_004       LT.HSC_006    
# Length:280         Min.   : 0.000   Min.   : 0.000   Min.   : 0.000   Min.   : 0.000   Min.   : 0.000  
# Class :character   1st Qu.: 0.000   1st Qu.: 0.000   1st Qu.: 0.000   1st Qu.: 0.000   1st Qu.: 0.000  
# Mode  :character   Median : 0.000   Median : 1.404   Median : 0.000   Median : 0.000   Median : 0.000  
# Mean   : 2.273   Mean   : 2.557   Mean   : 2.242   Mean   : 2.288   Mean   : 1.960  
# 3rd Qu.: 3.474   3rd Qu.: 5.369   3rd Qu.: 3.412   3rd Qu.: 3.213   3rd Qu.: 2.496  
# Max.   :11.367   Max.   :10.274   Max.   :11.433   Max.   :10.322   Max.   :11.321  
# LT.HSC_009       LT.HSC_010       LT.HSC_012      LT.HSC_013       LT.HSC_014       LT.HSC_015    
# Min.   : 0.000   Min.   : 0.000   Min.   :0.000   Min.   : 0.000   Min.   : 0.000   Min.   : 0.000  
# 1st Qu.: 0.000   1st Qu.: 0.000   1st Qu.:0.000   1st Qu.: 0.000   1st Qu.: 0.000   1st Qu.: 0.000  
# Median : 0.000   Median : 0.000   Median :1.528   Median : 1.871   Median : 0.000   Median : 0.000  
# Mean   : 2.093   Mean   : 2.331   Mean   :2.172   Mean   : 2.205   Mean   : 2.448   Mean   : 2.162  
# 3rd Qu.: 2.087   3rd Qu.: 3.116   3rd Qu.:2.824   3rd Qu.: 3.166   3rd Qu.: 3.312   3rd Qu.: 3.053  
# Max.   :10.937   Max.   :12.019   Max.   :9.758   Max.   :10.244   Max.   :11.792   Max.   :10.682  
# LT.HSC_016       LT.HSC_017       LT.HSC_018       LT.HSC_020       LT.HSC_021       LT.HSC_022    
# Min.   : 0.000   Min.   : 0.000   Min.   : 0.000   Min.   : 0.000   Min.   : 0.000   Min.   : 0.000  
# 1st Qu.: 0.000   1st Qu.: 0.000   1st Qu.: 0.000   1st Qu.: 0.000   1st Qu.: 0.000   1st Qu.: 0.000  
# Median : 1.657   Median : 0.000   Median : 0.000   Median : 1.197   Median : 1.001   Median : 0.000  
# Mean   : 2.378   Mean   : 2.110   Mean   : 2.152   Mean   : 2.362   Mean   : 2.356   Mean   : 1.983  
# 3rd Qu.: 3.265   3rd Qu.: 2.333   3rd Qu.: 2.207   3rd Qu.: 3.180   3rd Qu.: 3.419   3rd Qu.: 3.121  
# Max.   :11.149   Max.   :11.666   Max.   :11.290   Max.   :10.582   Max.   :11.116   Max.   :11.668  
# LT.HSC_023       LT.HSC_024    
# Min.   : 0.000   Min.   : 0.000  
# 1st Qu.: 0.000   1st Qu.: 0.000  
# Median : 1.625   Median : 0.000  
# Mean   : 2.219   Mean   : 2.009  
# 3rd Qu.: 2.859   3rd Qu.: 2.268  
# Max.   :11.804   Max.   :12.929  

# summarising expression levels (looking at epxression levels of cells) -----------------------------------------------

# # below is a pointrange plot Pointrange puts a dot 
# at the mean and a line between a minimum and a maximum 
# such as +/- one s.d. Not unlike a boxplot, but when you 
# need the boxes too be very narrow!

lthsc_summary_samp <- lthsc |>
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


lthsc_summary_samp |> 
  ggplot(aes(x = cell, y = mean)) +
  geom_pointrange(aes(ymin = mean - sd, ymax = mean + sd),
                  size = 0.5,  # Adjust the size here for the points
                  color = "darkgreen") +
  geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd),
                width = 0.1,  # Adjust the width of the lines here
                linewidth = 1.5,  # Adjust the size of the lines here
                color = "darkgreen") +
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

lthsc_summary_samp |> 
  ggplot(aes(x = reorder(cell, mean), y = mean)) +
  geom_pointrange(aes(ymin = mean - sd, 
                      ymax = mean + sd ),
                  size = 0.1)
# so this is a plot with the cells reordered instead of being alphabetical they 
# are now in order of smallest to biggest mean.

lthsc_summary_samp |> 
  ggplot(aes(x = reorder(cell, mean), y = mean)) +
  geom_pointrange(aes(ymin = mean - sd, ymax = mean + sd),
                  size = 0.5,  # Adjust the size here for the points
                  color = "green") +
  geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd),
                width = 0.1,  # Adjust the width of the lines here
                linewidth = 0.1,  # Adjust the size of the lines here
                color = "darkgreen") +
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

lthsc_summary_gene <- lthsc |>
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

lthsc_summary_gene |> 
  ggplot(aes(x = reorder(ensembl_gene_id, mean), y = mean)) +
  geom_pointrange(aes(ymin = mean - sd, 
                      ymax = mean + sd ),
                  size = 0.1)

lthsc_summary_gene |> 
  ggplot(aes(x = reorder(ensembl_gene_id, mean), y = mean)) +
  geom_pointrange(aes(ymin = mean - sd, ymax = mean + sd, color = "green"),
                  size = 0.1) +
  geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd, color = "darkgreen"),
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






