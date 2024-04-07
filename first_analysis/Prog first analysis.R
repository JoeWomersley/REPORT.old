setwd("/Users/marti/Desktop/Desktop/OneDrive/University of york third year/group project/core data analysis/R /omics data science /omic workshops ")
library(tidyverse)
prog <- read.csv("data/surfaceome_prog.csv")

# data distribution in expression levels per cells ------------------------
prog |> 
  pivot_longer(cols = -ensembl_gene_id,
               names_to = "cell",
               values_to = "expression") |> 
  ggplot(aes(x = expression)) +
  geom_histogram(fill = "blue", color = "black", binwidth = 0.5) +
  theme_minimal() +
  theme(panel.grid.major = element_blank(),  # Remove major grid lines
        panel.grid.minor = element_blank(),  # Remove minor grid lines
        axis.line.x = element_line(color = "black", size = 1),  # Adjust x-axis line
        axis.line.y = element_line(color = "black", size = 1),
        axis.text = element_text(size = 30), # Increase axis text size
        axis.title = element_text(size = 35),axis.line = element_line(size = 2), # Make axis lines thicker
        plot.title = element_text(size = 40, hjust = 0.5), plot.margin = margin(0.4, 0.4, 0.4, 0.4, "cm")) + scale_y_continuous(limits = c(0, 114300)) # Adjust y-axis line



# summarising for each prog the min and max expression level  --------

summary(prog)

summary(prog[1:20])
# ensembl_gene_id       Prog_001         Prog_002         Prog_003         Prog_004         Prog_006     
# Length:280         Min.   : 0.000   Min.   : 0.000   Min.   : 0.000   Min.   : 0.000   Min.   : 0.000  
# Class :character   1st Qu.: 0.000   1st Qu.: 0.000   1st Qu.: 0.000   1st Qu.: 0.000   1st Qu.: 0.000  
# Mode  :character   Median : 1.080   Median : 0.786   Median : 1.050   Median : 0.000   Median : 0.000  
# Mean   : 2.355   Mean   : 1.713   Mean   : 2.220   Mean   : 1.975   Mean   : 2.022  
# 3rd Qu.: 4.337   3rd Qu.: 1.666   3rd Qu.: 2.723   3rd Qu.: 2.049   3rd Qu.: 1.501  
# Max.   :10.360   Max.   :10.717   Max.   :11.340   Max.   :10.971   Max.   :13.609  
# Prog_007          Prog_008         Prog_009          Prog_010         Prog_011         Prog_012     
# Min.   : 0.0000   Min.   : 0.000   Min.   : 0.0000   Min.   : 0.000   Min.   : 0.000   Min.   : 0.000  
# 1st Qu.: 0.0000   1st Qu.: 0.000   1st Qu.: 0.0000   1st Qu.: 0.000   1st Qu.: 0.000   1st Qu.: 0.000  
# Median : 0.7008   Median : 0.000   Median : 0.6569   Median : 0.000   Median : 1.131   Median : 0.000  
# Mean   : 2.1447   Mean   : 2.009   Mean   : 1.9286   Mean   : 2.652   Mean   : 2.019   Mean   : 2.302  
# 3rd Qu.: 4.0009   3rd Qu.: 1.823   3rd Qu.: 1.5181   3rd Qu.: 6.087   3rd Qu.: 2.526   3rd Qu.: 3.588  
# Max.   :10.3219   Max.   :12.833   Max.   :11.8691   Max.   :11.938   Max.   :11.312   Max.   :10.666  
# Prog_013          Prog_014         Prog_015         Prog_016         Prog_017          Prog_018     
# Min.   : 0.0000   Min.   : 0.000   Min.   : 0.000   Min.   : 0.000   Min.   : 0.0000   Min.   : 0.000  
# 1st Qu.: 0.0000   1st Qu.: 0.000   1st Qu.: 0.000   1st Qu.: 0.000   1st Qu.: 0.0000   1st Qu.: 0.000  
# Median : 0.7286   Median : 0.000   Median : 0.000   Median : 0.000   Median : 0.6852   Median : 0.000  
# Mean   : 2.4302   Mean   : 2.032   Mean   : 2.595   Mean   : 2.260   Mean   : 2.6568   Mean   : 2.137  
# 3rd Qu.: 4.4002   3rd Qu.: 2.439   3rd Qu.: 6.317   3rd Qu.: 2.938   3rd Qu.: 6.0185   3rd Qu.: 2.616  
# Max.   :12.2556   Max.   :11.358   Max.   :11.891   Max.   :11.589   Max.   :10.5896   Max.   :12.538  
# Prog_019          Prog_020     
# Min.   : 0.0000   Min.   : 0.000  
# 1st Qu.: 0.0000   1st Qu.: 0.000  
# Median : 0.7909   Median : 0.000  
# Mean   : 2.6865   Mean   : 2.091  
# 3rd Qu.: 6.2670   3rd Qu.: 2.004  
# Max.   :11.9309   Max.   :11.129   

# summarising expression levels (looking at epxression levels of cells) -----------------------------------------------

# # below is a pointrange plot Pointrange puts a dot 
# at the mean and a line between a minimum and a maximum 
# such as +/- one s.d. Not unlike a boxplot, but when you 
# need the boxes too be very narrow!

prog_summary_samp <- prog |>
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

prog_summary_samp |> 
  ggplot(aes(x = cell, y = mean)) +
  geom_pointrange(aes(ymin = mean - sd, 
                      ymax = mean + sd ),
                  size = 0.1)

prog_summary_samp |> 
  ggplot(aes(x = cell, y = mean)) +
  geom_pointrange(aes(ymin = mean - sd, ymax = mean + sd),
                  size = 0.1,  # Adjust the size here for the points
                  color = "blue") +
  geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd),
                width = 0.1,  # Adjust the width of the lines here
                linewidth = 0.1,  # Adjust the size of the lines here
                color = "darkblue") +
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

prog_summary_samp |> 
  ggplot(aes(x = reorder(cell, mean), y = mean)) +
  geom_pointrange(aes(ymin = mean - sd, 
                      ymax = mean + sd ),
                  size = 0.1)
# so this is a plot with the cells reordered instead of being alphabetical they 
# are now in order of smallest to biggest mean.


prog_summary_samp |> 
  ggplot(aes(x = reorder(cell, mean), y = mean)) +
  geom_pointrange(aes(ymin = mean - sd, ymax = mean + sd),
                  size = 0.1,  # Adjust the size here for the points
                  color = "blue") +
  geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd),
                width = 0.1,  # Adjust the width of the lines here
                linewidth = 0.1,  # Adjust the size of the lines here
                color = "darkblue") +
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

prog_summary_gene <- prog |>
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

prog_summary_gene |> 
  ggplot(aes(x = reorder(ensembl_gene_id, mean), y = mean)) +
  geom_pointrange(aes(ymin = mean - sd, 
                      ymax = mean + sd ),
                  size = 0.1)

prog_summary_gene |> 
  ggplot(aes(x = reorder(ensembl_gene_id, mean), y = mean)) +
  geom_pointrange(aes(ymin = mean - sd, ymax = mean + sd, color = "blue"),
                  size = 0.1) +
  geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd, color = "darkblue"),
                width = 0.1) +
  theme_minimal() +
  theme(axis.line.x = element_line(color = "black", size = 1),  # Adjust x-axis line
        axis.line.y = element_line(color = "black", size = 1),
        axis.text = element_text(size = 30),  # Increase axis text size
        axis.title = element_text(size = 35),  # Increase axis title size
        axis.line = element_line(size = 2), 
        axis.text.x = element_blank()) +
  labs(x = "genes", y = "Log Mean Expression") +
  scale_color_identity()



write_csv(prog_summary_gene, "processed data/prog_summary_gene.csv")

# Note again that the variability between genes 
# (average expression between 0.02 and and 10.03) is far greater than 
# between cells (average expression from1.46 to 3.18) which is expected.


