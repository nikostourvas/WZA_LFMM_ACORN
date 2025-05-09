library(ggplot2)
library(data.table)
library(reshape2)
library(qvalue)
library(tidyr)

FDR_LEVEL=snakemake@params[[2]]

envfactors = readLines(snakemake@input[[1]])
envfactors = factor(envfactors, levels=envfactors)
print(envfactors)

WZA_res = list()
for (i in envfactors){
  WZA_res[[i]] = fread(paste0(snakemake@params[[1]], i, "_WZA_output.csv"))
}

# Calculate the genomic inflation factor (GIF) based on all the Z-scores.
gif = list()
for (i in envfactors){
  gif[[i]] = median((WZA_res[[i]]$Z)^2)/(qchisq(0.5, df = 1, lower.tail = FALSE)) 
}
# Save the GIFs to a file
gif_table = reshape2::melt(gif)
write.table(gif_table, file=paste0(snakemake@params[[1]], "GIF.csv"), sep=",", row.names=FALSE, quote=FALSE)

# Calculate the gif corrected p-values and q-values
print("Calculating p-values")
# for some reason, while testing with a toy dataset, WZA produced few NAs
# Until we investigate further, I replace those NAs with 1
WZA_res = lapply(WZA_res, function(x) {x$Z_pVal[is.na(x$Z_pVal)] = 1; return(x)})
for (i in envfactors){
  WZA_res[[i]]$Z_pVal_gif_adj = pchisq(WZA_res[[i]]$Z^2/gif[[i]], df = 1, lower.tail = FALSE)
}

print("Calculating q-values")
# FDR correction
WZA_res = lapply(WZA_res, function(x) {x$qvalue = qvalue(x$Z_pVal, fdr.level=FDR_LEVEL)$qvalues; return(x)})
WZA_res = lapply(WZA_res, function(x) {x$qvalue_gif_adj = qvalue(x$Z_pVal_gif_adj, fdr.level=FDR_LEVEL)$qvalues; return(x)})

# manipulate the index names to match the output from BayPass
WZA_res = lapply(WZA_res, function(x) {x$index = gsub("H2.3-Sc", "Qrob_H2.3.Sc", x$index); return(x)})
WZA_res = lapply(WZA_res, function(x) {x$index = gsub("Chr0", "", x$index); return(x)})
WZA_res = lapply(WZA_res, function(x) {x$index = gsub("Chr", "", x$index); return(x)})

# export each table to a file
for (i in envfactors){
  write.table(WZA_res[[i]], file=paste0(snakemake@params[[1]], i, "_WZA_output_fdr.csv"), 
              sep=",", row.names=FALSE, quote=FALSE)
}

############################
# Find significant windows
############################

WZA_res = lapply(WZA_res, function(x) x[ ,c(1,6,7,8,9,10)]) # keep only the columns we need for plotting
WZA_res = reshape2::melt(WZA_res, id.vars = c("index", "POS"))
# Make the dataset wide. Spread the "variable" column into multiple columns
WZA_res = spread(WZA_res, variable, value)
print(head(WZA_res))

print("string manipulation")
# remove from rows of the column WZA_res$index the suffix "_window_" and anything after that
WZA_res$CHR = as.factor(gsub("_window_.*", "", WZA_res$index))

colnames(WZA_res) = c("index", "POS", "envfactor", 
                      "pvalue", "gif_cor_pvalue", 
                      "qvalue", "qvalue_gif_adj", "CHR")
#write.table(WZA_res, file=paste0(snakemake@params[[1]], "WZA_test.csv"), sep=",", row.names=FALSE, quote=FALSE)

WZA_res$CHR = factor(WZA_res$CHR, levels=unique(WZA_res$CHR))

print("Finding significant windows")
WZA_q0.1 = WZA_res[WZA_res$qvalue_gif_adj < 0.1,]
fwrite(WZA_q0.1, snakemake@output[[3]])

WZA_q0.01 = WZA_res[WZA_res$qvalue_gif_adj < 0.01,]
fwrite(WZA_q0.01, snakemake@output[[4]])

WZA_q0.001 = WZA_res[WZA_res$qvalue_gif_adj < 0.001,]
fwrite(WZA_q0.001, snakemake@output[[5]])

# Make sure the data set is a data.table
WZA_res = as.data.table(WZA_res)
# The following filters should be applied separately for each environmental factor
WZA_top0.1 = WZA_res[, .SD[pvalue < quantile(pvalue, 0.1)], by = envfactor]
fwrite(WZA_top0.1, snakemake@output[[6]])

WZA_top0.01 = WZA_res[, .SD[pvalue < quantile(pvalue, 0.01)], by = envfactor]
fwrite(WZA_top0.01, snakemake@output[[7]])

WZA_top0.001 = WZA_res[, .SD[pvalue < quantile(pvalue, 0.001)], by = envfactor]
fwrite(WZA_top0.001, snakemake@output[[8]])

############################
# Manhattan plots
############################
print("Creating Manhattan plots")

# Sort the CHR column properly
# First remove prefixes from contig names
WZA_res$CHR = gsub("^Qrob_H2.3.Sc00000", "", WZA_res$CHR)
WZA_res$CHR = gsub("^Qrob_H2.3.Sc0000", "", WZA_res$CHR)
# Then convert to numeric
WZA_res$CHR = as.integer(WZA_res$CHR)
# Then convert back to factor with properly sorted levels
# If I don't do this, scale_color_manual will not work properly
WZA_res$CHR = factor(WZA_res$CHR, levels=as.character(sort(unique(WZA_res$CHR))))
print(unique(WZA_res$CHR))

# add a column with the calculated -log10(p-value)
WZA_res$score = -log10(WZA_res$gif_cor_pvalue)

p <- ggplot(WZA_res, aes(x=POS, y=score, color=CHR)) + 
  geom_point(alpha=1, size=0.8) +
  geom_point(data=WZA_res[WZA_res$qvalue_gif_adj==TRUE,], aes(x=POS, y=score), color="darkred", size=.8) +
  geom_point(data=WZA_res[WZA_res$qvalue0.001_gif_adj==TRUE,], aes(x=POS, y=score), color="red", size=.8) +
  facet_grid(factor(envfactor, levels=envfactors)~CHR, space = "free_x", scales = "free") +
  scale_color_manual(values = rep(c("blue", "lightblue"), 2000)) + 
  theme_bw() +
  theme(legend.position = "none",
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.spacing.x = unit(0, "lines"),
        panel.border = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(), 
        panel.grid.minor.y = element_blank()) +
  xlab("chromosome")+ ylab("-log10(p-value)")

ggsave(snakemake@output[[1]],width=10,height=60, units="in", dpi=300, limitsize=FALSE)

WZA_res$score = -log10(WZA_res$pvalue)

p <- ggplot(WZA_res, aes(x=POS, y=score, color=CHR)) + 
  geom_point(alpha=1, size=0.8) +
  geom_point(data=WZA_res[WZA_res$qvalue==TRUE,], aes(x=POS, y=score), color="darkred", size=.8) +
  geom_point(data=WZA_res[WZA_res$qvalue0.001==TRUE,], aes(x=POS, y=score), color="red", size=.8) +
  facet_grid(factor(envfactor, levels=envfactors)~CHR, space = "free_x", scales = "free") +
  scale_color_manual(values = rep(c("blue", "lightblue"), 2000)) + 
  theme_bw() +
  theme(legend.position = "none",
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.spacing.x = unit(0, "lines"),
        panel.border = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(), 
        panel.grid.minor.y = element_blank()) +
  xlab("chromosome")+ ylab("-log10(p-value)")

ggsave(snakemake@output[[2]],width=10,height=60, units="in", dpi=300, limitsize=FALSE)

############################
# P-value distribution plots
############################

# Loop over each envfactor
for (env in envfactors) {
  # Set up the file name for the PNG
  png_filename <- paste0(snakemake@params[[1]], "Pvalue_and_QQ_Plots_", env, ".png")
  
  # Open the PNG device
  png(filename = png_filename, width = 1050, height = 750, units = "px")
  
  # Set up the plotting area
  par(mfrow = c(1, 2), mar = c(5, 5, 4, 1))
  
  # Subset data for the current envfactor
  data_subset <- subset(WZA_res, envfactor == env)
  
  # Plot the p-value distribution
  hist(data_subset$pvalue, col = "grey", main = paste("P-value Distribution -", env), xlab = "P-values")
  
  # Plot the QQ-plot
  observed <- -log10(sort(data_subset$pvalue))
  expected <- -log10(ppoints(length(observed)))
  plot(expected, observed, xlab = "Expected -log10(P)", ylab = "Observed -log10(P)", main = paste("QQ-Plot -", env), pch = 19, cex = 1)
  abline(0, 1)
  
  # Add median p-value to the QQ-plot
  median_pvalue <- median(data_subset$pvalue)
  points(-log10(median_pvalue), -log10(median_pvalue), col = "red", pch = 17, cex = 2)
  legend("bottomright", legend = paste("Median P = ", round(median_pvalue, 4)), col = "red", pch = 17, cex = 1)
  
  # Close the PNG device
  dev.off()
}

# Loop over each envfactor
for (env in envfactors) {
  # Set up the file name for the PNG
  png_filename <- paste0(snakemake@params[[1]], "Pvalue_and_QQ_Plots_(GIF_corrected)_", env, ".png")
  
  # Open the PNG device
  png(filename = png_filename, width = 1050, height = 750, units = "px")
  
  # Set up the plotting area
  par(mfrow = c(1, 2), mar = c(5, 5, 4, 1))
  
  # Subset data for the current envfactor
  data_subset <- subset(WZA_res, envfactor == env)
  
  # Plot the p-value distribution
  hist(data_subset$gif_cor_pvalue, col = "grey", main = paste("P-value Distribution -", env), xlab = "P-values")
  
  # Plot the QQ-plot
  observed <- -log10(sort(data_subset$gif_cor_pvalue))
  expected <- -log10(ppoints(length(observed)))
  plot(expected, observed, xlab = "Expected -log10(P)", ylab = "Observed -log10(P)", main = paste("QQ-Plot -", env), pch = 19, cex = 1)
  abline(0, 1)
  
  # Add median p-value to the QQ-plot
  median_pvalue <- median(data_subset$gif_cor_pvalue)
  points(-log10(median_pvalue), -log10(median_pvalue), col = "red", pch = 17, cex = 2)
  legend("bottomright", legend = paste("Median P = ", round(median_pvalue, 4)), col = "red", pch = 17, cex = 1)
  
  # Close the PNG device
  dev.off()
}