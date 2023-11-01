

# Set the working directory to one level up from the active RStudio document
setwd(glue::glue("{dirname(rstudioapi::getActiveDocumentContext()$path)}/../"))


# Load required packages and custom functions
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
library(dplyr)
library(tidyr)
library(data.table)
library(ggplot2)
library(reshape2)
library(gprofiler2)

source("scripts/functions.R")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#   Fig.S1 - Overlap statistics    #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Loading data
load('data/dds.rda')
load('data/Results_DEGs_list.rda') # list of DEGs (contrats) dataframes
load('data/ListOfModules.rda')
load('data/NetworkGenes.rda')



# Applying hypergeometric test across contrasts and modules
all_genes <- nrow(dds)
results <- list()
contrasts <- Results_DEGs_list


for(contrast in names(contrasts)){
  df <- contrasts[[contrast]] %>% 
    dplyr::mutate(GeneID=rownames(contrasts[[contrast]]))
  df_DEGs <- df$GeneID
  for (mod in names(ListOfModules)){
    mod_df <- ListOfModules[[mod]]
    mod_genes <- mod_df$Gene
    result <- hypergeo_test(df_DEGs, mod_genes, all_genes)
    results[[paste0(contrast, "_vs_", mod)]] <- c(contrast = contrast, module = mod, 
                                                  overlap = result$overlap,
                                                  total_DEGs = result$total_DEGs,
                                                  module_Genes = result$module_Genes,
                                                  All_genes = result$All_genes,
                                                  p_value = result$p_value)
  }
}

# Convert results list to a data frame
results_df <- do.call(rbind, lapply(results, function(x) {
  data.frame(t(unlist(x)))
}))

# Multiple testing correction
results_df$p.adjusted <- p.adjust(results_df$p_value, method = "bonferroni")

# Convert characters to numerics
vars_to_convert <- c("overlap", "total_DEGs", "module_Genes", "All_genes", "p_value")
results_df[vars_to_convert] <- lapply(results_df[vars_to_convert], as.numeric)



results_df$significance <- sapply(results_df$p.adjusted, add_significance)

results_df <- results_df %>% 
  dplyr::mutate(percent_overlap_TotalDEGs = (overlap/total_DEGs)*100,
                percent_overlap_ModuleGenes = (overlap/module_Genes)*100)

# write.csv(results_df,"data/Table S5.csv")

# Reshape the data for plotting
long_df <- reshape2::melt(results_df, 
                          id.vars = c("contrast", "module", "significance", "p.adjusted", "overlap", "total_DEGs"), 
                          measure.vars = c("percent_overlap_TotalDEGs", "percent_overlap_ModuleGenes"), 
                          variable.name = "Percentage_overlap", 
                          value.name = "Count")
                          # id.vars = c("contrast", "module", "significance", "p.adjusted", "total_DEGs"), 
                          # measure.vars = c("module_Genes", "overlap"), 
                          # variable.name = "Gene_set", 
                          # value.name = "Count")


# Remove the "module_" prefix
long_df$module <- gsub("module_", "", long_df$module)

# Convert to numeric
long_df$module <- as.numeric(long_df$module)


long_df$module <- factor(long_df$module, levels = unique(long_df$module))

long_df$contrast <- factor(long_df$contrast, 
                           levels = c("intestine.IVM11vsCtrl", "intestine.IVM9vsCtrl", "anterior.IVM11vsCtrl", "anterior.IVM9vsCtrl"))


long_df$contrast <- sapply(as.character(long_df$contrast), rename_contrast)




long_df$Percentage_overlap <- factor(long_df$Percentage_overlap, levels = c("percent_overlap_TotalDEGs", "percent_overlap_ModuleGenes"))
# Plot the data
Overlap_Plot <- ggplot2::ggplot(long_df, ggplot2::aes(x = module, y = Count + 0.1, fill = Percentage_overlap)) +
                    ggplot2::theme_bw() +
                    ggplot2::geom_bar(stat = "identity", position = "dodge", width = 0.6) +
                    ggplot2::geom_text(data = subset(long_df, Percentage_overlap == "percent_overlap_TotalDEGs"),
                    ggplot2::aes(label=significance, y = Count), vjust=-0.25, position = position_dodge(0.6), size=5) +
                    ggplot2::facet_wrap(~contrast, scales = "free_x", labeller = label_parsed) + 
                    # scale_y_log10() +
                    ggplot2::scale_fill_manual(values = c("#619CFF", "#F8766D"), 
                                      guide = "legend", aesthetics = "fill", 
                                      name = "Overlap categories",
                                      labels = c("DEGs to Module", "Module to DEGs")) +
                    ggplot2::labs(x = "Module", y = "% overlap")+
                    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 0, size = 12),
                          axis.title = ggplot2::element_text(face="bold", size = 15),
                          strip.text.x = ggplot2::element_text(size = 12),
                          legend.title = ggplot2::element_text(face="bold"),
                          legend.text = ggplot2::element_text(size = 10))

# Save the plot
ggplot2::ggsave('figures/Fig.S1.png', Overlap_Plot, width=12, height=7, dpi=600)



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#    Fig.S2 to Fig.S7 - Module ORA       #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# filter out significant modules
significant_modules <- results_df %>%
  dplyr::filter(p.adjusted < 0.05) %>%
  dplyr::select(module) %>%
  dplyr::distinct()

# Now, convert the 'module' column of the significant_modules data frame to a vector
Sign_module_vector <- sort(significant_modules$module)




# Dictionary to match module names to figure names
fig_name_dict <- list(
  module_01 = "Fig.S2",
  module_02 = "Fig.S3",
  module_04 = "Fig.S4",
  module_05 = "Fig.S5",
  module_08 = "Fig.S6",
  module_10 = "Fig.S7"
)

# Dictionary to match module names to table names
table_name_dict <- list(
  module_01 = "Table S6",
  module_02 = "Table S7",
  module_04 = "Table S8",
  module_05 = "Table S9",
  module_08 = "Table S10",
  module_10 = "Table S11"
)


custom_bg <- NetworkGenes$GeneID

# Loop through each module in Sign_module_vector
for (module_name in Sign_module_vector) {
  if (module_name == "module_07"){
    next
  }
  # Extract the module number
  module_num <- sub("module_", "", module_name)
  
  # Get the current module from ListOfModules
  current_module <- ListOfModules[[module_name]]
  
  # Create a list of row names of the current module
  current_module <- data.table::fread("data/PgR160_g006NN.txt")
  rownames(current_module) <- current_module$Gene
  current_query <- list(rownames(current_module))
  
  # Run gprofiler2::gost with the current query
  current_gostres <- gprofiler2::gost(
    current_query,
    organism = "paunivprjna386823",
    ordered_query = FALSE,
    multi_query = FALSE,
    significant = TRUE,
    exclude_iea = FALSE,
    measure_underrepresentation = FALSE,
    evcodes = TRUE,
    user_threshold = 0.05,
    correction_method = "g_SCS",
    domain_scope = "custom",
    custom_bg = custom_bg,
    numeric_ns = "",
    highlight = TRUE,
    sources = c("GO:BP", "GO:MF", "GO:CC"),
    as_short_link = T
  )
  
  
  # Get the figure name from the dictionary
  fig_name <- fig_name_dict[[module_name]]  
  
  # Get the figure name from the dictionary
  table_name <- table_name_dict[[module_name]]  
  
  
  # Extract the result and write to a CSV file
  current_df <- current_gostres$result %>% 
    dplyr::select(source,	term_id,	highlighted, term_name,	term_size,	effective_domain_size,	p_value,	intersection_size,	intersection) %>% 
    dplyr::rename(driver_term = highlighted)
  
  readr::write_csv(current_df, file=paste0("data/", table_name, ".csv"))
  
  # Prepare data for plotting
  current_PlotDATA <- current_df %>%
    dplyr::select(source, term_name, p_value, intersection_size, driver_term) %>%
    dplyr::rename(Term = term_name, DEG_count = intersection_size) %>% 
    dplyr::filter(driver_term==TRUE)
  
  current_PlotDATA$source <- factor(current_PlotDATA$source, levels = unique(current_PlotDATA$source))
  
  current_PlotDATA$source <- recode(
    current_PlotDATA$source,
    `GO:CC` = "Cellular Component",
    `GO:BP` = "Biological Process",
    `GO:MF` = "Molecular Function"
  )
  
  # Create the plot
  ORA_plot <- ggplot2::ggplot(current_PlotDATA, ggplot2::aes(x = -log10(p_value), y = Term, group = source)) +
    ggplot2::geom_point(ggplot2::aes(group = source, color=source), size=3) +
    ggplot2::geom_text(ggplot2::aes(label = DEG_count), size = 4, vjust = -1) +
    ggplot2::scale_y_discrete(limits = rev(current_PlotDATA$Term)) +
    ggplot2::scale_color_manual(values = c("Cellular Component" = "#00BA38",
                                  "Biological Process" = "#F8766D",
                                  "Molecular Function" = "#619CFF"),
                       name = "Source") +
    ggplot2::scale_size_continuous(name = expression(log[10]("DEG count"))) +
    ggplot2::labs(x = expression(-log[10](italic(p))), y = "Term") +
    ggplot2::theme_bw() +
    ggplot2::theme(
      axis.title = ggplot2::element_text(size = 16), # increase axis title size
      axis.text = ggplot2::element_text(size = 13),  # increase axis text size
      legend.title = ggplot2::element_text(size = 16), # increase legend title size
      legend.text = ggplot2::element_text(size = 13), # increase legend text size
      panel.grid.major.x = ggplot2::element_line(color = "transparent"),
      panel.grid.minor.x = ggplot2::element_line(color = "transparent")

    )
    
  # Save the plot to a file
  ggplot2::ggsave(ORA_plot, file=paste0("figures/", fig_name, ".png"), width=12.87, height=7.48)
}



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#           Fig.S8 - Correlation of Network metrics           #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# load data
load("data/cor_mat.spearman.rda")


# Plot the correlation heatmap
Fig.S8  <- pheatmap::pheatmap(cor_mat.spearman$r, 
                            color = viridis::viridis(256, option = "D"), fontsize = 15,
                            scale = "none", limits = c(-0.2, 1), breaks = seq(-0.2, 1, length.out = 256),
                            angle_col = "45", cellwidth = 40, cellheight = 40, display_numbers = T, number_color = 'white',
                            fontsize_number = 15,
                            border_color = "white",
                            cutree_rows = 2, cutree_cols = 2, show_rownames=T, show_colnames=T)
png(filename="figures/Fig.S8.png", width=4000, height=3500, res=600)
Fig.S8 <- grid::grid.draw(Fig.S8)
dev.off()


# Show session information for debugging and reproducibility
sessionInfo()



