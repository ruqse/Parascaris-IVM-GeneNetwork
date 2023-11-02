
# Set the working directory to one level up from the active RStudio document
setwd(glue::glue("{dirname(rstudioapi::getActiveDocumentContext()$path)}/../"))


# Load required packages and custom functions
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

library(dplyr)
library(tidyr)
library(readr)
library(ggplot2)
library(gprofiler2)

source("scripts/functions.R")


# Load the dataframe of all genes with non-NA p-values
load("data/all_resDFnoNA.rda")
custom_bg <- rownames(all_resDFnoNA)

# Load the list object containing DEGs results for different contrasts
load("data/Results_DEGs_list.rda")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#    Fig.1 to Fig.3 - DEGs ORA            #
#                                         #
#    Table S1 to Table S3 - DEGs ORA      #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


# Define Figure and Table Name Mappings
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Define dictionaries to map contrast names to figure and table names
fig_name_dict <- list(
  intestine.IVM11vsCtrl = "Fig.1",
  intestine.IVM9vsCtrl = "Fig.2",
  anterior.IVM9vsCtrl = "Fig.3"
)

table_name_dict <- list(
  intestine.IVM11vsCtrl = "Table S2",
  intestine.IVM9vsCtrl = "Table S3",
  anterior.IVM9vsCtrl = "Table S4"
)

# Loop through each Contrast for ORA Analysis and Plotting
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
for (contrast_name in names(Results_DEGs_list)) {
  # Skip a specific contrast
  if (contrast_name == "anterior.IVM11vsCtrl"){
    next
  }
  
  # Extract data for the current contrast
  current_contrast <- Results_DEGs_list[[contrast_name]]
  current_query <- list(rownames(current_contrast))
  
  # Run Over-Representation Analysis (ORA) using gprofiler2

  current_ORA_res <- gprofiler2::gost(
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
    as_short_link = F
  )
  
  # Get the corresponding figure and table names
  fig_name <- fig_name_dict[[contrast_name]]
  table_name <- table_name_dict[[contrast_name]]
  
  # Prepare Data for CSV and Plot
  # Extract relevant columns and write to a CSV
  current_df <- current_ORA_res$result %>%
    dplyr::select(source,	term_id,	highlighted, term_name,	term_size,	effective_domain_size,	p_value,	intersection_size,	intersection) %>% 
    dplyr::rename(driver_term = highlighted)
  
  readr::write_csv(current_df, file=paste0("data/", table_name, ".csv"))
  
  
  # Prepare data specifically for plotting
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
  
  # Create and Save Plot
  ORA_plot <- ggplot2::ggplot(current_PlotDATA, ggplot2::aes(x = -log10(p_value), y = Term, group = source)) +
    ggplot2::geom_point(ggplot2::aes(group = source, color=source), size=5) +
    ggplot2::geom_text(ggplot2::aes(label = DEG_count), size = 5, vjust = -1) +
    ggplot2::scale_y_discrete(limits = rev(current_PlotDATA$Term)) +
    ggplot2::scale_color_manual(values = c("Cellular Component" = "#00BA38",
                                           "Biological Process" = "#F8766D",
                                           "Molecular Function" = "#619CFF"),
                                name = "Source") +
    ggplot2::scale_size_continuous(name = expression(log[10]("DEG count"))) +
    ggplot2::labs(x = expression(-log[10](italic(p))), y = "Term") +
    ggplot2::theme_bw() +
    ggplot2::theme(
      axis.title = ggplot2::element_text(size = 18), # increase axis title size
      axis.text = ggplot2::element_text(size = 15),  # increase axis text size
      legend.title = ggplot2::element_text(size = 18), # increase legend title size
      legend.text = ggplot2::element_text(size = 15), # increase legend text size
      panel.grid.major.x = ggplot2::element_line(color = "transparent"),
      panel.grid.minor.x = ggplot2::element_line(color = "transparent")
      
    )
  ggplot2::ggsave(ORA_plot, file=paste0("figures/", fig_name, ".png"), width=12.87, height=7.48)
}

# Show session information for debugging and reproducibility
sessionInfo()

