
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


# Load the dataframe of all genes with non-zero counts

df <- data.table::fread("data/Allgenes.csv") %>% 
              dplyr::select(Gene.stable.ID)

rownames(df) <- df$Gene.stable.ID
custom_bg <- rownames(df)

# Load the list object containing DEGs results for different contrasts
load("data/Results_DEGs_list.rda")



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#           Fig.1  - DEGs ORA             #
#                                         #
#         Table S2A & S2B - DEGs ORA      #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# intestinal upregulated genes

intestine11.down <- Results_DEGs_list$intestine.IVM11vsCtrl %>% 
  dplyr::filter(log2FoldChange<0) %>% 
  rownames(.)



intestine11.down_ORA_res <- gprofiler2::gost(
  intestine11.down, 
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
intestine11.down_df <- intestine11.down_ORA_res$result %>%
  dplyr::select(source,	term_id,	highlighted, term_name,	term_size,	effective_domain_size,	p_value,	intersection_size,	intersection) %>% 
  dplyr::rename(driver_term = highlighted)

readr::write_csv(intestine11.down_df, file="data/S2A Table.csv")


# Prepare data specifically for plotting
intestine11.down_PlotDATA <- intestine11.down_df %>%
  dplyr::select(source, term_name, p_value, intersection_size, driver_term) %>%
  dplyr::rename(Term = term_name, DEG_count = intersection_size) %>% 
  dplyr::filter(driver_term==TRUE)

intestine11.down_PlotDATA$source <- factor(intestine11.down_PlotDATA$source, levels = unique(intestine11.down_PlotDATA$source))

intestine11.down_PlotDATA$source <- recode(
  intestine11.down_PlotDATA$source,
  `GO:CC` = "Cellular Component",
  `GO:BP` = "Biological Process",
  `GO:MF` = "Molecular Function"
)

# Create and Save Plot
intestine11.down_ORA_plot <- ggplot2::ggplot(intestine11.down_PlotDATA, ggplot2::aes(x = round(-log10(p_value),1), y = Term, group = source)) +
  ggplot2::geom_point(ggplot2::aes(group = source, color=source), size=2) +
  ggplot2::geom_text(ggplot2::aes(label = DEG_count), size = 2, vjust = -1) +
  ggplot2::scale_y_discrete(limits = rev(intestine11.down_PlotDATA$Term)) +
  ggplot2::scale_color_manual(values = c("Cellular Component" = "#00BA38",
                                         "Biological Process" = "#F8766D",
                                         "Molecular Function" = "#619CFF"),
                              name = "Source") +
  ggplot2::scale_size_continuous(name = expression(log[10]("DEG count"))) +
  ggplot2::labs(x = expression(bold(-log[10](italic(p)))), y = "Term") +
  ggplot2::theme_bw() +
  ggplot2::theme(
    axis.title.x = element_text(size = 12),
    axis.text.x = element_text(size = 10),
    axis.title.y = element_text(size = 12, angle=90, face = "bold"),
    axis.text.y = element_text(size = 10)
    
  )

#  intestinal upregulated genes

intestine11.up <- Results_DEGs_list$intestine.IVM11vsCtrl %>% 
  dplyr::filter(log2FoldChange>0) %>% 
  rownames(.)


intestine11.up_ORA_res <- gprofiler2::gost(
  intestine11.up, 
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
intestine11.up_df <- intestine11.up_ORA_res$result %>%
  dplyr::select(source,	term_id,	highlighted, term_name,	term_size,	effective_domain_size,	p_value,	intersection_size,	intersection) %>% 
  dplyr::rename(driver_term = highlighted)

readr::write_csv(intestine11.up_df, file="data/S2B Table.csv")


# Prepare data specifically for plotting
intestine11.up_PlotDATA <- intestine11.up_df %>%
  dplyr::select(source, term_name, p_value, intersection_size, driver_term) %>%
  dplyr::rename(Term = term_name, DEG_count = intersection_size) %>% 
  dplyr::filter(driver_term==TRUE)

intestine11.up_PlotDATA$source <- factor(intestine11.up_PlotDATA$source, levels = unique(intestine11.up_PlotDATA$source))

intestine11.up_PlotDATA$source <- recode(
  intestine11.up_PlotDATA$source,
  `GO:CC` = "Cellular Component",
  `GO:BP` = "Biological Process",
  `GO:MF` = "Molecular Function"
)

# Create and Save Plot
intestine11.up_ORA_plot <- ggplot2::ggplot(intestine11.up_PlotDATA, ggplot2::aes(x = round(-log10(p_value),1), y = Term, group = source)) +
  ggplot2::geom_point(ggplot2::aes(group = source, color=source), size=2) +
  ggplot2::geom_text(ggplot2::aes(label = DEG_count), size = 2, vjust = -1) +
  ggplot2::scale_y_discrete(limits = rev(intestine11.up_PlotDATA$Term)) +
  ggplot2::scale_color_manual(values = c("Cellular Component" = "#00BA38",
                                         "Biological Process" = "#F8766D",
                                         "Molecular Function" = "#619CFF"),
                              name = "Source") +
  ggplot2::scale_size_continuous(name = expression(log[10]("DEG count"))) +
  ggplot2::labs(x = expression(bold(-log[10](italic(p)))), y = "Term") +
  ggplot2::theme_bw() +
  ggplot2::theme(
    axis.title.x = element_text(size = 12),
    axis.text.x = element_text(size = 10),
    axis.title.y = element_text(size = 12, angle=90, face = "bold"),
    axis.text.y = element_text(size = 10)
  )


# create a common theme for all plots

common_theme <- theme(
  legend.position = "bottom",
  legend.justification = c(1, 0), # Anchor point is bottom left
  legend.title = ggplot2::element_blank(),
  legend.text = element_text(size = 10),
  legend.box = "horizontal", # Ensure the legend keys are arranged horizontally
  panel.grid.major.x = element_line(color = "transparent"),
  panel.grid.minor.x = element_line(color = "transparent"),
  plot.margin = margin(t=10, r=10, b=10, l=10)
)


# Add the common theme adjustments to the intestine11.down_ORA_plot
intestine11.down_ORA_plot <- intestine11.down_ORA_plot + common_theme + 
  guides(color = guide_legend(nrow = 1)) # Ensure legend keys in one row


# Add the common theme adjustments to the intestine11.up_ORA_plot
intestine11.up_ORA_plot <- intestine11.up_ORA_plot + common_theme + 
  guides(color = guide_legend(nrow = 1)) # Ensure legend keys in one row



# Modify the plot A to remove the legend
intestine11.down_ORA_plot <- intestine11.down_ORA_plot + theme(legend.position = "none")

# Modify the plot B to keep the legend and move it to the bottom
intestine11.up_ORA_plot <- intestine11.up_ORA_plot + common_theme + 
  guides(color = guide_legend(nrow = 1)) +
  theme(legend.position = "bottom")

# Now, we'll create a legend from plot B and use it separately
legend_b <- cowplot::get_legend(intestine11.up_ORA_plot) 

# Adjust plot B to remove the legend (since we'll add it back separately)
intestine11.up_ORA_plot <- intestine11.up_ORA_plot + theme(legend.position = "none")

# Combine the plots and the legend using plot_grid
intestine_combo <- cowplot::plot_grid(
  cowplot::plot_grid(intestine11.down_ORA_plot, intestine11.up_ORA_plot, ncol = 2, align = "h", rel_widths = c(1, 1.058), labels = c("A)","B)")),
  legend_b, 
  ncol = 1, 
  rel_heights = c(1, 0.08), # Adjust space for the legend
  align = 'h',
  labels = c("","")
)

# Save the combined plot with the central legend for plot B
ggsave("figures/Fig1.tiff", intestine_combo, width = 7.5, height = 3.5, dpi = 600)



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#           Fig.2  - DEGs ORA             #
#                                         #
#         Table S3 - DEGs ORA             #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


intestine9.up <- Results_DEGs_list$intestine.IVM9vsCtrl %>% 
  dplyr::filter(log2FoldChange>0) %>% 
  rownames(.)



intestine9.up_ORA_res <- gprofiler2::gost(
  intestine9.up, 
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
intestine9.up_df <- intestine9.up_ORA_res$result %>%
  dplyr::select(source,	term_id,	highlighted, term_name,	term_size,	effective_domain_size,	p_value,	intersection_size,	intersection) %>% 
  dplyr::rename(driver_term = highlighted)

readr::write_csv(intestine9.up_df, file="data/S3 Table.csv")


# Prepare data specifically for plotting
intestine9.up_PlotDATA <- intestine9.up_df %>%
  dplyr::select(source, term_name, p_value, intersection_size, driver_term) %>%
  dplyr::rename(Term = term_name, DEG_count = intersection_size) %>% 
  dplyr::filter(driver_term==TRUE)

intestine9.up_PlotDATA$source <- factor(intestine9.up_PlotDATA$source, levels = unique(intestine9.up_PlotDATA$source))

intestine9.up_PlotDATA$source <- recode(
  intestine9.up_PlotDATA$source,
  `GO:CC` = "Cellular Component",
  `GO:BP` = "Biological Process",
  `GO:MF` = "Molecular Function"
)

# Create and Save Plot
intestine9.up_ORA_plot <- ggplot2::ggplot(intestine9.up_PlotDATA, ggplot2::aes(x = round(-log10(p_value),1), y = Term, group = source)) +
  ggplot2::geom_point(ggplot2::aes(group = source, color=source), size=2) +
  ggplot2::geom_text(ggplot2::aes(label = DEG_count), size = 2, vjust = -1) +
  ggplot2::scale_y_discrete(limits = rev(intestine9.up_PlotDATA$Term)) +
  ggplot2::scale_color_manual(values = c("Cellular Component" = "#00BA38",
                                         "Biological Process" = "#F8766D",
                                         "Molecular Function" = "#619CFF"),
                              name = "Source") +
  ggplot2::scale_size_continuous(name = expression(log[10]("DEG count"))) +
  ggplot2::labs(x = expression(bold(-log[10](italic(p)))), y = "Term") +
  ggplot2::theme_bw() +
  ggplot2::theme(
    axis.title.x = element_text(size = 12),
    axis.text.x = element_text(size = 10),
    axis.title.y = element_text(size = 12, angle=90, face = "bold"),
    axis.text.y = element_text(size = 10)
    
  )


intestine9.up_ORA_plot <- intestine9.up_ORA_plot + common_theme
ggsave("figures/Fig2.tiff", intestine9.up_ORA_plot, width = 7, height = 3.5, dpi = 600)



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#           Fig.3  - DEGs ORA             #
#                                         #
#         Table S4 - DEGs ORA             #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

anterior9.up <- Results_DEGs_list$anterior.IVM9vsCtrl %>% 
  dplyr::filter(log2FoldChange>0) %>% 
  rownames(.)


anterior9.up_ORA_res <- gprofiler2::gost(
  anterior9.up, 
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
anterior9.up_df <- anterior9.up_ORA_res$result %>%
  dplyr::select(source,	term_id,	highlighted, term_name,	term_size,	effective_domain_size,	p_value,	intersection_size,	intersection) %>% 
  dplyr::rename(driver_term = highlighted)

readr::write_csv(anterior9.up_df, file="data/S4 Table.csv")


# Prepare data specifically for plotting
anterior9.up_PlotDATA <- anterior9.up_df %>%
  dplyr::select(source, term_name, p_value, intersection_size, driver_term) %>%
  dplyr::rename(Term = term_name, DEG_count = intersection_size) %>% 
  dplyr::filter(driver_term==TRUE)

anterior9.up_PlotDATA$source <- factor(anterior9.up_PlotDATA$source, levels = unique(anterior9.up_PlotDATA$source))

anterior9.up_PlotDATA$source <- recode(
  anterior9.up_PlotDATA$source,
  `GO:CC` = "Cellular Component",
  `GO:BP` = "Biological Process",
  `GO:MF` = "Molecular Function"
)

# Create and Save Plot
anterior9.up_ORA_plot <- ggplot2::ggplot(anterior9.up_PlotDATA, ggplot2::aes(x = round(-log10(p_value),1), y = Term, group = source)) +
  ggplot2::geom_point(ggplot2::aes(group = source, color=source), size=2) +
  ggplot2::geom_text(ggplot2::aes(label = DEG_count), size = 2, vjust = -1) +
  ggplot2::scale_y_discrete(limits = rev(anterior9.up_PlotDATA$Term)) +
  ggplot2::scale_color_manual(values = c("Cellular Component" = "#00BA38",
                                         "Biological Process" = "#F8766D",
                                         "Molecular Function" = "#619CFF"),
                              name = "Source") +
  ggplot2::scale_size_continuous(name = expression(log[10]("DEG count"))) +
  ggplot2::labs(x = expression(bold(-log[10](italic(p)))), y = "Term") +
  ggplot2::theme_bw() +
  ggplot2::theme(
    axis.title.x = element_text(size = 12),
    axis.text.x = element_text(size = 10),
    axis.title.y = element_text(size = 12, angle=90, face = "bold"),
    axis.text.y = element_text(size = 10)
    
  )


anterior9.up_ORA_plot <- anterior9.up_ORA_plot + common_theme
ggsave("figures/Fig3.tiff", anterior9.up_ORA_plot, width = 7, height = 3.5, dpi = 600)



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#           Fig.5  - DEGs ORA             #
#                                         #
#         Table S13 - DEGs ORA             #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


Network_custom_bg <- data.table::fread("data/Network/Metrices_ConsesusNetwork.csv") 

Network_custom_bg <- Network_custom_bg$Gene
PgR047_g066NN <- data.table::fread("data/PgR047_g066NN.csv")
PgR047_g066NN <- PgR047_g066NN$Gene

PgR047_g066NN_ORA_res <- gprofiler2::gost(
  PgR047_g066NN, 
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
  custom_bg = Network_custom_bg,
  numeric_ns = "",
  highlight = TRUE,
  sources = c("GO:BP", "GO:MF", "GO:CC"),
  as_short_link = F
)

PgR047_g066NN_df <- PgR047_g066NN_ORA_res$result %>%
  dplyr::select(source,	term_id,	highlighted, term_name,	term_size,	effective_domain_size,	p_value,	intersection_size,	intersection) %>% 
  dplyr::rename(driver_term = highlighted)

readr::write_csv(PgR047_g066NN_df, file="data/S13 Table.csv")


# Prepare data specifically for plotting
PgR047_g066NN_PlotDATA <- PgR047_g066NN_df %>%
  dplyr::select(source, term_name, p_value, intersection_size, driver_term) %>%
  dplyr::rename(Term = term_name, DEG_count = intersection_size) %>% 
  dplyr::filter(driver_term==TRUE)

PgR047_g066NN_PlotDATA$source <- factor(PgR047_g066NN_PlotDATA$source, levels = unique(PgR047_g066NN_PlotDATA$source))

PgR047_g066NN_PlotDATA$source <- recode(
  PgR047_g066NN_PlotDATA$source,
  `GO:CC` = "Cellular Component",
  `GO:BP` = "Biological Process",
  `GO:MF` = "Molecular Function"
)

# Create and Save Plot
PgR047_g066NN_ORA_plot <- ggplot2::ggplot(PgR047_g066NN_PlotDATA, ggplot2::aes(x = round(-log10(p_value),1), y = Term, group = source)) +
  ggplot2::geom_point(ggplot2::aes(group = source, color=source), size=2) +
  ggplot2::geom_text(ggplot2::aes(label = DEG_count), size = 2, vjust = -1) +
  ggplot2::scale_y_discrete(limits = rev(PgR047_g066NN_PlotDATA$Term)) +
  ggplot2::scale_color_manual(values = c("Cellular Component" = "#00BA38",
                                         "Biological Process" = "#F8766D",
                                         "Molecular Function" = "#619CFF"),
                              name = "Source") +
  ggplot2::scale_size_continuous(name = expression(log[10]("DEG count"))) +
  ggplot2::scale_x_continuous(limits = c(3, 4), breaks = seq(3, 4, by = 0.2)) +  # Setting custom x-axis scale with one decimal place
  ggplot2::labs(x = expression(bold(-log[10](italic(p)))), y = "Term") +
  ggplot2::theme_bw() +
  ggplot2::theme(
    axis.title.x = element_text(size = 12),
    axis.text.x = element_text(size = 10),
    axis.title.y = element_text(size = 12, angle=90, face = "bold"),
    axis.text.y = element_text(size = 10)
  )


PgR047_g066NN_ORA_plot <- PgR047_g066NN_ORA_plot + common_theme
ggsave("figures/Fig5.tiff", PgR047_g066NN_ORA_plot, width = 7, height = 3.5, dpi = 600)





# Show session information for debugging and reproducibility
sessionInfo()
