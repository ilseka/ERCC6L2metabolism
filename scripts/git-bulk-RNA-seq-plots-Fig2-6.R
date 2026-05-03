# RNASeq dotplots figs 2, 5
library(tidyverse)

# pathway genes 
pathway_genes <- list(
  ferroptosis = c("TFRC", "FTH1", "FTL", "TP53", "SLC7A11", "GPX4", "ACSL4", "ALOX5", "SAT1", "NCOA4", "AIFM2"),
  redox = c("SOD1", "SOD2", "PRDX2", "PRDX3", "CAT", "TXN", "TXNRD1", "PRDX1", "PRDX6", "GSR"),
  niacin = c("NAPRT", "SIRT2", "PARP1", "NAMPT", "SIRT1", "NMNAT1", "PARP2"),
  glutathione = c("GCLC", "GCLM", "GSS", "GSR", "GPX4", "GSTP1", "GSTO1", "SLC7A11", "GLRX", "GGT1"),
  glutamine = c("GLS", "SLC1A5", "IDH2", "GLS2", "GLUL", "SLC38A1", "GLUD1"),
  test = c("NRF1", "NRF2", "ZNRF1", "ZNRF2", "NFE2L2", "NFE2", "KEAP", "CUL3", "MAF"),
  anti_apoptotic = c("BCL2", "MCL1", "BCL2L2", "BCL2A1", "BCL2L1", "BAG3", "XIAP", "BIRC5", "BCL10"),
  ferroptosis_genes = c("NQO1","HMOX1","FTH1","FTL","HERC2","SLC40A1","ABCB6","FECH","PIR","MT1G",
                        "SLC7A11","GCL","GSS","GSR","GPX4","AIFM2","MGST1","ALDH1A1","ALDH3A1","G6PD",
                        "NRF1","NFE2L2","ZNRF1","ZNRF2","GCLC","GCLM","GSTP1","GSTO1","GLRX","GGT1",
                        "TFRC","TP53","ACSL4","ALOX5","SAT1","NCOA4"), 
  glycolysis = c(
    "HK1", "PFKP", "GAPDH", "ENO1", 
    "GPI", "ALDOA", "TPI1", "PGK1", "PGAM1"
  ),
  pyruvate_lactate = c(
    "LDHA", "LDHB", "PDHA1", "PDHB", 
    "PDK1", "PDK2", "PDK3", "PDK4", 
    "MPC1", "MPC2", "DLAT", "DLD"
  ),
  tca = c(
    "OGDH","CS", "ACO2", "IDH2", "IDH3A", "IDH3B", "SUCLG1", "SUCLG2", "FH", 
    "SDHA", "SDHB", "SDHC", "SDHD", 
    "MDH2", "MDH1", "GOT1", "GOT2"
  ),
  glucose_transport = c(
    "SLC2A1", "SLC2A3", "SLC2A4", 
    "SLC2A2", "SLC2A10", "SLC2A8", 
    "SLC5A1", "SLC5A2"
  )
  
)

all_genes <- unique(unlist(pathway_genes))

# Map each gene to its pathway + store the order within each pathway
gene_pathway_map <- enframe(pathway_genes, name = "Pathway", value = "Gene") %>%
  unnest(Gene) %>%
  group_by(Pathway) %>%
  mutate(Gene_order = row_number()) %>%  # keeps your specified order per pathway
  ungroup()


#  Load metadata
fibro_meta <- read_tsv("metadata"
) %>%
  mutate(
    Condition = case_when(
      Erkki == "yes" ~ "ED",
      str_detect(str_to_lower(Condition), "sds") ~ "SDS",
      TRUE ~ "Control"
    ),
    Counts_names = str_replace_all(Counts_names, "_", ".")
  )


# Load and prepare counts
fibro_counts <- read_tsv("counts data",
  col_names = TRUE
)

colnames(fibro_counts)[1] <- "Gene"
colnames(fibro_counts) <- str_replace_all(colnames(fibro_counts), "_", ".")

fibro_counts_filtered <- fibro_counts %>%
  filter(Gene %in% all_genes)

fibro_counts_filtered[,-1] <- lapply(fibro_counts_filtered[,-1], as.numeric)


# join metadata, and filter
fibro_long <- fibro_counts_filtered %>%
  pivot_longer(-Gene, names_to = "Counts_names", values_to = "expression") %>%
  left_join(fibro_meta, by = "Counts_names") %>%
  left_join(gene_pathway_map, by = "Gene") %>%  # brings in Pathway + Gene_order
  mutate(
    Sample_type = "Fibroblast",
    Condition = case_when(
      Erkki == "yes" ~ "ED",
      str_detect(str_to_lower(Condition), "sds") ~ "SDS",
      TRUE ~ "Control"
    )
  ) %>%
  filter(!is.na(Condition), !is.na(Pathway))


# Summarize expression (by Condition only)
fibro_summary <- fibro_long %>%
  group_by(Pathway, Gene, Gene_order, Condition) %>%
  summarize(
    avg_expr = mean(expression, na.rm = TRUE),
    pct_expr = mean(expression > 1, na.rm = TRUE) * 100,
    .groups = "drop"
  ) %>%
  group_by(Gene) %>%
  mutate(zscore = scale(avg_expr)[,1]) %>%
  ungroup()


# Global limits so all pathway plots share the same scales
z_min <- min(fibro_summary$zscore, na.rm = TRUE)
z_max <- max(fibro_summary$zscore, na.rm = TRUE)
pct_max <- max(fibro_summary$pct_expr, na.rm = TRUE)


#  Plot per pathway, preserving gene order from pathway_genes
pathway_plots <- fibro_summary %>%
  group_split(Pathway, .keep = TRUE) %>%
  map(~ {
    df <- .x
    pathway_name <- unique(df$Pathway)
    
    # Make the gene factor levels exactly match your list order (within this pathway)
    df <- df %>%
      arrange(Gene_order) %>%
      mutate(Gene = factor(Gene, levels = unique(Gene)))
    
    
    ggplot(df, aes(x = Condition, y = Gene)) +
      geom_point(aes(size = pct_expr, color = zscore)) +
      scale_color_gradient2(
        low = "blue", mid = "white", high = "red", midpoint = 0,
        limits = c(z_min, z_max), oob = scales::squish
      ) +
      scale_size(range = c(1, 6), limits = c(0, pct_max)) +
      theme_minimal(base_size = 12) +
      labs(
        title = paste("Pathway:", pathway_name),
        x = "Condition", y = "Gene",
        size = "% Expressed", color = "Z-score"
      ) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
  })

# Assign plot names
pathway_names <- fibro_summary %>% distinct(Pathway) %>% pull()
names(pathway_plots) <- pathway_names


#  View plots
pathway_plots$ferroptosis
pathway_plots$redox
pathway_plots$niacin
pathway_plots$glutathione
pathway_plots$glutamine
pathway_plots$anti_apoptotic
pathway_plots$glycolysis
pathway_plots$ferroptosis_genes
pathway_plots$tca
pathway_plots$glucose_transport
pathway_plots$pyruvate_lactate

