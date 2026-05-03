
# Bulk RNA-seq: Cell-cycle phase scores Fig 1
# Plots: S score, G2M score, Cycling index, G1/G0-like
# Colors: Control / ERCC6L2 / SDS


library(tidyverse)
library(stringr)

# --------- Colors ----------
pal <- c(
  "Control" = "darkolivegreen",
  "ERCC6L2" = "orange2",
  "SDS"     = "plum3"
)

# S-phase gene list 
s_phase_genes <- c(
  "MCM2","MCM3","MCM4","MCM5","MCM6","MCM7",
  "CDC6","CDT1","CDC45","ORC1","ORC2","ORC5",
  "GINS1","GINS2","GINS3","GINS4",
  "PCNA","RFC2","RFC3","RFC4","RFC5",
  "RPA1","RPA2","RPA3",
  "POLA1","POLA2","POLD1","POLD2","POLE",
  "PRIM1","PRIM2",
  "RRM1","RRM2","TYMS","UNG","FEN1","EXO1",
  "CLSPN","TIPIN","TIMELESS","RAD51","RAD51AP1",
  "ATR","CHEK1","CHEK2",
  "SLBP","NASP","ASF1B","HIST1H1C","HIST1H2AC","HIST1H4C",
  "GMNN","CDCA7","DTL","UHRF1","HELLS",
  "CCNE1","CCNE2"
)

# G2M gene list 
g2m_genes <- c(
  "HMGB2","CDK1","NUSAP1","UBE2C","BIRC5","TPX2","TOP2A","NDC80","CKS2",
  "MKI67","CENPF","TACC3","SMC4","CCNB2","AURKB","BUB1","CDC20",
  "KIF11","ANLN","AURKA","CENPE"
)


# Z-score per gene across samples 
expr_z <- t(scale(t(expr_use)))

# Genes present
s_used   <- intersect(s_phase_genes, rownames(expr_z))
g2m_used <- intersect(g2m_genes, rownames(expr_z))

message("S genes used: ", length(s_used), " / ", length(s_phase_genes))
message("G2M genes used: ", length(g2m_used), " / ", length(g2m_genes))

#  Compute scores
scores_df <- tibble(
  Sample = colnames(expr_z),
  S_score   = colMeans(expr_z[s_used, , drop = FALSE], na.rm = TRUE),
  G2M_score = colMeans(expr_z[g2m_used, , drop = FALSE], na.rm = TRUE)
) %>%
  mutate(
    Cycling_index = (S_score + G2M_score) / 2,
    G1_G0_like = -Cycling_index
  ) %>%
  left_join(meta_use, by = c("Sample" = "Counts_names")) %>%
  filter(!is.na(Condition))

# Long format for plotting 
scores_long <- scores_df %>%
  select(Sample, Condition, S_score, G2M_score, Cycling_index, G1_G0_like) %>%
  pivot_longer(
    cols = c(S_score, G2M_score, Cycling_index, G1_G0_like),
    names_to = "Metric",
    values_to = "Score"
  ) %>%
  mutate(
    Metric = factor(
      Metric,
      levels = c("S_score", "G2M_score", "Cycling_index", "G1_G0_like"),
      labels = c("S phase", "G2/M phase", "Cycling index", "G1/G0-like")
    )
  )

# Plot: show each sample + mean ± 95% CI 
p_scores <- ggplot(scores_long, aes(x = Condition, y = Score, color = Condition)) +
  geom_jitter(width = 0.12, height = 0, size = 2.2, alpha = 0.75, show.legend = FALSE) +
  stat_summary(fun.data = mean_cl_normal, geom = "pointrange", size = 0.9, show.legend = FALSE) +
  facet_wrap(~ Metric, ncol = 2, scales = "free_y") +
  scale_color_manual(values = pal) +
  labs(
    title = "Cell-cycle phase scores (bulk RNA-seq; z-scored gene averages)",
    x = NULL,
    y = "Score (gene-wise z, averaged)"
  ) +
  theme_bw(base_size = 16) +
  theme(
    plot.title = element_text(face = "bold"),
    strip.text = element_text(face = "bold"),
    axis.text.x = element_text(angle = 20, hjust = 1)
  )

print(p_scores)


# Group comparisons (Welch t-test)

compare_scores <- function(df, score_col) {
  df %>%
    select(Condition, all_of(score_col)) %>%
    filter(Condition %in% c("Control", "ERCC6L2", "SDS")) %>%
    group_modify(~{
      tibble(
        comparison = c("ERCC6L2 vs Control", "SDS vs Control"),
        p_value = c(
          t.test(.x[[score_col]] ~ .x$Condition,
                 subset = .x$Condition %in% c("Control","ERCC6L2"))$p.value,
          t.test(.x[[score_col]] ~ .x$Condition,
                 subset = .x$Condition %in% c("Control","SDS"))$p.value
        )
      )
    })
}

stats_S   <- compare_scores(scores_df, "S_score")   %>% mutate(Metric = "S phase")
stats_G2M <- compare_scores(scores_df, "G2M_score") %>% mutate(Metric = "G2/M phase")
stats_Cyc <- compare_scores(scores_df, "Cycling_index") %>% mutate(Metric = "Cycling index")

stats_all <- bind_rows(stats_S, stats_G2M, stats_Cyc)

print(stats_all)

