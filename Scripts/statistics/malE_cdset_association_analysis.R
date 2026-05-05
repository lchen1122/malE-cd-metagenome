#========================================================
# 0. packages
#========================================================
library(tidyverse)

#========================================================
# 1. path
#========================================================
data_dir <- "F:/其他/skx_cd/"

long_file   <- file.path(data_dir, "family_abundance.long.tsv")
matrix_file <- file.path(data_dir, "family_abundance.matrix.tsv")
reads_file  <- file.path(data_dir, "sample_total_reads.tsv")

out_dir <- file.path(data_dir, "R_results")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

#========================================================
# 2. read data
#========================================================
ab_long <- readr::read_tsv(long_file, show_col_types = FALSE)
ab_mat  <- readr::read_tsv(matrix_file, show_col_types = FALSE)
reads_df <- readr::read_tsv(reads_file, show_col_types = FALSE)

#========================================================
# 3. metadata
#========================================================
meta <- ab_mat %>%
  dplyr::select(sample) %>%
  dplyr::mutate(
    group = dplyr::case_when(
      stringr::str_starts(sample, "control_") ~ "control",
      stringr::str_starts(sample, "pollution_") ~ "pollution",
      TRUE ~ "unknown"
    ),
    group = factor(group, levels = c("control", "pollution"))
  ) %>%
  dplyr::left_join(reads_df, by = "sample")

#========================================================
# 4. helper: robust family matching
#========================================================
pick_family <- function(x, candidates) {
  hit <- candidates[candidates %in% colnames(x)]
  if (length(hit) == 0) return(NULL)
  hit[1]
}

family_names <- setdiff(colnames(ab_mat), "sample")
print(family_names)

malE_col <- pick_family(ab_mat, c("MalE", "malE", "male"))
cadA_col <- pick_family(ab_mat, c("CadA", "cadA"))
czcA_col <- pick_family(ab_mat, c("CzcA", "czcA"))
czcB_col <- pick_family(ab_mat, c("CzcB", "czcB"))
czcC_col <- pick_family(ab_mat, c("CzcC", "czcC"))
mntB_col <- pick_family(ab_mat, c("mntB", "MntB"))
mntC_col <- pick_family(ab_mat, c("mntC", "MntC"))
mntH_col <- pick_family(ab_mat, c("mntH", "MntH"))
smtA_col <- pick_family(ab_mat, c("smtA", "SmtA"))
tolC_col <- pick_family(ab_mat, c("TolC", "tolC"))
zupT_col <- pick_family(ab_mat, c("ZupT", "zupT"))

#========================================================
# 5. build final analysis table
#========================================================
dat <- ab_mat %>%
  dplyr::left_join(meta, by = "sample") %>%
  dplyr::mutate(
    Cd_resistance_set = 0 +
      if (!is.null(cadA_col)) .data[[cadA_col]] else 0 +
      if (!is.null(czcA_col)) .data[[czcA_col]] else 0 +
      if (!is.null(czcB_col)) .data[[czcB_col]] else 0 +
      if (!is.null(czcC_col)) .data[[czcC_col]] else 0 +
      if (!is.null(mntB_col)) .data[[mntB_col]] else 0 +
      if (!is.null(mntC_col)) .data[[mntC_col]] else 0 +
      if (!is.null(mntH_col)) .data[[mntH_col]] else 0 +
      if (!is.null(smtA_col)) .data[[smtA_col]] else 0 +
      if (!is.null(tolC_col)) .data[[tolC_col]] else 0 +
      if (!is.null(zupT_col)) .data[[zupT_col]] else 0,
    
    MalE_abundance = if (!is.null(malE_col)) .data[[malE_col]] else NA_real_
  ) %>%
  dplyr::mutate(
    MalE_log10 = log10(MalE_abundance + 1),
    Cd_resistance_set_log10 = log10(Cd_resistance_set + 1)
  )

print(table(dat$group))
print(nrow(dat))

# 保存分析表
readr::write_tsv(dat, file.path(out_dir, "analysis_table.tsv"))

#========================================================
# 6. Spearman correlations
#========================================================
run_cor <- function(df, xname, yname, group_name) {
  x <- df[[xname]]
  y <- df[[yname]]
  ok <- complete.cases(x, y)
  
  if (sum(ok) < 3) {
    return(tibble(
      group = group_name,
      x = xname,
      y = yname,
      rho = NA_real_,
      p_value = NA_real_,
      n = sum(ok)
    ))
  }
  
  ct <- cor.test(x[ok], y[ok], method = "spearman")
  
  tibble(
    group = group_name,
    x = xname,
    y = yname,
    rho = unname(ct$estimate),
    p_value = ct$p.value,
    n = sum(ok)
  )
}

cor_res <- bind_rows(
  run_cor(dat, "MalE_abundance", "Cd_resistance_set", "all"),
  run_cor(dplyr::filter(dat, group == "control"), "MalE_abundance", "Cd_resistance_set", "control"),
  run_cor(dplyr::filter(dat, group == "pollution"), "MalE_abundance", "Cd_resistance_set", "pollution")
)

readr::write_tsv(cor_res, file.path(out_dir, "spearman_correlations.tsv"))
print(cor_res)

#========================================================
# 7. interaction model
#========================================================
fit_cdset <- lm(Cd_resistance_set_log10 ~ MalE_log10 * group, data = dat)

sink(file.path(out_dir, "interaction_model.txt"))
cat("Model: Cd_resistance_set_log10 ~ MalE_log10 * group\n\n")
print(summary(fit_cdset))
sink()

summary(fit_cdset)

#========================================================
# 8.  p 值格式化函数
#========================================================
fmt_p <- function(p) {
  dplyr::case_when(
    is.na(p) ~ "p = NA",
    p < 0.001 ~ "p < 0.001",
    TRUE ~ paste0("p = ", sprintf("%.3f", p))
  )
}

#========================================================
# 9. 分组计算 Spearman rho 和 p
#========================================================
cor_lab <- dat %>%
  dplyr::group_by(group) %>%
  dplyr::summarise(
    rho = suppressWarnings(
      cor(MalE_abundance, Cd_resistance_set, method = "spearman", use = "complete.obs")
    ),
    p = suppressWarnings(
      cor.test(MalE_abundance, Cd_resistance_set, method = "spearman")$p.value
    ),
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    label = dplyr::case_when(
      group == "control" ~ paste0("control: \u03C1 = ", sprintf("%.2f", rho), ", ", fmt_p(p)),
      group == "pollution" ~ paste0("pollution: \u03C1 = ", sprintf("%.2f", rho), ", ", fmt_p(p)),
      TRUE ~ paste0(group, ": \u03C1 = ", sprintf("%.2f", rho), ", ", fmt_p(p))
    )
  )

#========================================================
# 10. 提取 interaction p
#========================================================
fit_cdset <- lm(Cd_resistance_set_log10 ~ MalE_log10 * group, data = dat)

coef_tab <- summary(fit_cdset)$coefficients
p_int <- coef_tab["MalE_log10:grouppollution", "Pr(>|t|)"]

int_lab <- paste0("interaction: ", fmt_p(p_int))

#========================================================
# 11. 指定左上角文字位置
#========================================================
x_min  <- min(dat$MalE_log10, na.rm = TRUE)
x_max  <- max(dat$MalE_log10, na.rm = TRUE)
y_min  <- min(dat$Cd_resistance_set_log10, na.rm = TRUE)
y_max  <- max(dat$Cd_resistance_set_log10, na.rm = TRUE)

x_pos <- x_min + 0.03 * (x_max - x_min)
y1 <- y_max - 0.04 * (y_max - y_min)
y2 <- y_max - 0.11 * (y_max - y_min)
y3 <- y_max - 0.18 * (y_max - y_min)

cor_lab <- cor_lab %>%
  dplyr::arrange(group) %>%
  dplyr::mutate(
    x = x_pos,
    y = c(y1, y2)[seq_len(n())]
  )

int_df <- data.frame(
  x = x_pos,
  y = y3,
  label = int_lab
)

#========================================================
# 12. 作图：散点 + 线 + 左上角标注
#========================================================
p_final_anno <- ggplot(dat, aes(x = MalE_log10, y = Cd_resistance_set_log10, color = group)) +
  geom_point(size = 3, alpha = 0.85) +
  geom_smooth(method = "lm", se = FALSE, linewidth = 0.9) +
  geom_text(
    data = cor_lab,
    aes(x = x, y = y, label = label, color = group),
    inherit.aes = FALSE,
    hjust = 0,
    size = 4.2,
    show.legend = FALSE
  ) +
  geom_text(
    data = int_df,
    aes(x = x, y = y, label = label),
    inherit.aes = FALSE,
    hjust = 0,
    color = "black",
    size = 4.2
  ) +
  labs(
    x = "log10(MalE abundance + 1)",
    y = "log10(Curated Cd-resistance gene set abundance + 1)",
    title = "MalE vs curated Cd-resistance gene set"
  ) +
  theme_bw(base_size = 12)

print(p_final_anno)
