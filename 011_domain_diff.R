rm(list = ls())

library(tidyfst)
pkg_load(tidyverse,fs,arrow,DBI,duckdb)

con<- dbConnect(
  duckdb(),
  dbdir = "data/sel_data.duckdb"
)

tbl(con,"review_citation") %>% 
  select(UT,domain,pr_duration) %>% 
  collect() -> dat

# ============================================================
# 单因素方差分析：不同 domain 间 pr_duration 的差异
# 数据框 dat 结构：UT | domain | pr_duration
# ============================================================

tidyfst::pkg_load(
  "dplyr",
  "ggplot2",
  "multcomp",       # Tukey HSD 及 CLD 字母提取
  "multcompView",   # multcompLetters
  "rcompanion"
)

# ---------------- 1. 单因素 ANOVA ----------------
dat$domain <- as.factor(dat$domain)

aov_fit <- aov(pr_duration ~ domain, data = dat)
aov_summary <- summary(aov_fit)
print(aov_summary)

# 提取 F 值、自由度、p 值，用于图上标注
f_val   <- aov_summary[[1]]["domain", "F value"]
df1     <- aov_summary[[1]]["domain", "Df"]
df2     <- aov_summary[[1]]["Residuals", "Df"]
p_val   <- aov_summary[[1]]["domain", "Pr(>F)"]
p_label <- ifelse(p_val < 0.001, "p < 0.001", paste0("p = ", signif(p_val, 3)))

anova_label <- sprintf("One-way ANOVA: F(%d, %d) = %.2f, %s",
                       df1, df2, f_val, p_label)

# （可选）方差齐性检验，若不齐考虑 Welch ANOVA + Games-Howell
# car::leveneTest(pr_duration ~ domain, data = dat)

# ---------------- 2. Tukey HSD 事后比较 + CLD 字母 ----------------
tukey_res <- TukeyHSD(aov_fit)

cld <- multcompLetters4(aov_fit, tukey_res)
cld_df <- data.frame(
  domain = names(cld$domain$Letters),
  letter = cld$domain$Letters,
  row.names = NULL
)

# ---------------- 3. 汇总统计表（均值、SD、SE、n、95%CI） ----------------
summary_tbl <- dat %>%
  group_by(domain) %>%
  summarise(
    n      = n(),
    mean   = mean(pr_duration, na.rm = TRUE),
    sd     = sd(pr_duration, na.rm = TRUE),
    se     = sd / sqrt(n),
    ci_low  = mean - qt(0.975, n - 1) * se,
    ci_high = mean + qt(0.975, n - 1) * se,
    .groups = "drop"
  ) %>%
  left_join(cld_df, by = "domain") %>%
  arrange(mean)

print(summary_tbl)

# 固定 domain 的因子顺序（按均值从低到高），用于绘图
summary_tbl$domain <- factor(summary_tbl$domain, levels = summary_tbl$domain)

# ---- 4. 可视化：柱状图 + 误差线 + 数值标签 + CLD 字母 ----

# 离散型调色板：Set2 风格的柔和但清晰配色（10 类，色相区分度高）
domain_colors <- c(
  "#4E79A7", "#F28E2B", "#E15759", "#76B7B2", "#59A14F",
  "#EDC948", "#B07AA1", "#FF9DA7", "#9C755F", "#BAB0AC"
)

summary_tbl <- summary_tbl %>%
  mutate(label = paste0(floor(mean + 0.5), " ",letter)) # 四舍五入

fig_caption <- "Figure 1. Mean peer review duration by academic domain. Labels show the mean value followed by the Tukey HSD post-hoc group letter; domains sharing the same letter do not differ significantly (p > 0.05)."

p <- ggplot(summary_tbl, aes(x = domain, y = mean, fill = domain)) +
  geom_col(width = 0.65, color = "white", linewidth = 0.3) +
  geom_text(aes(label = label),
            hjust = -0.15, size = 4, fontface = "bold", color = "grey20") +
  scale_fill_manual(values = domain_colors) +
  coord_flip() +
  labs(
    x = NULL,
    y = "Mean peer review duration (days)"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    legend.position = "none",
    axis.text.y = element_text(size = 11),
    axis.title.x = element_text(margin = margin(t = 15)),
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_blank()
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15)))

print(p)
cat("\n", fig_caption, "\n")

ggsave("fig/01_domain_anova_barplot.png", p, width = 6, height = 4, dpi = 300)

# info

# > cat("\n", fig_caption, "\n")
# 
# Figure 1. Mean peer review duration by academic domain. Labels show the mean value followed by the Tukey HSD post-hoc group letter; domains sharing the same letter do not differ significantly (p > 0.05). 
# > print(aov_summary)
# Df    Sum Sq   Mean Sq F value Pr(>F)    
# domain            9 1.509e+09 167653075   18264 <2e-16 ***
#   Residuals   5646211 5.183e+10      9179                   
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

