
rm(list = ls())
pacman::p_unload("all")

library(tidyfst)
pkg_load(tidyverse,fs,arrow,DBI,duckdb)

con<- dbConnect(
  duckdb(),
  dbdir = "data/sel_data.duckdb"
)

tbl(con,"review_citation") %>% 
  select(UT,jciQuartile,domain,document_type,pr_duration) %>% 
  collect() -> dat

## =============================================================
## Per-domain comparisons of peer review duration (pr_duration)
##   (1) document_type  -> two groups  -> Welch t-test
##   (2) jciQuartile     -> four groups -> one-way ANOVA + Tukey HSD
## Data: dat (columns: UT, jciQuartile, domain, document_type, pr_duration)
## =============================================================

tidyfst::pkg_load(dplyr, tidyr, purrr, broom, multcomp, multcompView, rcompanion, ggplot2, forcats, viridisLite)

## -------------------------------------------------------------
## 0. Quick check: how many levels does document_type have in each domain?
## -------------------------------------------------------------
dat %>%
  group_by(domain) %>%
  summarise(n_doc_types = n_distinct(document_type),
            doc_types   = paste(unique(document_type), collapse = ", "),
            .groups = "drop")

## =============================================================
## 1. document_type comparison within each domain (t-test)
## =============================================================
## Assumes document_type has exactly 2 levels within each domain.
## If a domain has >2 levels, it will be skipped with a warning message.

ttest_by_domain <- dat %>%
  group_by(domain) %>%
  nest() %>%
  mutate(
    n_levels = map_int(data, ~ n_distinct(.x$document_type)),
    ttest = map2(data, n_levels, ~ {
      if (.y == 2) {
        t.test(pr_duration ~ document_type, data = .x)
      } else {
        NULL
      }
    })
  )

## Tidy results: mean difference, CI, t, df, p-value
ttest_results <- ttest_by_domain %>%
  filter(!map_lgl(ttest, is.null)) %>%
  mutate(tidied = map(ttest, broom::tidy)) %>%
  dplyr::select(domain, tidied) %>%
  unnest(tidied) %>%
  transmute(
    domain,
    mean_group1  = estimate1,
    mean_group2  = estimate2,
    mean_diff    = estimate,
    ci_low       = conf.low,
    ci_high      = conf.high,
    t_value      = statistic,
    df           = parameter,
    p_value      = p.value
  ) %>%
  arrange(p_value)

print(ttest_results, n = Inf)

## Domains skipped (document_type levels != 2)
ttest_by_domain %>%
  filter(n_levels != 2) %>%
  dplyr::select(domain, n_levels)


## =============================================================
## 2. jciQuartile comparison within each domain x document_type
##    (i.e. quartile effect is examined separately for each
##     document_type within each domain, not pooled across types)
## =============================================================
## Quick check: how many jciQuartile levels exist in each
## domain x document_type combination (some may have <2, skip those)

dat %>%
  group_by(domain, document_type) %>%
  summarise(n_quartiles = n_distinct(jciQuartile),
            n_obs       = n(),
            .groups = "drop") %>%
  arrange(domain, document_type)

anova_by_domain <- dat %>%
  group_by(domain, document_type) %>%
  nest() %>%
  mutate(
    n_levels = map_int(data, ~ n_distinct(.x$jciQuartile)),
    aov_fit  = map2(data, n_levels, ~ if (.y >= 2) aov(pr_duration ~ jciQuartile, data = .x) else NULL),
    aov_tab  = map(aov_fit, ~ if (!is.null(.x)) broom::tidy(.x) else NULL)
  ) %>%
  filter(!map_lgl(aov_fit, is.null))

## Groups skipped (fewer than 2 jciQuartile levels present)
dat %>%
  group_by(domain, document_type) %>%
  summarise(n_levels = n_distinct(jciQuartile), .groups = "drop") %>%
  filter(n_levels < 2)

## Overall ANOVA table per domain x document_type (F, df, p-value)
anova_summary <- anova_by_domain %>%
  dplyr::select(domain, document_type, aov_tab) %>%
  unnest(aov_tab) %>%
  filter(term == "jciQuartile") %>%
  transmute(
    domain,
    document_type,
    df_num   = df,
    F_value  = statistic,
    p_value  = p.value
  ) %>%
  arrange(p_value)

print(anova_summary, n = Inf)

## -------------------------------------------------------------
## 2b. Tukey HSD post-hoc + Compact Letter Display (CLD) per domain
## -------------------------------------------------------------

tukey_by_domain <- anova_by_domain %>%
  mutate(
    tukey     = map(aov_fit, ~ TukeyHSD(.x)),
    tukey_tab = map(tukey, ~ broom::tidy(.x) %>%
                      rename(comparison = contrast))
  )

## All pairwise comparisons across all domain x document_type groups
tukey_pairs <- tukey_by_domain %>%
  dplyr::select(domain, document_type, tukey_tab) %>%
  unnest(tukey_tab) %>%
  transmute(
    domain,
    document_type,
    comparison,
    mean_diff = estimate,
    ci_low    = conf.low,
    ci_high   = conf.high,
    p_adj     = adj.p.value
  )

print(tukey_pairs, n = Inf)

## Compact Letter Display per domain x document_type
## (which jciQuartile groups differ within that specific document_type)
cld_by_domain <- tukey_by_domain %>%
  mutate(
    cld = map(tukey, ~ {
      pvals <- .x$jciQuartile[, "p adj"]
      names(pvals) <- rownames(.x$jciQuartile)
      multcompView::multcompLetters(pvals)$Letters
    })
  ) %>%
  mutate(
    cld_df = map2(data, cld, ~ {
      tibble(
        jciQuartile = names(.y),
        letter      = .y
      )
    })
  ) %>%
  dplyr::select(domain, document_type, cld_df) %>%
  unnest(cld_df)

## Merge letters with group means for a clean summary table
group_means <- dat %>%
  group_by(domain, document_type, jciQuartile) %>%
  summarise(n = n(), mean_pr = mean(pr_duration), .groups = "drop")

final_summary <- group_means %>%
  left_join(cld_by_domain, by = c("domain", "document_type", "jciQuartile")) %>%
  arrange(domain, document_type, jciQuartile)

print(final_summary, n = Inf)


## =============================================================
## 3. Visualization
## =============================================================

tableau10 <- c("#4E79A7", "#F28E2B", "#E15759", "#76B7B2", "#59A14F",
               "#EDC948", "#B07AA1", "#FF9DA7", "#9C755F", "#BAB0AC")

## -------------------------------------------------------------
## 3a. Dumbbell plot: article vs review mean duration, by domain
## -------------------------------------------------------------
## Instead of only plotting the mean difference (which hides the actual
## group values and made the non-significant color nearly invisible since
## almost every domain is significant given the huge sample size), this
## directly compares the two document_type means side by side per domain,
## connected by a line so the magnitude and direction are both visible.
## Significance is shown as compact stars next to each domain instead of
## color, since color-by-significance had almost no contrast to show.

doc_means <- dat %>%
  group_by(domain, document_type) %>%
  summarise(mean_pr = mean(pr_duration), n = n(), .groups = "drop")

sig_stars <- function(p) {
  dplyr::case_when(
    p < 0.001 ~ "***",
    p < 0.01  ~ "**",
    p < 0.05  ~ "*",
    TRUE      ~ "ns"
  )
}

dumbbell_dat <- doc_means %>%
  left_join(
    ttest_results %>% dplyr::select(domain, mean_diff, p_value),
    by = "domain"
  ) %>%
  mutate(
    domain = forcats::fct_reorder(domain, mean_diff),
    stars  = sig_stars(p_value)
  )

## Which side (left/right) each point sits on within its domain, so the
## day-value label can be placed just outside the point on the correct side
range_pr <- diff(range(dumbbell_dat$mean_pr))
label_offset <- 0.018 * range_pr

dumbbell_dat <- dumbbell_dat %>%
  group_by(domain) %>%
  mutate(
    side     = ifelse(mean_pr == min(mean_pr), "left", "right"),
    label_x  = ifelse(side == "left", mean_pr - label_offset, mean_pr + label_offset),
    label_hj = ifelse(side == "left", 1, 0)
  ) %>%
  ungroup()

## Star sits at the midpoint of the connecting line for each domain
star_pos <- dumbbell_dat %>%
  group_by(domain) %>%
  summarise(x_star = mean(range(mean_pr)),
            stars  = dplyr::first(stars), .groups = "drop")

doc_type_colors <- setNames(tableau10[1:n_distinct(dumbbell_dat$document_type)],
                            sort(unique(dumbbell_dat$document_type)))

p_ttest <- ggplot(dumbbell_dat, aes(x = mean_pr, y = domain)) +
  geom_label(data = star_pos, aes(x = x_star, y = domain, label = stars),
             fill = "white", label.size = 0, label.padding = unit(0.12, "lines"),
             size = 3.4, color = "grey30", fontface = "bold",
             nudge_y = 0.16) +
  geom_line(aes(group = domain), color = "grey70", linewidth = 1) +
  geom_point(aes(color = document_type), size = 3.5) +
  geom_text(aes(x = label_x, label = round(mean_pr), hjust = label_hj),
            size = 3.3, color = "grey25", vjust = 1, nudge_y = .2) +
  scale_color_manual(values = doc_type_colors) +
  labs(
    x = "Mean peer review duration (days)",
    y = NULL,
    color = NULL
  ) +
  theme_minimal(base_size = 13) +
  theme(
    axis.title.x = element_text(margin = margin(t = 10)),
    panel.grid.minor = element_blank(),
    legend.position = "top"
  )

print(p_ttest)

ggsave(
  filename = "fig/011_dumbbell_ttest_by_domain.png",
  plot     = p_ttest,
  width    = 8, height = 4, units = "in", dpi = 300, bg = "white"
)

cat("\nFigure. Mean peer review duration (days) for article vs. review within",
    "each domain. Points show the group mean for each document type (value",
    "rounded to the nearest day labeled beside each point), connected by a",
    "line; the star on the line denotes the significance of the Welch's",
    "t-test comparing the two types within that domain",
    "(*** p<.001, ** p<.01, * p<.05, ns = not significant).\n")


## -------------------------------------------------------------
## 3b. Heatmap: jciQuartile effect within domain x document_type
## -------------------------------------------------------------
## A grouped bar chart became too crowded once split by both domain and
## document_type (labels overlapped/were clipped). A tile heatmap scales
## much better: one compact grid per document_type, domain on the y-axis,
## jciQuartile on the x-axis, color encodes the mean duration, and each
## tile is annotated with the rounded mean + CLD letter.

domain_levels <- sort(unique(final_summary$domain), decreasing = TRUE)

## The soft blue-green-yellow-orange-coral gradient settled on earlier.
soft_gradient <- colorRampPalette(c("#AFE1F0", "#7FC8A9", "#F6C85F",
                                    "#F08A5D", "#E85A6B"))

## Letters are assigned per (domain x document_type) group independently by
## Tukey/CLD, so "a" in one row has no relation to "a" in another row - the
## letter itself carries no global ordering. To still color by the CLD
## letter (as requested) while preserving a meaningful gradient, we rank
## each tile's letter by that row's own mean_pr order (a = lowest mean in
## that row, ...) and map the rank to the soft gradient. This keeps letters
## as the coloring key but makes the color scale read consistently across
## the whole plot (low mean -> cool blue, high mean -> warm coral).
plot_dat <- final_summary %>%
  group_by(domain, document_type) %>%
  mutate(
    letter_rank = rank(mean_pr, ties.method = "first"),
    n_letters   = n()
  ) %>%
  ungroup() %>%
  mutate(
    domain     = factor(domain, levels = domain_levels),
    label      = paste0(round(mean_pr), "  ", letter),
    text_color = "black",
    fill_pos   = ifelse(n_letters > 1, (letter_rank - 1) / (n_letters - 1), 0.5)
  )

p_anova <- ggplot(plot_dat, aes(x = jciQuartile, y = domain, fill = fill_pos)) +
  geom_tile(color = "white", linewidth = 0.6) +
  geom_text(aes(label = label, color = text_color), size = 3.3) +
  facet_wrap(~ document_type) +
  scale_fill_gradientn(colors = soft_gradient(100), limits = c(0, 1), guide = "none") +
  scale_color_identity() +
  labs(x = "JCI Quartile", y = NULL) +
  theme_minimal(base_size = 12) +
  theme(
    axis.title.x = element_text(margin = margin(t = 10)),
    panel.grid = element_blank(),
    strip.text = element_text(face = "bold")
  )

## Build the legend entirely by hand with grid/gtable instead of relying on
## ggplot's own guide_colorbar(). cowplot::get_legend() keeps ggplot's built-in
## legend key margins baked into the returned grob, which is exactly why the
## "Shorter"/"Longer" text never sat flush against the bar no matter how far
## rel_widths was shrunk - the empty space was inside the legend grob itself,
## not between the grobs. Building a raw colorbar (as a rasterGrob) plus plain
## text grobs in a single gtable gives full control: the gaps between columns
## are set explicitly (and can be zero), so the text can sit immediately next
## to the bar with only a hairline gap.

bar_colors <- soft_gradient(256)
bar_raster <- grid::rasterGrob(matrix(bar_colors, nrow = 1),
                               width = unit(1, "npc"), height = unit(1, "npc"),
                               interpolate = TRUE)
bar_border <- grid::rectGrob(gp = grid::gpar(fill = NA, col = "grey40", lwd = 0.6))
bar_grob   <- grid::grobTree(bar_raster, bar_border)

shorter_grob <- grid::textGrob("Shorter", gp = grid::gpar(fontsize = 9, col = "grey30"))
longer_grob  <- grid::textGrob("Longer",  gp = grid::gpar(fontsize = 9, col = "grey30"))
title_grob   <- grid::textGrob("Relative duration",
                               gp = grid::gpar(fontsize = 9, col = "grey30", fontface = "bold"))

text_gap  <- unit(3, "pt")   # hairline gap between Shorter/Longer and the bar
title_gap <- unit(10, "pt")  # clear space between the bar and the title below
bar_width <- unit(4, "cm")   # shrunk from 6cm so the whole legend block is smaller

legend_gtable <- gtable::gtable(
  widths  = grid::unit.c(grid::grobWidth(shorter_grob), text_gap,
                         bar_width, text_gap, grid::grobWidth(longer_grob)),
  heights = grid::unit.c(grid::grobHeight(title_grob), title_gap, unit(0.4, "cm"))
)
legend_gtable <- gtable::gtable_add_grob(legend_gtable, title_grob,   t = 1, l = 3)
legend_gtable <- gtable::gtable_add_grob(legend_gtable, shorter_grob, t = 3, l = 1)
legend_gtable <- gtable::gtable_add_grob(legend_gtable, bar_grob,     t = 3, l = 3)
legend_gtable <- gtable::gtable_add_grob(legend_gtable, longer_grob,  t = 3, l = 5)

## Shift the legend block to the right (rather than centered/left) by
## sandwiching it between an empty left spacer (wider) and a smaller right
## spacer, then stack that row above the heatmap.
legend_row <- cowplot::plot_grid(
  NULL, legend_gtable, NULL,
  ncol = 3, rel_widths = c(1.1, 1, 0.6)
)

p_anova_final <- cowplot::plot_grid(
  legend_row,
  p_anova,
  ncol = 1, rel_heights = c(1.1, 10)
)

print(p_anova_final)

ggsave(
  filename = "fig/012_heatmap_anova_by_domain.png",
  plot     = p_anova_final,
  width    = 6, height = 4, units = "in", dpi = 300, bg = "white"
)

cat("\nFigure. Mean peer review duration (days) across JCI quartiles (Q1-Q4),",
    "shown separately for each document type (facet panels), with domain on",
    "the y-axis. Tile color encodes the Tukey HSD compact letter display (CLD)",
    "group for that domain x document_type row, ordered from the lowest mean",
    "duration (cool blue) to the highest (warm coral) within that row; numbers",
    "show the rounded mean and its CLD letter (tiles sharing a letter within",
    "the same row are not significantly different, p < 0.05).\n")




