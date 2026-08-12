library(tidyfst)
pkg_load(tidyverse,fs,arrow,DBI,duckdb)

con <- dbConnect(
  duckdb(),
  dbdir = "data/sel_data.duckdb"
)

tbl(con,"review_citation") %>% 
  transmute(domain, cnp = citation_normalized_percentile, pr_duration) %>% 
  collect() -> dat

# ============================================================================
# 审稿周期与引用影响力：两断点分段线性回归
# 总体分析 + 分领域比较
#
# 复现性策略：使用缓存机制存储已估计的断点
# 第一次运行：估计断点 → 缓存存储
# 后续运行：直接读取缓存 → 完全可复现
# ============================================================================

library(tidyfst)

tidyfst::pkg_load(
  dplyr, tidyr, purrr, tibble, stringr,
  ggplot2, scales, patchwork,
  segmented, qs2, gt, conflicted
)

conflicted::conflict_prefer("filter", "dplyr")
conflicted::conflict_prefer("select", "dplyr")
conflicted::conflict_prefer("lag", "dplyr")

# ---- 1. 参数集中管理 ------------------------------------------------------

PARS <- list(
  duration_lo = 1,
  duration_hi = 730,
  
  psi_subsample = 150000L,
  psi_init = c(21, 102),
  
  min_domain_n = 20000L,
  min_segment_n = 2000L,
  
  facet_nrow = 2L,
  facet_ncol = 5L,
  
  slope_per_days = 10,
  
  fig1_w = 9, fig1_h = 7,
  fig2_w = 14, fig2_h = 8,
  fig_dpi = 400
)

PARS$n_facet <- PARS$facet_nrow * PARS$facet_ncol

# ---- 2. 输出目录 ----------------------------------------------------------

if (!dir.exists("data")) {
  stop("未找到 data 目录。请先在工作目录下手动创建 data。", call. = FALSE)
}

out_dir <- file.path("data", "segmented_by_domain")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

cache_dir <- file.path(out_dir, ".cache")
if (!dir.exists(cache_dir)) dir.create(cache_dir, recursive = TRUE)

# ---- 3. 配色方案 ----------------------------------------------------------

SEG_COLS <- c(
  "Segment 1" = "#1B6CA8",
  "Segment 2" = "#2E9E5B",
  "Segment 3" = "#D1495B"
)

BRK_COLS <- c(
  "Break 1" = "#7B4FA8",
  "Break 2" = "#E08214"
)

GREY_PT <- "grey55"

# ---- 4. 读入数据 ----------------------------------------------------------

stopifnot(
  exists("dat"),
  all(c("domain", "cnp", "pr_duration") %in% names(dat))
)

analysis_dat <- dat |>
  filter(
    !is.na(cnp),
    !is.na(pr_duration),
    !is.na(domain),
    pr_duration >= PARS$duration_lo,
    pr_duration <= PARS$duration_hi
  ) |>
  mutate(
    domain = stringr::str_squish(as.character(domain)),
    row_id = dplyr::row_number()
  )

message(
  "分析样本量：", format(nrow(analysis_dat), big.mark = ","),
  "（周期限定在 ", PARS$duration_lo, "-", PARS$duration_hi, " 天）"
)

# ---- 5. 核心函数 ----------------------------------------------------------

take_even_subsample <- function(d, n_target) {
  n <- nrow(d)
  if (n <= n_target) return(d)
  idx <- round(seq(1, n, length.out = n_target))
  d[unique(idx), , drop = FALSE]
}

# 5.2 估计断点，支持缓存机制
estimate_breaks <- function(d, cache_file = NULL) {
  # 检查缓存：如果缓存存在，直接读取
  if (!is.null(cache_file) && file.exists(cache_file)) {
    cached <- readRDS(cache_file)
    return(cached)
  }
  
  # 没有缓存时计算
  sub <- take_even_subsample(d, PARS$psi_subsample)
  base_fit <- lm(cnp ~ pr_duration, data = sub)
  
  seg_fit <- tryCatch(
    segmented::segmented(
      base_fit,
      seg.Z = ~ pr_duration,
      psi   = PARS$psi_init,
      control = segmented::seg.control(
        it.max = 100,
        tol = 1e-7,
        display = FALSE
      )
    ),
    error = function(e) NULL
  )
  
  if (is.null(seg_fit) || is.null(seg_fit$psi)) return(NULL)
  
  psi <- sort(as.numeric(seg_fit$psi[, "Est."]))
  if (length(psi) != 2 || any(!is.finite(psi))) return(NULL)
  
  # 保存到缓存
  if (!is.null(cache_file)) {
    saveRDS(psi, cache_file)
  }
  
  psi
}

# 5.3 固定断点后在全量数据上拟合，并导出三段斜率
fit_fixed_breaks <- function(d, psi) {
  d2 <- d |>
    mutate(
      h1 = pmax(pr_duration - psi[1], 0),
      h2 = pmax(pr_duration - psi[2], 0)
    )
  
  fit <- lm(cnp ~ pr_duration + h1 + h2, data = d2)
  
  b  <- coef(fit)
  V  <- vcov(fit)
  df <- df.residual(fit)
  
  L <- rbind(
    "Segment 1" = c(0, 1, 0, 0),
    "Segment 2" = c(0, 1, 1, 0),
    "Segment 3" = c(0, 1, 1, 1)
  )
  
  slope    <- as.numeric(L %*% b)
  slope_se <- sqrt(diag(L %*% V %*% t(L)))
  t_stat   <- slope / slope_se
  p_val    <- 2 * pt(abs(t_stat), df = df, lower.tail = FALSE)
  
  n_seg <- c(
    sum(d2$pr_duration <= psi[1]),
    sum(d2$pr_duration > psi[1] & d2$pr_duration <= psi[2]),
    sum(d2$pr_duration > psi[2])
  )
  
  k <- PARS$slope_per_days
  
  slopes <- tibble::tibble(
    segment     = rownames(L),
    range_lo    = c(min(d2$pr_duration), psi[1], psi[2]),
    range_hi    = c(psi[1], psi[2], max(d2$pr_duration)),
    slope_daily = slope,
    slope_per10 = slope * k,
    se_per10    = slope_se * k,
    t_stat      = t_stat,
    p_value     = p_val,
    stars       = sig_stars(p_val),
    n_obs       = n_seg,
    stable      = ifelse(n_seg >= PARS$min_segment_n, "yes", "no")
  )
  
  list(fit = fit, psi = psi, slopes = slopes, data_n = nrow(d2))
}

sig_stars <- function(p) {
  dplyr::case_when(
    is.na(p)   ~ "",
    p < 0.001  ~ "***",
    p < 0.01   ~ "**",
    p < 0.05   ~ "*",
    TRUE       ~ "n.s."
  )
}

effect_size <- function(fit, d, psi) {
  grid <- tibble::tibble(
    pr_duration = seq(min(d$pr_duration), max(d$pr_duration), length.out = 500)
  ) |>
    mutate(
      h1  = pmax(pr_duration - psi[1], 0),
      h2  = pmax(pr_duration - psi[2], 0),
      fit = predict(fit, newdata = pick(everything()))
    )
  
  q <- quantile(d$pr_duration, c(0.05, 0.95), na.rm = TRUE)
  inner <- grid |> filter(pr_duration >= q[1], pr_duration <= q[2])
  
  sd_y <- sd(d$cnp, na.rm = TRUE)
  
  tibble::tibble(
    sd_cnp         = sd_y,
    span_full      = max(grid$fit) - min(grid$fit),
    span_full_sd   = (max(grid$fit) - min(grid$fit)) / sd_y,
    span_inner     = max(inner$fit) - min(inner$fit),
    span_inner_sd  = (max(inner$fit) - min(inner$fit)) / sd_y,
    r_squared      = summary(fit)$r.squared
  )
}

make_curve <- function(fit, d, psi, label = NA_character_) {
  tibble::tibble(
    pr_duration = seq(min(d$pr_duration), max(d$pr_duration), length.out = 600)
  ) |>
    mutate(
      h1 = pmax(pr_duration - psi[1], 0),
      h2 = pmax(pr_duration - psi[2], 0)
    ) |>
    mutate(fit = predict(fit, newdata = pick(everything()))) |>
    mutate(
      segment = dplyr::case_when(
        pr_duration <= psi[1] ~ "Segment 1",
        pr_duration <= psi[2] ~ "Segment 2",
        TRUE                  ~ "Segment 3"
      ),
      group_label = label
    ) |>
    select(pr_duration, fit, segment, group_label)
}

bin_means <- function(d, n_bin = 45, min_n = 500, label = NA_character_) {
  d |>
    mutate(bin = cut(pr_duration, breaks = n_bin)) |>
    group_by(bin) |>
    summarise(
      x = mean(pr_duration, na.rm = TRUE),
      y = mean(cnp, na.rm = TRUE),
      n = dplyr::n(),
      .groups = "drop"
    ) |>
    filter(n >= min_n) |>
    mutate(group_label = label) |>
    select(x, y, n, group_label)
}

# ---- 6. 总体分析 ----------------------------------------------------------

message("\n=== 总体分析 ===")

# 使用缓存文件存储总体断点
cache_overall <- file.path(cache_dir, "breaks_overall.rds")
psi_all <- estimate_breaks(analysis_dat, cache_file = cache_overall)

if (is.null(psi_all)) {
  stop("总体断点估计失败，请调整 psi_init 初值。", call. = FALSE)
}

message("断点估计：day ", paste(round(psi_all, 1), collapse = " / day "))

res_all   <- fit_fixed_breaks(analysis_dat, psi_all)
eff_all   <- effect_size(res_all$fit, analysis_dat, psi_all)
curve_all <- make_curve(res_all$fit, analysis_dat, psi_all, "Overall")
bins_all  <- bin_means(analysis_dat, n_bin = 60, min_n = 500, "Overall")

print(res_all$slopes)
print(eff_all)

# ---- 7. 分领域分析 --------------------------------------------------------

message("\n=== 分领域分析 ===")

domain_size <- analysis_dat |>
  count(domain, name = "n_obs") |>
  arrange(desc(n_obs))

domain_keep <- domain_size |>
  filter(n_obs >= PARS$min_domain_n) |>
  slice_head(n = PARS$n_facet) |>
  pull(domain)

message("进入分面图的领域（", length(domain_keep), " 个）：")
message(paste(" -", domain_keep, collapse = "\n"))

if (nrow(domain_size) > length(domain_keep)) {
  dropped <- setdiff(domain_size$domain, domain_keep)
  message("未纳入分面图的领域：", paste(dropped, collapse = ", "))
}

fit_one_domain <- function(dm) {
  d <- analysis_dat |> filter(domain == dm)
  
  # 每个领域的断点使用独立的缓存文件
  cache_domain <- file.path(cache_dir, paste0("breaks_", chartr(" ", "_", dm), ".rds"))
  psi <- estimate_breaks(d, cache_file = cache_domain)
  
  if (is.null(psi)) {
    message("  [", dm, "] 断点估计失败，改用总体断点固定拟合。")
    psi <- psi_all
    psi_source <- "overall"
  } else {
    psi_source <- "own"
  }
  
  res <- fit_fixed_breaks(d, psi)
  eff <- effect_size(res$fit, d, psi)
  
  list(
    domain     = dm,
    psi        = psi,
    psi_source = psi_source,
    slopes     = res$slopes |> mutate(domain = dm, psi_source = psi_source),
    effects    = eff |> mutate(domain = dm),
    curve      = make_curve(res$fit, d, psi, dm),
    bins       = bin_means(d, n_bin = 40, min_n = 200, dm),
    breaks_tbl = tibble::tibble(
      domain    = dm,
      break_id  = c("Break 1", "Break 2"),
      break_day = psi
    ),
    fit        = res$fit
  )
}

domain_res <- purrr::map(domain_keep, fit_one_domain)
names(domain_res) <- domain_keep

slopes_dom  <- purrr::map_dfr(domain_res, "slopes")
effects_dom <- purrr::map_dfr(domain_res, "effects")
curve_dom   <- purrr::map_dfr(domain_res, "curve")
bins_dom    <- purrr::map_dfr(domain_res, "bins")
breaks_dom  <- purrr::map_dfr(domain_res, "breaks_tbl")

domain_levels <- domain_size |>
  filter(domain %in% domain_keep) |>
  pull(domain)

curve_dom  <- curve_dom  |> mutate(group_label = factor(group_label, domain_levels))
bins_dom   <- bins_dom   |> mutate(group_label = factor(group_label, domain_levels))
breaks_dom <- breaks_dom |> mutate(domain = factor(domain, domain_levels))
slopes_dom <- slopes_dom |> mutate(domain = factor(domain, domain_levels))

# ---- 8. 绘图主题 ----------------------------------------------------------

theme_paper <- function(base_size = 13) {
  theme_minimal(base_size = base_size) +
    theme(
      panel.grid.minor   = element_blank(),
      panel.grid.major   = element_line(colour = "grey92", linewidth = 0.3),
      panel.border       = element_rect(colour = "grey80", fill = NA, linewidth = 0.4),
      strip.background   = element_rect(fill = "grey95", colour = NA),
      strip.text         = element_text(face = "bold", size = base_size - 1),
      plot.title         = element_text(face = "bold", size = base_size + 3),
      plot.subtitle      = element_text(colour = "grey30", size = base_size),
      plot.caption       = element_text(colour = "grey35", hjust = 0, size = base_size - 3),
      legend.position    = "bottom",
      legend.title       = element_text(face = "bold"),
      axis.title         = element_text(size = base_size)
    )
}

# ---- 9. 图一：总体 --------------------------------------------------------

brk_all <- tibble::tibble(
  break_id  = c("Break 1", "Break 2"),
  break_day = psi_all
)

lab_all <- res_all$slopes |>
  mutate(
    x = (pmax(range_lo, min(analysis_dat$pr_duration)) + range_hi) / 2,
    txt = sprintf("%s: %+.4f%s", c("S1", "S2", "S3"), slope_per10, stars)
  ) |>
  rowwise() |>
  mutate(
    y = curve_all |>
      filter(pr_duration >= range_lo, pr_duration <= range_hi) |>
      pull(fit) |>
      mean()
  ) |>
  ungroup()

y_rng_all <- range(c(curve_all$fit, bins_all$y))
y_pad_top_all <- diff(y_rng_all) * 0.28
y_pad_bot_all <- diff(y_rng_all) * 0.16

p_overall <- ggplot() +
  geom_vline(
    data = brk_all,
    aes(xintercept = break_day, colour = break_id),
    linetype = "dashed", linewidth = 0.8, show.legend = TRUE
  ) +
  geom_point(
    data = bins_all, aes(x = x, y = y),
    colour = GREY_PT, size = 1.6, alpha = 0.75
  ) +
  geom_line(
    data = curve_all,
    aes(x = pr_duration, y = fit, colour = segment, group = segment),
    linewidth = 1.5
  ) +
  geom_label(
    data = brk_all,
    aes(x = break_day, y = y_rng_all[1] - y_pad_bot_all * 0.55,
        label = paste0("day ", round(break_day)), colour = break_id),
    fill = "white", label.size = 0.3, size = 3.9, fontface = "bold",
    show.legend = FALSE
  ) +
  geom_label(
    data = lab_all,
    aes(x = x, y = y, label = txt, colour = segment),
    fill = "white", label.size = 0.25, size = 4.1,
    vjust = -0.9, fontface = "bold", show.legend = FALSE
  ) +
  scale_colour_manual(
    values = c(SEG_COLS, BRK_COLS),
    breaks = c(names(SEG_COLS), names(BRK_COLS)),
    name   = NULL
  ) +
  scale_x_continuous(labels = scales::comma, expand = expansion(mult = c(0.02, 0.04))) +
  coord_cartesian(
    ylim = c(y_rng_all[1] - y_pad_bot_all, y_rng_all[2] + y_pad_top_all),
    clip = "off"
  ) +
  labs(
    x = "Peer review duration (days)",
    y = "Normalized citation impact"
  ) +
  guides(colour = guide_legend(nrow = 1, override.aes = list(linewidth = 1.6, label = ""))) +
  theme_paper(13) +
  theme(plot.margin = margin(t = 20, r = 20, b = 10, l = 10),
        axis.title.x = element_text(margin = margin(t = 15, b = 5)),  # 上边距15，下边距5
        axis.title.y = element_text(margin = margin(r = 15, l = 5)))  # 右边距15，左边距5

ggsave(
  file.path(out_dir, "fig1_overall_segmented.png"),
  p_overall,
  width = PARS$fig1_w, height = PARS$fig1_h, dpi = PARS$fig_dpi, bg = "white"
)

# ---- 10. 图二：分领域 2x5 分面 ----

panel_rng <- curve_dom |>
  bind_rows(bins_dom |> transmute(group_label, fit = y)) |>
  group_by(group_label) |>
  summarise(lo = min(fit, na.rm = TRUE), hi = max(fit, na.rm = TRUE), .groups = "drop") |>
  mutate(pad = (hi - lo) * 0.30)

x_max <- max(curve_dom$pr_duration)

slope_block <- slopes_dom |>
  mutate(tag = c("S1", "S2", "S3")[match(segment, names(SEG_COLS))]) |>
  group_by(domain) |>
  arrange(segment, .by_group = TRUE) |>
  mutate(line_no = dplyr::row_number()) |>
  ungroup() |>
  rename(group_label = domain) |>
  left_join(panel_rng, by = "group_label") |>
  mutate(
    txt = sprintf("%s %+.4f%s", tag, slope_per10, stars),
    x   = x_max * 0.97,
    y   = hi + pad * (1.02 - 0.30 * (line_no - 1))
  )

break_lab <- breaks_dom |>
  rename(group_label = domain) |>
  left_join(panel_rng, by = "group_label") |>
  mutate(
    txt = paste0("d", round(break_day)),
    y   = lo - pad * ifelse(break_id == "Break 1", 0.16, 0.42)
  )

p_domain <- ggplot() +
  geom_vline(
    data = breaks_dom |> rename(group_label = domain),
    aes(xintercept = break_day, colour = break_id),
    linetype = "dashed", linewidth = 0.7
  ) +
  geom_point(
    data = bins_dom, aes(x = x, y = y),
    colour = GREY_PT, size = 1.0, alpha = 0.7
  ) +
  geom_line(
    data = curve_dom,
    aes(x = pr_duration, y = fit, colour = segment, group = segment),
    linewidth = 1.25
  ) +
  geom_text(
    data = slope_block,
    aes(x = x, y = y, label = txt, colour = segment),
    hjust = 1, size = 3.15, fontface = "bold", show.legend = FALSE
  ) +
  geom_text(
    data = break_lab,
    aes(x = break_day, y = y, label = txt, colour = break_id),
    hjust = -0.12, size = 3.0, fontface = "bold", show.legend = FALSE
  ) +
  scale_colour_manual(
    values = c(SEG_COLS, BRK_COLS),
    breaks = c(names(SEG_COLS), names(BRK_COLS)),
    name   = NULL
  ) +
  scale_x_continuous(labels = scales::comma, expand = expansion(mult = c(0.03, 0.06))) +
  scale_y_continuous(expand = expansion(mult = c(0.16, 0.16))) +
  facet_wrap(
    ~ group_label,
    nrow   = PARS$facet_nrow,
    ncol   = PARS$facet_ncol,
    scales = "free_y"
  ) +
  labs(
    x = "Peer review duration (days)",
    y = "Normalized citation impact"
  ) +
  guides(colour = guide_legend(nrow = 1, override.aes = list(linewidth = 1.6))) +
  theme_paper(13) +
  theme(
    panel.spacing.x = unit(0.9, "lines"),
    panel.spacing.y = unit(1.3, "lines"),
    axis.title.x = element_text(margin = margin(t = 15, b = 5)),  # 上边距15，下边距5
    axis.title.y = element_text(margin = margin(r = 15, l = 5))   # 右边距15，左边距5
  )

ggsave(
  file.path(out_dir, "fig2_domain_comparison.png"),
  p_domain,
  width = PARS$fig2_w, height = PARS$fig2_h, dpi = PARS$fig_dpi, bg = "white"
)

# ggsave(
#   file.path(out_dir, "fig2_domain_comparison.pdf"),
#   p_domain,
#   width = PARS$fig2_w, height = PARS$fig2_h, device = cairo_pdf
# )
# 
# ggsave(
#   file.path(out_dir, "fig1_overall_segmented.pdf"),
#   p_overall,
#   width = PARS$fig1_w, height = PARS$fig1_h, device = cairo_pdf
# )

# ---- 11. 结果表格 ---------------------------------------------------------

slopes_out <- bind_rows(
  res_all$slopes |> mutate(domain = "Overall", psi_source = "own"),
  slopes_dom |> mutate(domain = as.character(domain))
) |>
  select(
    domain, segment, range_lo, range_hi,
    slope_per10, se_per10, t_stat, p_value, stars,
    n_obs, stable, psi_source
  )

breaks_out <- bind_rows(
  brk_all |> mutate(domain = "Overall"),
  breaks_dom |> mutate(domain = as.character(domain))
) |>
  select(domain, break_id, break_day)

effects_out <- bind_rows(
  eff_all |> mutate(domain = "Overall"),
  effects_dom |> mutate(domain = as.character(domain))
) |>
  select(domain, sd_cnp, span_full, span_full_sd, span_inner, span_inner_sd, r_squared)

write.csv(slopes_out,  file.path(out_dir, "segment_slopes.csv"),  row.names = FALSE)
write.csv(breaks_out,  file.path(out_dir, "break_points.csv"),    row.names = FALSE)
write.csv(effects_out, file.path(out_dir, "effect_sizes.csv"),    row.names = FALSE)

gt_slopes <- slopes_out |>
  gt::gt(groupname_col = "domain") |>
  gt::fmt_number(c(range_lo, range_hi), decimals = 0) |>
  gt::fmt_number(c(slope_per10, se_per10), decimals = 4) |>
  gt::fmt_number(t_stat, decimals = 1) |>
  gt::fmt_scientific(p_value, decimals = 2) |>
  gt::fmt_number(n_obs, decimals = 0, use_seps = TRUE) |>
  gt::cols_label(
    segment     = "Segment",
    range_lo    = "From (day)",
    range_hi    = "To (day)",
    slope_per10 = paste0("Slope / ", PARS$slope_per_days, "d"),
    se_per10    = "SE",
    t_stat      = "t",
    p_value     = "p",
    stars       = "",
    n_obs       = "n",
    stable      = "Stable",
    psi_source  = "Breaks"
  ) |>
  gt::tab_header(
    title    = "Segment slopes from two-break piecewise linear regression",
    subtitle = paste0(
      "Slopes are the change in normalized citation impact per ",
      PARS$slope_per_days, " additional days of review"
    )
  ) |>
  gt::tab_source_note(
    "Stable = no indicates fewer than 2,000 papers in that segment; treat the slope as unreliable."
  )

gt::gtsave(gt_slopes, file.path(out_dir, "segment_slopes.html"))
gt::gtsave(gt_slopes, file.path(out_dir, "segment_slopes.rtf"))

# ---- 12. 用 qs2 保存对象 -------------------------------------------------

qs_threads <- max(1L, parallel::detectCores() - 1L)

slim_lm <- function(m) {
  m$model         <- NULL
  m$residuals     <- NULL
  m$fitted.values <- NULL
  m$effects       <- NULL
  m$qr$qr         <- NULL
  m$weights       <- NULL
  m$prior.weights <- NULL
  m
}

model_bundle <- list(
  pars        = PARS,
  overall     = list(
    psi     = psi_all,
    fit     = slim_lm(res_all$fit),
    slopes  = res_all$slopes,
    effects = eff_all,
    n       = res_all$data_n
  ),
  by_domain   = purrr::map(domain_res, function(x) {
    list(
      domain     = x$domain,
      psi        = x$psi,
      psi_source = x$psi_source,
      fit        = slim_lm(x$fit),
      slopes     = x$slopes,
      effects    = x$effects
    )
  }),
  domain_size = domain_size
)

qs2::qs_save(
  model_bundle,
  file.path(out_dir, "segmented_models.qs2"),
  nthreads = qs_threads
)

plot_bundle <- list(
  curve_all  = curve_all,
  bins_all   = bins_all,
  brk_all    = brk_all,
  curve_dom  = curve_dom,
  bins_dom   = bins_dom,
  breaks_dom = breaks_dom,
  slopes_out = slopes_out,
  seg_cols   = SEG_COLS,
  brk_cols   = BRK_COLS
)

qs2::qs_save(
  plot_bundle,
  file.path(out_dir, "plot_data.qs2"),
  nthreads = qs_threads
)

message("\n完成。所有结果已写入：", normalizePath(out_dir))
message("图形：fig1_overall_segmented.(png|pdf)、fig2_domain_comparison.(png|pdf)")
message("表格：segment_slopes.(csv|html|rtf)、break_points.csv、effect_sizes.csv")
message("对象：segmented_models.qs2、plot_data.qs2")
message("\n缓存位置：", normalizePath(cache_dir))
message("后续运行将自动读取缓存中的断点，确保完全复现。")
