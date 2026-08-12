
rm(list = ls())

library(tidyfst)
pkg_load(tidyverse,fs,arrow,DBI,duckdb,scales)

con<- dbConnect(
  duckdb(),
  dbdir = "data/sel_data.duckdb"
)

tbl(con,"review_citation") %>% 
  select(UT,retracted,pr_duration) %>% 
  summarise(
    retracted = sum(retracted,na.rm = T),
    n = n(),
    .by = pr_duration
  ) %>% 
  arrange(pr_duration) %>% 
  collect() -> retracted_data

# ============================================================
# 审稿周期与撤稿率关系：分箱可视化 + "拐点"（曲线变平点）检测
# Relationship between peer-review duration and retraction rate:
# binned visualization + detection of the "flattening point" of the smooth curve
#
# 数据：retracted_data，列为 pr_duration（天，逐天）、retracted（撤稿数）、n（论文数）
# Data: retracted_data, columns are pr_duration (days, day-by-day),
#       retracted (number of retracted papers), n (total number of papers)
# ============================================================

library(tidyverse)
library(mgcv)      # 拟合平滑曲线（logistic GAM） / fit smooth curve via logistic GAM
library(scales)

# ---------- 0. 读入并检查数据 / Load and check data ----------
# retracted_data <- readRDS("retracted_data.rds")
stopifnot(all(c("pr_duration","retracted","n") %in% names(retracted_data)))

dat_day <- retracted_data %>%
  filter(pr_duration >= 0, n > 0) %>%
  arrange(pr_duration)

# ============================================================
# 1. 用 logistic GAM 在"逐天"数据上拟合平滑曲线
#    —— 这一步只是为了得到一条稳定、可求导的曲线，
#       和你之前散点图里的 geom_smooth 是同一条逻辑，
#       但用二项分布 + 样本量加权，避免长尾小样本天数的噪声主导拟合。
#
# 1. Fit a smooth curve with a logistic GAM on the day-by-day data
#    -- This step is only to obtain a stable, differentiable curve,
#       following the same logic as geom_smooth() in your earlier scatter plot,
#       but using a binomial family + sample-size weighting so that noisy,
#       small-sample days in the long tail don't dominate the fit.
# ============================================================

gam_fit <- gam(
  cbind(retracted, n - retracted) ~ s(pr_duration, k = 30, bs = "tp"),
  family = binomial(link = "logit"),
  data   = dat_day,
  method = "REML"
)

# 限定检测范围：避免极端尾部（样本量稀疏）产生的曲线抖动被误判为拐点
# Restrict the detection range: avoid curve jitter from the sparse extreme tail
# being mistaken for the flattening point
cum_n   <- cumsum(dat_day$n) / sum(dat_day$n)
lo_cut  <- dat_day$pr_duration[which(cum_n >= 0.01)[1]]
hi_cut  <- dat_day$pr_duration[which(cum_n >= 0.90)[1]]

grid <- tibble(pr_duration = seq(lo_cut, hi_cut, by = 1))

pred <- predict(gam_fit, newdata = grid, type = "link", se.fit = TRUE)
grid <- grid %>%
  mutate(
    logit_fit = pred$fit,
    prob_fit  = plogis(logit_fit)
  )

# ============================================================
# 2. 检测"红线斜率由负转为 0"的点
#    定义：曲线（概率尺度）的一阶导数由负变为 ~0（且此后不再显著为负）
#    做法：
#      a) 数值微分得到 d(prob)/d(day)
#      b) 找到导数首次从负值区间进入"接近 0 的稳定平台"的位置
#         —— 用一个更稳健的判据：导数序列首次达到
#            其自身负向峰值的 5%（即已衰减 95%）且此后保持在该阈值以内
#
# 2. Detect the point where the red curve's slope turns from negative to ~0
#    Definition: the first derivative of the curve (on the probability scale)
#    goes from negative to approximately zero (and stays non-negative afterward)
#    Approach:
#      a) Numerically differentiate to get d(prob)/d(day)
#      b) Find where the derivative first moves from the negative region
#         into a "stable near-zero plateau"
#         -- using a more robust criterion: the first point where the
#            derivative reaches 5% of its own negative peak (i.e., has decayed
#            by 95%) and stays within that threshold afterward
# ============================================================

d1 <- diff(grid$prob_fit) / diff(grid$pr_duration)
d1_at <- grid$pr_duration[-1] - 0.5   # 导数对应的中点位置 / midpoint position corresponding to each derivative value

# 只关注前半段的"下降期"（导数为负的主体部分）
# Focus only on the early "declining phase" (the main negative-derivative segment)
neg_peak <- min(d1, na.rm = TRUE)     # 最陡的下降斜率（负值最大） / steepest downward slope (most negative value)
thresh   <- neg_peak * 0.05           # 衰减到峰值的 5% 视为"基本变平" / decayed to 5% of the peak is considered "essentially flat"

# 找到导数首次从 neg_peak 附近回升、并越过 thresh（即从很负变得接近 0）
# 且此后一段窗口内（比如 20 天）导数都不再跌破 thresh，避免抖动误判
# Find the first point where the derivative rises from near neg_peak and
# crosses thresh (i.e., moves from very negative to near zero), and stays
# above thresh for a following window (e.g., 20 days), to avoid false
# positives caused by jitter
window_days <- 20
flatten_day <- NA_real_
for (i in seq_along(d1)) {
  if (d1[i] >= thresh) {
    idx_end <- which(d1_at >= d1_at[i] + window_days)[1]
    if (is.na(idx_end)) idx_end <- length(d1)
    if (all(d1[i:idx_end] >= thresh * 1.5, na.rm = TRUE) || i == length(d1)) {
      flatten_day <- d1_at[i]
      break
    }
  }
}

# 兜底：如果上面的稳健判据没找到，直接退化为曲线最小值点（导数恰好为 0 之处）
# Fallback: if the robust criterion above finds nothing, fall back to the
# curve's minimum point (where the derivative is exactly zero)
if (is.na(flatten_day)) {
  flatten_day <- grid$pr_duration[which.min(grid$prob_fit)]
}

cat(sprintf("检测到的拐点（曲线由降转平）大约在第 %.0f 天 / Detected flattening point (curve turns from declining to flat): approximately day %.0f\n", flatten_day, flatten_day))

# # 可选：画一下导数曲线，直观检查这个判据是否合理
# # Optional: plot the derivative curve to visually check whether this criterion is reasonable
# p_deriv <- ggplot(tibble(day = d1_at, deriv = d1), aes(day, deriv)) +
#   geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
#   geom_hline(yintercept = thresh, linetype = "dotted", color = "firebrick") +
#   geom_line(color = "steelblue", linewidth = 0.8) +
#   geom_vline(xintercept = flatten_day, color = "purple", linetype = "dashed") +
#   labs(x = "Peer review duration (days)",
#        y = "d(predicted retraction probability) / d(day)",
#        caption = sprintf("Dotted line = 5%% of peak negative slope threshold; purple = detected flattening point (day %.0f)", flatten_day)) +
#   theme_bw(base_size = 13)
# p_deriv
# ggsave("derivative_check.png", p_deriv, width = 7, height = 5, dpi = 300)

# ============================================================
# 3. 分箱（等频，按累计样本量），用于可视化点和误差棒
# 3. Binning (equal-frequency, based on cumulative sample size),
#    used for the visualization points and error bars
# ============================================================

n_bins <- 30

dat_binned <- dat_day %>%
  mutate(
    cum_n  = cumsum(n),
    bin_id = ceiling(cum_n / (sum(n) / n_bins))
  ) %>%
  group_by(bin_id) %>%
  summarise(
    day_mid    = weighted.mean(pr_duration, w = n),
    day_min    = min(pr_duration),
    day_max    = max(pr_duration),
    n_sum      = sum(n),
    retr_sum   = sum(retracted),
    rate       = retr_sum / n_sum,
    .groups = "drop"
  ) %>%
  rowwise() %>%
  mutate(
    # Wilson 95% 置信区间，比正态近似在小样本/极端比例下更稳健
    # Wilson 95% confidence interval, more robust than the normal approximation
    # for small samples or extreme proportions
    ci = list(binom::binom.confint(retr_sum, n_sum, methods = "wilson")),
    ci_low  = ci$lower,
    ci_high = ci$upper
  ) %>%
  ungroup() %>%
  select(-ci)

# ============================================================
# 4. 主图：分箱散点 + GAM 平滑曲线 + 拐点标注
# 4. Main plot: binned scatter points + GAM smooth curve + flattening-point annotation
# ============================================================

gam_curve <- grid %>%
  mutate(
    se       = pred$se.fit,
    prob_low  = plogis(logit_fit - 1.96 * se),
    prob_high = plogis(logit_fit + 1.96 * se)
  )

p_main <- ggplot() +
  geom_ribbon(data = gam_curve, aes(x = pr_duration, ymin = prob_low, ymax = prob_high),
              fill = "firebrick", alpha = 0.15) +
  geom_line(data = gam_curve, aes(x = pr_duration, y = prob_fit),
            color = "firebrick", linewidth = 1) +
  geom_errorbar(data = dat_binned, aes(x = day_mid, ymin = ci_low, ymax = ci_high),
                width = 0, color = "grey40", alpha = 0.6) +
  geom_point(data = dat_binned, aes(x = day_mid, y = rate, size = n_sum),
             color = "navy", alpha = 0.55) +
  geom_vline(xintercept = flatten_day, color = "purple", linetype = "dashed", linewidth = 0.8) +
  annotate("text", x = flatten_day, y = max(dat_binned$ci_high, na.rm = TRUE) * 0.97,
           label = sprintf("day %.0f", flatten_day), color = "purple",
           hjust = -0.15, fontface = "bold", size = 4.2) +
  scale_y_continuous(labels = percent_format(accuracy = 0.001)) +
  scale_size_continuous(name = "Number of papers\nin bin", labels = comma) +
  labs(
    x = "Peer review duration (days)",
    y = "Retraction rate"
  ) +
  theme_bw(base_size = 13) +
  theme(
    plot.caption = element_text(hjust = 0, size = 9, color = "grey30"),
    legend.position = "right",
    axis.title.x = element_text(margin = margin(t = 15, b = 5)),  # 上边距15，下边距5
    axis.title.y = element_text(margin = margin(r = 15, l = 5))   # 右边距15，左边距5
  )

print(p_main)
# ggsave("retraction_rate_flattening_point.png", p_main, width = 9, height = 6, dpi = 400)
# ggsave("retraction_rate_flattening_point.pdf", p_main, width = 9, height = 6)

# ============================================================
# 5. 截断版（排除极端长尾，便于展示主要区间细节）
# 5. Truncated version (excludes the extreme long tail, for clearer detail
#    on the main range)
# ============================================================

# 截断依据：以“论文篇数”为权重（而非天数本身）计算审稿周期的分位数，
# 取第 90% 分位数作为横轴上限——即保留审稿周期最短、样本最密集的 90% 的论文，
# 把剩下 10% 分布在稀疏长尾里的论文排除在图外，从而避免长尾拖长横轴、掩盖主体区间的细节。
# 注意：这里的 cut_q（90% 分位数）与第 1 节里 GAM 拟合/拐点检测用的 hi_cut（同样是 90% 分位数）
# 数值上一致，但含义和用途是两回事：hi_cut 是为了让拐点检测不被稀疏尾部的抖动干扰，
# 这里的 cut_q 只是为了让展示图更聚焦、裁掉长尾，不影响任何模型拟合或拐点计算。
# 如果以后两者需要分别调整（比如拐点检测范围和图形展示范围取不同的分位数），
# 直接修改本节的 0.90 即可，不会影响第 1 节的检测结果。
#
# Truncation basis: compute the quantile of peer-review duration weighted by
# the number of papers (not by day count itself), and take the 90th percentile
# as the x-axis upper limit -- i.e., keep the 90% of papers with the shortest
# review durations / densest sample coverage, and exclude the remaining 10%
# that are scattered across the sparse long tail. This prevents the long tail
# from stretching the x-axis and obscuring detail in the main range.
# Note: this cut_q (90th percentile) happens to share the same numeric value
# as hi_cut (also 90th percentile) used in Section 1 for the GAM fit /
# flattening-point detection, but they serve different purposes: hi_cut keeps
# the breakpoint detection from being disturbed by tail jitter, while cut_q
# here is purely a display choice for a more focused plot and does not affect
# any model fitting or flattening-point calculation. If the two ever need to
# diverge (e.g., a different quantile for detection vs. display), just edit
# the 0.90 here independently -- it won't change Section 1's results.
cut_q <- quantile(rep(dat_day$pr_duration, times = dat_day$n), 0.90)
# 如果 90% 分位数仍然偏长（即长尾依然可见），可以尝试调低到 0.85、0.80 等更小的分位数
# If the 90% quantile is still too long (i.e., the tail is still visible),
# try lowering it further to 0.85, 0.80, or an even smaller quantile

p_trunc <- p_main +
  coord_cartesian(xlim = c(0, cut_q))

print(p_trunc)
ggsave("fig/03_retraction_rate_flattening_point_truncated.png", p_trunc, width = 9, height = 6, dpi = 400)

caption_text <- sprintf(
  "Dashed line: detected flattening point (day %.0f), defined as where the smooth curve's slope\ndecays to ~5%% of its steepest negative value and remains flat thereafter.\nPoints: quantile-based bins (equal cumulative sample size); point size = number of papers per bin;\nerror bars = Wilson 95%% CI. Red line/band: logistic GAM fit with 95%% CI.",
  flatten_day
)

cat(caption_text, "\n")


