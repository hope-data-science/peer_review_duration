rm(list = ls())

library(tidyfst)
pkg_load(tidyverse,fs,arrow,DBI,duckdb,scales)

con<- dbConnect(
  duckdb(),
  dbdir = "data/sel_data.duckdb"
)

# tbl(con,"review_citation") %>% 
#   filter(document_type == "Article") %>% 
#   select(DOI,received_dt,pr_duration,journal,COVID_mark) -> dat_raw
# 
# 
# tbl(con,"review_citation") %>% 
#   filter(document_type == "Article")
#   filter(domain %in% c("Medicine")) %>% 
#   collect() -> all_data
# 
# all_data %>% 
#   select(DOI,received_dt,pr_duration) -> dat

# We defined 11 March 2020—the date on which the World Health Organization 
# characterized COVID-19 as a pandemic—as the primary interruption point.
# https://www.who.int/news/item/27-04-2020-who-timeline---covid-19

# We use 30 November 2022 as the main intervention date 
# because it marks the public launch of ChatGPT, 
# the first widely accessible general-purpose LLM; 
# its effects on peer-review duration are expected to 
# emerge gradually as researchers, 
# editors, and reviewers adopted the tool.

# =============================================================================
# COVID 对同行评审时长的影响：核心回归 + 领域异质性
# -----------------------------------------------------------------------------
# 研究对象：最终被接收的 article。结果仅代表"最终接收论文的投稿至接收时长"。
# 时间轴：received_dt（投稿日）。
#
# 本脚本分为两个独立部分，可分别运行：
#   PART 1 —— 核心回归。回答两个问题：
#              (1) COVID 是否影响审稿周期，影响多大？
#              (2) COVID 主题论文受到的影响是否更大？
#              期刊固定效应作为设定的一部分直接纳入，不再单独论证。
#              输出：Core_results.xlsx（4 张表）
#   PART 2 —— 领域异质性。分领域重估，输出可视化。
#              输出：Figure_domain.png / .pdf + Domain_results.xlsx
#
# 前置假定：期刊固定效应的必要性、样本筛选规则的稳健性，已由前期预实验确认，
#           此处直接采用，不再重复检验。
# =============================================================================

# ============================ 0. 环境与参数 ==================================
required_pkgs <- c("dplyr", "readr", "lubridate", "stringr", "fixest",
                   "tidyr", "openxlsx", "qs2", "ggplot2", "forcats")
missing_pkgs <- required_pkgs[!vapply(required_pkgs, requireNamespace,
                                      logical(1), quietly = TRUE)]
if (length(missing_pkgs) > 0) {
  install.packages(missing_pkgs, repos = "https://cloud.r-project.org")
}
invisible(lapply(required_pkgs, library, character.only = TRUE))

OUTPUT_DIR <- "data/covid_review_results"
# 每次运行前清空同名目录，避免旧结果残留与新结果混淆
if (dir.exists(OUTPUT_DIR)) {
  unlink(OUTPUT_DIR, recursive = TRUE, force = TRUE)
}
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

# --- 关键参数（前期预实验已确定，此处固定） ---
COVID_START           <- as.Date("2020-03-11")  # WHO 宣布大流行
DATA_EXTRACTION_DATE  <- as.Date("2026-01-01")  # 数据抓取日；如有确切日期请替换
FOLLOWUP_DAYS         <- 730L                   # 充分随访窗口：2 年
MAX_DURATION_DAYS     <- 1460L                  # 剔除超过 4 年的极端时长
MIN_JOURNAL_PAPERS_YR <- 20L                    # 稳定期刊：每活跃年 >= 20 篇
MIN_JOURNAL_YEARS_PRE <- 2L                     # 疫情前至少 2 个活跃年
MIN_JOURNAL_YEARS_POST<- 2L                     # 疫情后至少 2 个活跃年
MIN_DOMAIN_PAPERS     <- 5000L                  # 领域分析的最小样本量门槛

followup_cutoff <- DATA_EXTRACTION_DATE - FOLLOWUP_DAYS

# ============================ 1. 读入数据 =====================================
tbl(con, "review_citation") %>%
  filter(document_type == "Article") %>%
  select(DOI, received_dt, pr_duration, journal, COVID_mark, domain) %>%
  collect() -> dat_raw

required_cols <- c("DOI", "received_dt", "pr_duration", "journal",
                   "COVID_mark", "domain")
missing_cols <- setdiff(required_cols, names(dat_raw))
if (length(missing_cols) > 0) {
  stop("数据缺少字段：", paste(missing_cols, collapse = ", "))
}

dat0 <- dat_raw %>%
  mutate(
    DOI         = as.character(DOI),
    journal     = str_squish(as.character(journal)),
    domain      = str_squish(as.character(domain)),
    received_dt = as.Date(received_dt),
    pr_duration = as.numeric(pr_duration),
    COVID_mark  = as.integer(COVID_mark)
  )

# ============================ 2. 预处理 =======================================
# 逐步记录样本量，供方法学部分引用
flow <- tibble(step = "原始 article 记录", n = nrow(dat0))

dat1 <- dat0 %>%
  filter(!is.na(DOI), DOI != "", !is.na(received_dt), !is.na(pr_duration),
         !is.na(journal), journal != "", COVID_mark %in% c(0L, 1L))
flow <- add_row(flow, step = "核心变量完整", n = nrow(dat1))

dat2 <- dat1 %>%
  arrange(DOI, pr_duration) %>%
  distinct(DOI, .keep_all = TRUE)
flow <- add_row(flow, step = "DOI 去重", n = nrow(dat2))

dat3 <- dat2 %>%
  filter(pr_duration > 0, pr_duration <= MAX_DURATION_DAYS)
flow <- add_row(flow, step = paste0("时长在 1-", MAX_DURATION_DAYS, " 天"),
                n = nrow(dat3))

# 充分随访：只保留投稿早于"抓取日 - 730 天"的记录，
# 避免近期投稿中"只有审得快的才已被接收"造成的选择性偏差
dat4 <- dat3 %>%
  filter(received_dt <= followup_cutoff)
flow <- add_row(flow, step = paste0("充分随访（投稿 <= ", followup_cutoff, "）"),
                n = nrow(dat4))

# 构造时间与暴露变量
dat4 <- dat4 %>%
  mutate(
    submit_month     = floor_date(received_dt, "month"),
    submit_year      = year(received_dt),
    submit_month_num = month(received_dt),          # 1-12 月季节性
    post_covid       = as.integer(received_dt >= COVID_START),
    months_from_covid= interval(COVID_START, received_dt) %/% months(1),
    log_duration     = log1p(pr_duration),
    # COVID 分段项：允许"先剧变、后适应、再趋稳/反转"的非单调过程
    covid_level      = post_covid,
    covid_0_3m       = pmax(0, pmin(months_from_covid, 3)),
    covid_3_12m      = pmax(0, pmin(months_from_covid - 3, 9)),
    covid_12_24m     = pmax(0, pmin(months_from_covid - 12, 12)),
    covid_after_24m  = pmax(0, months_from_covid - 24),
    # 主题论文的疫情后额外差异
    covid_mark_post  = COVID_mark * post_covid
  )

# 稳定期刊样本：排除疫情前后期刊进入/退出造成的构成变化
journal_year <- dat4 %>%
  count(journal, submit_year, name = "n_papers") %>%
  mutate(
    period = ifelse(as.Date(paste0(submit_year, "-01-01")) < COVID_START,
                    "pre", "post"),
    active = n_papers >= MIN_JOURNAL_PAPERS_YR
  )

stable_journals <- journal_year %>%
  filter(active) %>%
  count(journal, period, name = "n_active_years") %>%
  pivot_wider(names_from = period, values_from = n_active_years,
              values_fill = 0) %>%
  filter(pre >= MIN_JOURNAL_YEARS_PRE, post >= MIN_JOURNAL_YEARS_POST) %>%
  pull(journal)

dat <- dat4 %>% filter(journal %in% stable_journals)
flow <- add_row(flow, step = "稳定期刊主样本", n = nrow(dat))

# 疫情前线性趋势（疫情后冻结在断点水平，作为反事实外推的基准）
dat <- dat %>%
  mutate(
    month_index = interval(min(submit_month), submit_month) %/% months(1),
    pre_trend_month = pmin(
      month_index,
      max(month_index[received_dt < COVID_START], na.rm = TRUE)
    )
  )

cat("=== 样本构建完成 ===\n")
cat("充分随访投稿截止日：", as.character(followup_cutoff), "\n")
cat("稳定期刊数量：", length(stable_journals), "\n")
cat("最终分析样本：", nrow(dat), "\n\n")
print(flow)

# ========================= 3. 通用工具函数 ====================================
COVID_TERMS <- "covid_level + covid_0_3m + covid_3_12m + covid_12_24m + covid_after_24m"
PIECE_NAMES <- c("covid_level", "covid_0_3m", "covid_3_12m",
                 "covid_12_24m", "covid_after_24m")

# 构造第 k 个月的分段项取值向量（用于累计效应的线性组合）
piece_values <- function(k) {
  c(covid_level     = 1,
    covid_0_3m      = min(k, 3),
    covid_3_12m     = max(min(k - 3, 9), 0),
    covid_12_24m    = max(min(k - 12, 12), 0),
    covid_after_24m = max(k - 24, 0))
}

# 对任意线性组合做点估计、标准误、置信区间与 p 值
lincom <- function(model, weights) {
  b <- coef(model)
  V <- vcov(model)
  L <- setNames(rep(0, length(b)), names(b))
  keep <- intersect(names(weights), names(b))
  L[keep] <- weights[keep]
  est <- sum(L * b)
  se  <- sqrt(as.numeric(t(L) %*% V %*% L))
  z   <- est / se
  tibble(
    log_est = est,
    se      = se,
    pct     = 100 * (exp(est) - 1),
    lo_pct  = 100 * (exp(est - 1.96 * se) - 1),
    hi_pct  = 100 * (exp(est + 1.96 * se) - 1),
    p_value = 2 * pnorm(-abs(z))
  )
}

# 累计效应：第 k 个月相对疫情前趋势外推的偏离
cumulative_effect <- function(model, k, extra = NULL) {
  w <- piece_values(k)
  if (!is.null(extra)) w <- c(w, extra)
  lincom(model, w)
}

# 数值格式化：系数(标准误) + 显著性星号
star <- function(p) {
  ifelse(p < 0.01, "***", ifelse(p < 0.05, "**", ifelse(p < 0.1, "*", "")))
}
fmt_coef <- function(est, se, p, digits = 4) {
  sprintf("%.*f%s (%.*f)", digits, est, star(p), digits, se)
}
fmt_pct <- function(pct, lo, hi, p) {
  sprintf("%.2f%s [%.2f, %.2f]", pct, star(p), lo, hi)
}

# =============================================================================
#                        PART 1  核心回归（不涉及领域）
# =============================================================================

# ---- 3.1 模型 ----
# M1：COVID 主效应。期刊 FE + 月份季节性 FE + 疫情前趋势。
m1 <- feols(
  as.formula(paste0("log_duration ~ pre_trend_month + ", COVID_TERMS,
                    " | journal + submit_month_num")),
  data = dat, cluster = ~journal
)

# M2：主模型。加入 COVID_mark 及其与疫情后的交互，识别主题论文的额外差异。
m2 <- feols(
  as.formula(paste0("log_duration ~ pre_trend_month + ", COVID_TERMS,
                    " + COVID_mark + covid_mark_post | journal + submit_month_num")),
  data = dat, cluster = ~journal
)

# M3：M2 的双向聚类版本，允许同一投稿月的共同冲击。作为标准误稳健性附列。
m3 <- feols(
  as.formula(paste0("log_duration ~ pre_trend_month + ", COVID_TERMS,
                    " + COVID_mark + covid_mark_post | journal + submit_month_num")),
  data = dat, cluster = ~journal + submit_month
)

qs2::qs_save(list(m1 = m1, m2 = m2, m3 = m3),
             file.path(OUTPUT_DIR, "core_models.qs"))

# ---- 表 1：描述统计（疫情前后 × 是否 COVID 主题） ----
T1 <- dat %>%
  mutate(
    Period = ifelse(post_covid == 1, "Post-COVID", "Pre-COVID"),
    Topic  = ifelse(COVID_mark == 1, "COVID-related", "Non-COVID-related")
  ) %>%
  group_by(Period, Topic) %>%
  summarise(
    N            = n(),
    Mean_days    = mean(pr_duration),
    SD_days      = sd(pr_duration),
    Q1_days      = quantile(pr_duration, 0.25),
    Median_days  = median(pr_duration),
    Q3_days      = quantile(pr_duration, 0.75),
    .groups = "drop"
  ) %>%
  arrange(desc(Period), Topic) %>%
  mutate(across(c(Mean_days, SD_days, Q1_days, Median_days, Q3_days),
                ~ round(.x, 1)))

# ---- 表 2：主回归系数 ----
coef_row <- function(model, term, label) {
  b <- coef(model); V <- vcov(model)
  if (!term %in% names(b)) return(tibble(Variable = label, v = "—"))
  est <- b[[term]]; se <- sqrt(V[term, term])
  p <- 2 * pnorm(-abs(est / se))
  tibble(Variable = label, v = fmt_coef(est, se, p))
}

term_labels <- tibble(
  term = c("pre_trend_month", "covid_level", "covid_0_3m", "covid_3_12m",
           "covid_12_24m", "covid_after_24m", "COVID_mark", "covid_mark_post"),
  label = c("Pre-pandemic monthly trend",
            "Immediate level change at COVID onset",
            "Monthly slope, months 0-3",
            "Monthly slope, months 3-12",
            "Monthly slope, months 12-24",
            "Monthly slope, after month 24",
            "COVID-related topic (pre-pandemic baseline)",
            "COVID-related topic x Post-COVID")
)

build_col <- function(model, colname) {
  body <- lapply(seq_len(nrow(term_labels)), function(i) {
    coef_row(model, term_labels$term[i], term_labels$label[i])
  }) %>% bind_rows()
  foot <- tibble(
    Variable = c("Journal fixed effects", "Calendar-month fixed effects",
                 "Observations", "R2", "Within R2"),
    v = c("Yes", "Yes",
          format(nobs(model), big.mark = ","),
          sprintf("%.4f", fitstat(model, "r2")$r2),
          sprintf("%.4f", fitstat(model, "wr2")$wr2))
  )
  out <- bind_rows(body, foot)
  names(out)[2] <- colname
  out
}

c1 <- build_col(m1, "(1) COVID only")
c2 <- build_col(m2, "(2) Main model")
c3 <- build_col(m3, "(3) Two-way cluster")

T2 <- c1 %>%
  left_join(c2, by = "Variable") %>%
  left_join(c3, by = "Variable")

# ---- 表 3：COVID 累计动态效应（结论一：影响有多大） ----
key_months <- c(0, 1, 3, 6, 12, 18, 24, 30, 36, 48)
T3 <- lapply(key_months, function(k) {
  cumulative_effect(m2, k) %>%
    mutate(Months_since_onset = k, .before = 1)
}) %>%
  bind_rows() %>%
  transmute(
    Months_since_onset,
    Change_pct   = round(pct, 2),
    CI_lower_pct = round(lo_pct, 2),
    CI_upper_pct = round(hi_pct, 2),
    P_value      = signif(p_value, 3),
    Significance = star(p_value)
  )

# ---- 表 4：COVID 主题论文的额外影响（结论二） ----
# Panel A：主题差异的三个关键对比
pa_pre  <- lincom(m2, c(COVID_mark = 1))
pa_post <- lincom(m2, c(COVID_mark = 1, covid_mark_post = 1))
pa_diff <- lincom(m2, c(covid_mark_post = 1))

T4a <- bind_rows(
  pa_pre  %>% mutate(Contrast = "COVID-topic vs non-COVID-topic, pre-pandemic"),
  pa_post %>% mutate(Contrast = "COVID-topic vs non-COVID-topic, post-pandemic"),
  pa_diff %>% mutate(Contrast = "Change in the topic gap (interaction)")
) %>%
  transmute(
    Panel = "A. Topic gap",
    Contrast,
    Estimate_log = round(log_est, 4),
    SE           = round(se, 4),
    Change_pct   = round(pct, 2),
    CI_lower_pct = round(lo_pct, 2),
    CI_upper_pct = round(hi_pct, 2),
    P_value      = signif(p_value, 3),
    Significance = star(p_value)
  )

# Panel B：两类论文各自的累计效应（便于直接比较幅度）
panel_b_months <- c(3, 12, 24)
T4b <- lapply(panel_b_months, function(k) {
  non_covid <- cumulative_effect(m2, k) %>%
    mutate(Contrast = paste0("Non-COVID-topic papers, month ", k))
  covid_tp <- cumulative_effect(m2, k,
                                extra = c(COVID_mark = 1, covid_mark_post = 1)) %>%
    mutate(Contrast = paste0("COVID-topic papers, month ", k))
  bind_rows(non_covid, covid_tp)
}) %>%
  bind_rows() %>%
  transmute(
    Panel = "B. Cumulative effect by topic",
    Contrast,
    Estimate_log = round(log_est, 4),
    SE           = round(se, 4),
    Change_pct   = round(pct, 2),
    CI_lower_pct = round(lo_pct, 2),
    CI_upper_pct = round(hi_pct, 2),
    P_value      = signif(p_value, 3),
    Significance = star(p_value)
  )

T4 <- bind_rows(T4a, T4b)

# ---- 写入 Core_results.xlsx ----
core_readme <- tibble(
  Sheet = c("T1_Descriptives", "T2_Main_regression",
            "T3_Dynamic_effects", "T4_Topic_effects", "S0_Sample_flow"),
  Suggested_place = c("Table 1", "Table 2", "Table 3", "Table 4",
                      "Methods / Appendix"),
  Supports_conclusion = c(
    "样本概况",
    "结论一与结论二的系数来源",
    "结论一：COVID 影响的幅度与时间演化",
    "结论二：COVID 主题论文受到的影响更大",
    "方法学：样本构建过程"
  ),
  Note_for_manuscript = c(
    "Values are submission-to-acceptance days among eventually accepted articles.",
    "Dependent variable is log(1 + days). SE clustered by journal. Column (3) uses two-way clustering. *** p<0.01, ** p<0.05, * p<0.1.",
    "Cumulative model-implied deviation from the extrapolated pre-pandemic trend, in percent. 95% CI from journal-clustered SE.",
    "Panel A reports the topic gap; the post-pandemic gap is a joint linear combination, not a single coefficient. Panel B reports each group's cumulative deviation.",
    "Stepwise sample construction under the stable-journal and adequate-follow-up rules."
  )
)

wb1 <- createWorkbook()
hs <- createStyle(textDecoration = "bold", halign = "center",
                  fgFill = "#F2F2F2", border = "TopBottom")
add_sheet <- function(wb, name, df) {
  addWorksheet(wb, name)
  writeData(wb, name, df, headerStyle = hs)
  freezePane(wb, name, firstRow = TRUE)
  setColWidths(wb, name, cols = seq_len(ncol(df)), widths = "auto")
}
add_sheet(wb1, "README", core_readme)
add_sheet(wb1, "T1_Descriptives", T1)
add_sheet(wb1, "T2_Main_regression", T2)
add_sheet(wb1, "T3_Dynamic_effects", T3)
add_sheet(wb1, "T4_Topic_effects", T4)
add_sheet(wb1, "S0_Sample_flow", flow)
saveWorkbook(wb1, file.path(OUTPUT_DIR, "Core_results.xlsx"), overwrite = TRUE)

cat("\n=== PART 1 完成：Core_results.xlsx ===\n")
cat("COVID 即时效应：", sprintf("%.2f%%", T3$Change_pct[T3$Months_since_onset == 0]), "\n")
cat("第 12 个月累计：", sprintf("%.2f%%", T3$Change_pct[T3$Months_since_onset == 12]), "\n")
cat("疫情后主题差距：", sprintf("%.2f%% (p = %.3g)", pa_post$pct, pa_post$p_value), "\n\n")

# =============================================================================
#                     PART 2  领域异质性（可视化）
# =============================================================================

# ---- 4.1 领域筛选与分领域建模 ----
domain_ok <- dat %>%
  filter(!is.na(domain), domain != "") %>%
  count(domain, name = "n") %>%
  filter(n >= MIN_DOMAIN_PAPERS) %>%
  pull(domain)

cat("纳入领域分析的 domain 数量：", length(domain_ok), "\n")

fit_domain <- function(d) {
  sub <- dat %>% filter(domain == d)
  # 领域内仍保留期刊固定效应，比较的是同一期刊内不同投稿时点的变化
  has_topic <- sum(sub$COVID_mark) >= 100 &&
    sum(sub$COVID_mark * sub$post_covid) >= 50
  f <- if (has_topic) {
    paste0("log_duration ~ pre_trend_month + ", COVID_TERMS,
           " + COVID_mark + covid_mark_post | journal + submit_month_num")
  } else {
    paste0("log_duration ~ pre_trend_month + ", COVID_TERMS,
           " | journal + submit_month_num")
  }
  m <- try(feols(as.formula(f), data = sub, cluster = ~journal), silent = TRUE)
  if (inherits(m, "try-error")) return(NULL)
  
  e12 <- cumulative_effect(m, 12)
  e24 <- cumulative_effect(m, 24)
  gap <- if (has_topic) {
    lincom(m, c(COVID_mark = 1, covid_mark_post = 1))
  } else {
    tibble(log_est = NA, se = NA, pct = NA, lo_pct = NA,
           hi_pct = NA, p_value = NA)
  }
  
  tibble(
    domain      = d,
    n_papers    = nrow(sub),
    n_journals  = n_distinct(sub$journal),
    eff12_pct   = e12$pct, eff12_lo = e12$lo_pct,
    eff12_hi    = e12$hi_pct, eff12_p = e12$p_value,
    eff24_pct   = e24$pct, eff24_lo = e24$lo_pct,
    eff24_hi    = e24$hi_pct, eff24_p = e24$p_value,
    gap_pct     = gap$pct, gap_lo = gap$lo_pct,
    gap_hi      = gap$hi_pct, gap_p = gap$p_value,
    topic_model = has_topic
  )
}

domain_res <- lapply(domain_ok, fit_domain) %>%
  bind_rows() %>%
  arrange(eff12_pct)

qs2::qs_save(domain_res, file.path(OUTPUT_DIR, "domain_results.qs"))

# ---- 4.2 领域异质性的联合检验 ----
# 全样本单一模型中加入 domain × COVID 即时项的交互，检验领域间差异是否显著
dat_het <- dat %>%
  filter(domain %in% domain_ok) %>%
  mutate(domain_f = factor(domain))

m_het <- feols(
  as.formula(paste0(
    "log_duration ~ pre_trend_month + ", COVID_TERMS,
    " + i(domain_f, covid_level) | journal + submit_month_num"
  )),
  data = dat_het, cluster = ~journal
)

het_terms <- grep("domain_f", names(coef(m_het)), value = TRUE)
het_test <- if (length(het_terms) > 1) {
  w <- try(fixest::wald(m_het, keep = "domain_f", print = FALSE), silent = TRUE)
  if (inherits(w, "try-error")) {
    tibble(Test = "Domain x COVID interaction", Statistic = NA,
           df = NA, P_value = NA)
  } else {
    tibble(Test = "Joint test: domain x COVID onset interactions = 0",
           Statistic = round(w$stat, 3), df = w$df1,
           P_value = signif(w$p, 3))
  }
} else {
  tibble(Test = "Domain x COVID interaction", Statistic = NA, df = NA, P_value = NA)
}
print(het_test)

# ---- 4.3 可视化 ----
# 主图：森林图形式呈现各领域第 12 个月的累计效应，按效应大小排序。
# 不在图内写标题，标题与说明另存文本文件。
plot_dat <- domain_res %>%
  mutate(
    domain_lab = fct_reorder(domain, eff12_pct),
    sig = ifelse(eff12_p < 0.05, "p < 0.05", "n.s.")
  )

overall12 <- T3$Change_pct[T3$Months_since_onset == 12]

p_domain <- ggplot(plot_dat, aes(x = eff12_pct, y = domain_lab)) +
  geom_vline(xintercept = 0, colour = "grey40", linewidth = 0.4) +
  geom_vline(xintercept = overall12, linetype = "dashed",
             colour = "#C0392B", linewidth = 0.5) +
  geom_errorbarh(aes(xmin = eff12_lo, xmax = eff12_hi),
                 height = 0, linewidth = 0.6, colour = "#34495E") +
  geom_point(aes(fill = sig), shape = 21, size = 2.8,
             colour = "#2C3E50", stroke = 0.5) +
  scale_fill_manual(values = c("p < 0.05" = "#2874A6", "n.s." = "white"),
                    name = NULL) +
  labs(x = "Change in peer review duration at month 12 (%)", y = NULL) +
  theme_classic(base_size = 11) +
  theme(
    legend.position = "top",
    axis.text.y = element_text(size = 9),
    panel.grid.major.y = element_line(colour = "grey93", linewidth = 0.3)
  )

ggsave(file.path(OUTPUT_DIR, "Figure_domain.png"), p_domain,
       width = 7.5, height = max(4, 0.32 * nrow(plot_dat) + 1.5), dpi = 300)
ggsave(file.path(OUTPUT_DIR, "Figure_domain.pdf"), p_domain,
       width = 7.5, height = max(4, 0.32 * nrow(plot_dat) + 1.5))

# 图题与说明单独存放，便于直接粘贴进论文
fig_caption <- c(
  "Figure 1. Heterogeneity in the association between COVID-19 and peer review duration across research domains.",
  "",
  "Caption:",
  "Points show the estimated cumulative change in peer review duration (submission-to-acceptance time) 12 months after the onset of COVID-19, relative to the extrapolated pre-pandemic trend, estimated separately within each research domain. Horizontal bars denote 95% confidence intervals based on journal-clustered standard errors. Filled points indicate estimates significant at the 5% level; hollow points are not statistically distinguishable from zero. The solid vertical line marks no change; the dashed red line marks the pooled all-domain estimate. Each domain-specific model includes journal fixed effects, calendar-month fixed effects, and a pre-pandemic linear time trend. Domains with fewer than 5,000 eligible articles were excluded.",
  "",
  paste0("Pooled all-domain estimate at month 12: ",
         sprintf("%.2f%%", overall12)),
  paste0("Joint test of domain heterogeneity: ",
         ifelse(is.na(het_test$P_value[1]), "not available",
                paste0("p = ", het_test$P_value[1])))
)
writeLines(fig_caption, file.path(OUTPUT_DIR, "Figure_domain_caption.txt"))

# ---- 4.4 领域结果表（供撰写文字时取数） ----
T5 <- domain_res %>%
  transmute(
    Domain            = domain,
    N_papers          = n_papers,
    N_journals        = n_journals,
    Month12_pct       = round(eff12_pct, 2),
    Month12_CI        = sprintf("[%.2f, %.2f]", eff12_lo, eff12_hi),
    Month12_sig       = star(eff12_p),
    Month24_pct       = round(eff24_pct, 2),
    Month24_CI        = sprintf("[%.2f, %.2f]", eff24_lo, eff24_hi),
    Month24_sig       = star(eff24_p),
    Topic_gap_pct     = round(gap_pct, 2),
    Topic_gap_CI      = ifelse(is.na(gap_lo), "—",
                               sprintf("[%.2f, %.2f]", gap_lo, gap_hi)),
    Topic_gap_sig     = ifelse(is.na(gap_p), "—", star(gap_p))
  )

domain_readme <- tibble(
  Sheet = c("T5_Domain_effects", "T5b_Heterogeneity_test"),
  Suggested_place = c("Figure 1 的数值来源 / 附录表", "正文一句话引用"),
  Supports_conclusion = c(
    "结论三：不同领域受 COVID 影响程度不同",
    "结论三的统计凭据：领域间差异是否显著"
  ),
  Note_for_manuscript = c(
    "Domain-specific estimates of the cumulative change in peer review duration at months 12 and 24, relative to the extrapolated pre-pandemic trend. Each model includes journal and calendar-month fixed effects; SE clustered by journal.",
    "Joint Wald test of all domain x COVID-onset interaction terms in a pooled model."
  )
)

wb2 <- createWorkbook()
add_sheet(wb2, "README", domain_readme)
add_sheet(wb2, "T5_Domain_effects", T5)
add_sheet(wb2, "T5b_Heterogeneity_test", het_test)
saveWorkbook(wb2, file.path(OUTPUT_DIR, "Domain_results.xlsx"), overwrite = TRUE)

cat("\n=== PART 2 完成 ===\n")
cat("Figure_domain.png / .pdf\n")
cat("Figure_domain_caption.txt\n")
cat("Domain_results.xlsx\n")
print(T5)
