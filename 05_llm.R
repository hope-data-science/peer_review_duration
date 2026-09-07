rm(list = ls())

# =============================================================================
# LLM (ChatGPT) 公开发布与同行评审时长
#
# 设计（已定稿）：
#   - 分析对象：最终被接收的 article；因变量 log(1 + 投稿至接收天数)
#   - 断点：2022-11-30（ChatGPT 公开发布日）
#   - 暴露：post 指示变量 + "发布后第几月"经限制性立方样条展开
#   - 固定效应：期刊 FE（硬性要求）+ 日历月 FE（季节性）
#   - 发布前线性趋势，发布后冻结，作为反事实外推基准
#   - 标准误：聚类到期刊；附双向聚类（期刊 + 投稿月）列
#
# 本脚本只输出 3 张表 + 3 张图，编号与正文一致：
#   Table 1  规格比较 + 主回归系数      -> LLM_results.xlsx / Table1
#   Table 2  逐月累计效应              -> LLM_results.xlsx / Table2
#   Table 3  领域效应 + 异质性联合检验  -> LLM_results.xlsx / Table3_domain
#                                                          / Table3_tests
#   Figure 1 全样本逐月轨迹            -> Figure1_pooled_trajectory.png/.pdf
#   Figure 2 各领域第 13 月森林图       -> Figure2_domain_forest.png/.pdf
#   Figure 3 各领域轨迹小多图           -> Figure3_domain_panel.png/.pdf
#
# 另输出 Figure_captions.txt（三张图的图题草稿）与 run_log.txt（关键诊断）
# =============================================================================

# ============================ 0. 环境与参数 ==================================
required_pkgs <- c("dplyr", "tibble", "tidyr", "lubridate", "stringr",
                   "fixest", "splines", "openxlsx", "ggplot2", "forcats",
                   "scales", "DBI", "duckdb", "qs2", "MASS")
missing_pkgs <- required_pkgs[!vapply(required_pkgs, requireNamespace,
                                      logical(1), quietly = TRUE)]
if (length(missing_pkgs) > 0) {
  install.packages(missing_pkgs, repos = "https://cloud.r-project.org")
}
invisible(lapply(required_pkgs, library, character.only = TRUE))

OUTPUT_DIR <- "data/llm_review_final"
if (dir.exists(OUTPUT_DIR)) unlink(OUTPUT_DIR, recursive = TRUE, force = TRUE)
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

LLM_START              <- as.Date("2022-11-30")  # ChatGPT 公开发布日
DATA_END_DATE          <- as.Date("2024-12-31")  # 数据实际截止日
FOLLOWUP_DAYS          <- 365L                   # 随访缓冲
MAX_DURATION_DAYS      <- 1460L
MIN_JOURNAL_PAPERS_YR  <- 20L
MIN_JOURNAL_YEARS_PRE  <- 2L
MIN_JOURNAL_YEARS_POST <- 1L
MIN_DOMAIN_PAPERS      <- 5000L
SPLINE_DF              <- 3L   # 发布后暴露时长的样条自由度

followup_cutoff    <- DATA_END_DATE - FOLLOWUP_DAYS
max_reliable_month <- floor(as.numeric(followup_cutoff - LLM_START) / 30.44)

log_lines <- character(0)
say <- function(...) {
  msg <- paste0(...)
  cat(msg, "\n")
  log_lines <<- c(log_lines, msg)
}

say("充分随访投稿截止日：", as.character(followup_cutoff))
say("发布后可靠可观测窗口（月）：", max_reliable_month)

# ============================ 1. 读入数据 =====================================
con <- dbConnect(duckdb::duckdb(), dbdir = "data/sel_data.duckdb")
on.exit(try(dbDisconnect(con, shutdown = TRUE), silent = TRUE), add = TRUE)

dat_raw <- dbGetQuery(con, "
  SELECT DOI, received_dt, pr_duration, journal, domain
  FROM review_citation
  WHERE document_type = 'Article'
")

required_cols <- c("DOI", "received_dt", "pr_duration", "journal", "domain")
missing_cols  <- setdiff(required_cols, names(dat_raw))
if (length(missing_cols) > 0) {
  stop("数据缺少字段：", paste(missing_cols, collapse = ", "))
}

dat0 <- dat_raw %>%
  mutate(
    DOI         = as.character(DOI),
    journal     = str_squish(as.character(journal)),
    domain      = str_squish(as.character(domain)),
    received_dt = as.Date(received_dt),
    pr_duration = as.numeric(pr_duration)
  )

# ============================ 2. 样本构建 =====================================
flow <- tibble(step = "原始 article 记录", n = nrow(dat0))

dat1 <- dat0 %>%
  filter(!is.na(DOI), DOI != "", !is.na(received_dt), !is.na(pr_duration),
         !is.na(journal), journal != "")
flow <- add_row(flow, step = "核心变量完整", n = nrow(dat1))

dat2 <- dat1 %>% arrange(DOI, pr_duration) %>% distinct(DOI, .keep_all = TRUE)
flow <- add_row(flow, step = "DOI 去重", n = nrow(dat2))

dat3 <- dat2 %>% filter(pr_duration > 0, pr_duration <= MAX_DURATION_DAYS)
flow <- add_row(flow, step = paste0("时长 1-", MAX_DURATION_DAYS, " 天"),
                n = nrow(dat3))

dat4 <- dat3 %>% filter(received_dt <= followup_cutoff)
flow <- add_row(flow, step = paste0("充分随访（投稿 <= ", followup_cutoff, "）"),
                n = nrow(dat4))

dat4 <- dat4 %>%
  mutate(
    submit_month     = floor_date(received_dt, "month"),
    submit_year      = year(received_dt),
    submit_month_num = month(received_dt),
    post_llm         = as.integer(received_dt >= LLM_START),
    months_from_llm  = interval(LLM_START, received_dt) %/% months(1),
    log_duration     = log1p(pr_duration)
  )

journal_year <- dat4 %>%
  count(journal, submit_year, name = "n_papers") %>%
  mutate(
    period = ifelse(as.Date(paste0(submit_year, "-01-01")) < LLM_START,
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

# 发布前线性趋势（发布后冻结在断点水平）
dat <- dat %>%
  mutate(
    month_index     = interval(min(submit_month), submit_month) %/% months(1),
    pre_trend_month = pmin(month_index,
                           max(month_index[received_dt < LLM_START], na.rm = TRUE)),
    exposure_month  = ifelse(post_llm == 1, pmax(months_from_llm, 0), 0)
  )

say("稳定期刊数量：", length(stable_journals))
say("最终分析样本：", format(nrow(dat), big.mark = ","))
for (i in seq_len(nrow(flow))) {
  say("  [样本流] ", flow$step[i], " = ", format(flow$n[i], big.mark = ","))
}

# ===================== 3. 暴露时长的三种函数形式 ==============================
# 样条基以"发布后样本"的暴露月分位数定节点，再对全样本求值，
# 发布前各行置 0，使系数可直接解读为"相对外推发布前趋势的偏离"。
post_exposure <- dat$exposure_month[dat$post_llm == 1]
spline_knots  <- as.numeric(quantile(post_exposure,
                                     probs = seq(0, 1, length.out = SPLINE_DF + 1)[-c(1, SPLINE_DF + 1)]))
spline_bound  <- range(post_exposure)

make_basis <- function(x) {
  B <- ns(x, knots = spline_knots, Boundary.knots = spline_bound)
  B <- matrix(as.numeric(B), nrow = length(x))
  colnames(B) <- paste0("spl", seq_len(ncol(B)))
  B
}

B_all <- make_basis(dat$exposure_month)
B_all[dat$post_llm == 0, ] <- 0
spline_cols <- colnames(B_all)
dat <- bind_cols(dat, as_tibble(B_all))

dat <- dat %>%
  mutate(
    exp_linear = exposure_month,
    exp_log    = ifelse(post_llm == 1, log1p(exposure_month), 0)
  )

say("样条内部节点（暴露月）：", paste(round(spline_knots, 2), collapse = ", "))

FE      <- " | journal + submit_month_num"
BASE    <- "log_duration ~ pre_trend_month + post_llm"
F_LIN   <- paste0(BASE, " + exp_linear", FE)
F_LOG   <- paste0(BASE, " + exp_log", FE)
F_SPL   <- paste0(BASE, " + ", paste(spline_cols, collapse = " + "), FE)

m_lin <- feols(as.formula(F_LIN), data = dat, cluster = ~journal)
m_log <- feols(as.formula(F_LOG), data = dat, cluster = ~journal)
m_spl <- feols(as.formula(F_SPL), data = dat, cluster = ~journal)
m_spl_2w <- feols(as.formula(F_SPL), data = dat,
                  cluster = ~journal + submit_month)

qs2::qs_save(list(m_lin = m_lin, m_log = m_log, m_spl = m_spl,
                  m_spl_2w = m_spl_2w,
                  spline_knots = spline_knots, spline_bound = spline_bound),
             file.path(OUTPUT_DIR, "models.qs"))

# ---- 通用联合 Wald 检验 ----
# 不使用 fixest::wald()：当 keep 模式匹配不到存活系数时它返回原子 NA，
# 后续取 $stat 会直接报错。这里直接用系数向量与协方差矩阵手工计算，
# 并自动剔除因共线性被丢弃的 NA 系数，奇异时改用广义逆。
wald_joint <- function(model, terms, label = NA_character_) {
  b <- coef(model)
  b <- b[!is.na(b)]
  sel <- intersect(terms, names(b))
  if (length(sel) == 0) {
    return(tibble(Test = label, F_statistic = NA_real_, df1 = NA_integer_,
                  df2 = NA_integer_, P_value = NA_real_,
                  Note = "no estimable coefficients matched"))
  }
  bs <- b[sel]
  V  <- vcov(model)[sel, sel, drop = FALSE]
  Vi <- try(solve(V), silent = TRUE)
  singular <- inherits(Vi, "try-error")
  if (singular) Vi <- MASS::ginv(V)
  q  <- length(sel)
  W  <- as.numeric(t(bs) %*% Vi %*% bs)
  df2 <- tryCatch(fixest::degrees_freedom(model, type = "t"),
                  error = function(e) NA_integer_)
  Fstat <- W / q
  p <- if (is.na(df2)) pchisq(W, df = q, lower.tail = FALSE) else
    pf(Fstat, df1 = q, df2 = df2, lower.tail = FALSE)
  tibble(Test = label, F_statistic = round(Fstat, 3), df1 = q,
         df2 = df2, P_value = signif(p, 4),
         Note = if (singular) "generalised inverse used (singular VCOV)" else "")
}

# 按正则从模型中挑出系数名
pick_terms <- function(model, pattern) {
  b <- coef(model)
  grep(pattern, names(b[!is.na(b)]), value = TRUE)
}

# 线性约束的 Wald 检验：样条附加项是否联合为零
wald_shape <- wald_joint(m_spl, spline_cols[-1],
                         "Wald test of linear restriction")

say("AIC: linear=", round(AIC(m_lin), 1),
    " log=", round(AIC(m_log), 1),
    " spline=", round(AIC(m_spl), 1))
say("BIC: linear=", round(BIC(m_lin), 1),
    " log=", round(BIC(m_log), 1),
    " spline=", round(BIC(m_spl), 1))
say("样条附加项 Wald F=", wald_shape$F_statistic,
    " df1=", wald_shape$df1, " p=", wald_shape$P_value)

# ========================= 4. 通用推断工具 ====================================
lincom <- function(model, weights) {
  b <- coef(model); V <- vcov(model)
  L <- setNames(rep(0, length(b)), names(b))
  keep <- intersect(names(weights), names(b))
  L[keep] <- weights[keep]
  est <- sum(L * b)
  se  <- sqrt(as.numeric(t(L) %*% V %*% L))
  z   <- est / se
  tibble(log_est = est, se = se,
         pct    = 100 * (exp(est) - 1),
         lo_pct = 100 * (exp(est - 1.96 * se) - 1),
         hi_pct = 100 * (exp(est + 1.96 * se) - 1),
         p_value = 2 * pnorm(-abs(z)))
}

# 第 k 月的线性组合权重：post 指示 + 样条基在 k 处的取值
weights_at <- function(k, prefix = "") {
  bk <- as.numeric(make_basis(k))
  w  <- c(1, bk)
  names(w) <- paste0(prefix, c("post_llm", spline_cols))
  w
}

effect_at <- function(model, k, prefix = "") lincom(model, weights_at(k, prefix))

# 相邻月增量（用于判断末端是否趋稳）
increment_at <- function(model, k, prefix = "") {
  w <- weights_at(k, prefix) - weights_at(k - 1, prefix)
  lincom(model, w)
}

star <- function(p) ifelse(is.na(p), "",
                           ifelse(p < 0.01, "***",
                                  ifelse(p < 0.05, "**",
                                         ifelse(p < 0.1, "*", ""))))
fmt_coef <- function(est, se, p, digits = 4) {
  sprintf("%.*f%s (%.*f)", digits, est, star(p), digits, se)
}
coef_cell <- function(model, term, digits = 4) {
  b <- coef(model)
  if (!term %in% names(b)) return("—")
  se <- sqrt(vcov(model)[term, term])
  fmt_coef(b[[term]], se, 2 * pnorm(-abs(b[[term]] / se)), digits)
}

# ===================== 5. Table 1：规格比较 + 主回归 ==========================
term_map <- tibble(
  term  = c("pre_trend_month", "post_llm", "exp_linear", "exp_log", spline_cols),
  label = c("Pre-release monthly trend",
            "Level shift at release",
            "Exposure time, constant monthly slope",
            "Exposure time, logarithmic slope",
            paste0("Exposure-time spline basis ", seq_along(spline_cols)))
)

build_col <- function(model, colname) {
  body <- term_map %>%
    rowwise() %>%
    mutate(v = coef_cell(model, term)) %>%
    ungroup() %>%
    transmute(Variable = label, v)
  foot <- tibble(
    Variable = c("Journal fixed effects", "Calendar-month fixed effects",
                 "Standard errors clustered by",
                 "Observations", "AIC", "BIC", "Within R2"),
    v = c("Yes", "Yes",
          if (identical(colname, "(4) Spline, two-way SE")) "Journal & submission month" else "Journal",
          format(nobs(model), big.mark = ","),
          sprintf("%.1f", AIC(model)),
          sprintf("%.1f", BIC(model)),
          sprintf("%.4f", fitstat(model, "wr2")$wr2))
  )
  out <- bind_rows(body, foot)
  names(out)[2] <- colname
  out
}

Table1 <- build_col(m_lin, "(1) Linear") %>%
  left_join(build_col(m_log, "(2) Logarithmic"), by = "Variable") %>%
  left_join(build_col(m_spl, "(3) Spline"), by = "Variable") %>%
  left_join(build_col(m_spl_2w, "(4) Spline, two-way SE"), by = "Variable")

Table1 <- bind_rows(
  Table1,
  tibble(Variable = "Wald test of linear restriction (spline terms = 0)",
         `(1) Linear` = "—", `(2) Logarithmic` = "—",
         `(3) Spline` = sprintf("F = %.3f, p %s", wald_shape$F_statistic,
                                ifelse(wald_shape$P_value < 0.001, "< 0.001",
                                       paste0("= ", wald_shape$P_value))),
         `(4) Spline, two-way SE` = "—")
)

# 描述统计并入 Table 1 脚注所需的数字
desc <- dat %>%
  mutate(Period = ifelse(post_llm == 1, "Post-release", "Pre-release")) %>%
  group_by(Period) %>%
  summarise(N = n(), Mean = mean(pr_duration), SD = sd(pr_duration),
            Q1 = quantile(pr_duration, .25), Median = median(pr_duration),
            Q3 = quantile(pr_duration, .75), .groups = "drop")

desc_note <- paste0(
  "Unadjusted submission-to-acceptance days: pre-release mean ",
  sprintf("%.1f", desc$Mean[desc$Period == "Pre-release"]),
  " (SD ", sprintf("%.1f", desc$SD[desc$Period == "Pre-release"]),
  "; median ", sprintf("%.0f", desc$Median[desc$Period == "Pre-release"]),
  ", IQR ", sprintf("%.0f-%.0f", desc$Q1[desc$Period == "Pre-release"],
                    desc$Q3[desc$Period == "Pre-release"]),
  "; N = ", format(desc$N[desc$Period == "Pre-release"], big.mark = ","),
  "); post-release mean ",
  sprintf("%.1f", desc$Mean[desc$Period == "Post-release"]),
  " (SD ", sprintf("%.1f", desc$SD[desc$Period == "Post-release"]),
  "; median ", sprintf("%.0f", desc$Median[desc$Period == "Post-release"]),
  ", IQR ", sprintf("%.0f-%.0f", desc$Q1[desc$Period == "Post-release"],
                    desc$Q3[desc$Period == "Post-release"]),
  "; N = ", format(desc$N[desc$Period == "Post-release"], big.mark = ","), ")."
)
say(desc_note)

# ===================== 6. Table 2 + Figure 1：逐月效应 ========================
months_seq <- 0:max_reliable_month

Table2 <- lapply(months_seq, function(k) {
  e <- effect_at(m_spl, k)
  d <- if (k == 0) tibble(inc_pct = NA_real_, inc_p = NA_real_) else {
    ii <- increment_at(m_spl, k)
    tibble(inc_pct = ii$pct, inc_p = ii$p_value)
  }
  bind_cols(tibble(Month_since_release = k), e, d)
}) %>% bind_rows()

Table2_out <- Table2 %>%
  transmute(
    Month_since_release,
    Change_pct    = round(pct, 2),
    CI_95         = sprintf("[%.2f, %.2f]", lo_pct, hi_pct),
    P_value       = signif(p_value, 3),
    Sig           = star(p_value),
    Increment_pct = ifelse(is.na(inc_pct), NA, round(inc_pct, 2)),
    Increment_sig = star(inc_p)
  )

end_eff <- Table2 %>% filter(Month_since_release == max_reliable_month)
m12_eff <- Table2 %>% filter(Month_since_release == 12)
end_inc <- Table2 %>% filter(Month_since_release == max_reliable_month)
first_sig <- Table2 %>% filter(p_value < 0.05, pct > 0) %>%
  slice_min(Month_since_release, n = 1, with_ties = FALSE)

say("第 12 月累计效应：", sprintf("%.2f%%", m12_eff$pct))
say("第 ", max_reliable_month, " 月累计效应：",
    sprintf("%.2f%% [%.2f, %.2f]", end_eff$pct, end_eff$lo_pct, end_eff$hi_pct))
if (nrow(first_sig) == 1) {
  say("首个显著为正的月份：month ", first_sig$Month_since_release)
}
say("末月环比增量：", sprintf("%.2f%% (p = %s)", end_inc$inc_pct,
                       signif(end_inc$inc_p, 3)))

fig1 <- ggplot(Table2, aes(x = Month_since_release, y = pct)) +
  geom_hline(yintercept = 0, colour = "grey40", linewidth = 0.4) +
  geom_ribbon(aes(ymin = lo_pct, ymax = hi_pct), fill = "#2874A6",
              alpha = 0.18) +
  geom_line(colour = "#1B4F72", linewidth = 0.9) +
  geom_point(colour = "#1B4F72", size = 1.8) +
  scale_x_continuous(breaks = months_seq) +
  labs(x = "Months since public release of ChatGPT (30 November 2022)",
       y = "Deviation in peer review duration\nfrom extrapolated pre-release trend (%)") +
  theme_classic(base_size = 11) +
  theme(panel.grid.major.y = element_line(colour = "grey93", linewidth = 0.3))

ggsave(file.path(OUTPUT_DIR, "Figure1_pooled_trajectory.png"), fig1,
       width = 7.2, height = 4.4, dpi = 300)
ggsave(file.path(OUTPUT_DIR, "Figure1_pooled_trajectory.pdf"), fig1,
       width = 7.2, height = 4.4)

# =============== 7. Table 3 + Figure 2/3：领域异质性 ==========================
domain_ok <- dat %>%
  filter(!is.na(domain), domain != "") %>%
  count(domain, name = "n") %>%
  filter(n >= MIN_DOMAIN_PAPERS) %>%
  arrange(desc(n)) %>%
  pull(domain)

say("纳入领域分析的 domain 数量：", length(domain_ok))

dat_het <- dat %>% filter(domain %in% domain_ok) %>%
  mutate(domain_f = factor(domain))

# 分领域独立估计（共用同一组样条节点，保证跨领域可比）
fit_domain <- function(d) {
  sub <- dat_het %>% filter(domain == d)
  m <- try(feols(as.formula(F_SPL), data = sub, cluster = ~journal),
           silent = TRUE)
  if (inherits(m, "try-error")) {
    say("  [跳过] ", d, "：模型未收敛")
    return(NULL)
  }
  e6  <- effect_at(m, 6)
  eEnd <- effect_at(m, max_reliable_month)
  traj <- lapply(months_seq, function(k) {
    effect_at(m, k) %>% mutate(domain = d, month = k, .before = 1)
  }) %>% bind_rows()
  list(
    summary = tibble(
      domain     = d,
      n_papers   = nrow(sub),
      n_journals = n_distinct(sub$journal),
      level_shift_pct = 100 * (exp(coef(m)[["post_llm"]]) - 1),
      m6_pct = e6$pct, m6_lo = e6$lo_pct, m6_hi = e6$hi_pct, m6_p = e6$p_value,
      end_pct = eEnd$pct, end_lo = eEnd$lo_pct, end_hi = eEnd$hi_pct,
      end_p = eEnd$p_value
    ),
    traj = traj
  )
}

dom_fits <- lapply(domain_ok, fit_domain)
dom_fits <- dom_fits[!vapply(dom_fits, is.null, logical(1))]

domain_res <- bind_rows(lapply(dom_fits, `[[`, "summary")) %>%
  arrange(desc(end_pct))
domain_traj <- bind_rows(lapply(dom_fits, `[[`, "traj"))

qs2::qs_save(list(domain_res = domain_res, domain_traj = domain_traj),
             file.path(OUTPUT_DIR, "domain_results.qs"))

# --- 异质性联合检验：手工构造交互列，避免嵌套 i() 报错 ---
het_terms <- c("post_llm", spline_cols)
F_HET <- paste0(
  "log_duration ~ pre_trend_month + ", paste(het_terms, collapse = " + "),
  " + ", paste(paste0("domain_f:", het_terms), collapse = " + "), FE
)
m_het <- feols(as.formula(F_HET), data = dat_het, cluster = ~journal)

# 交互项系数名形如 "domain_fMedicine:post_llm" / "post_llm:domain_fMedicine"，
# 两种顺序都要能匹配到。
terms_all   <- pick_terms(m_het, "domain_f")
terms_level <- terms_all[grepl("post_llm", terms_all)]
terms_shape <- terms_all[grepl(paste(spline_cols, collapse = "|"), terms_all)]

say("交互项系数个数：合计 ", length(terms_all),
    "；水平跳变 ", length(terms_level),
    "；轨迹形状 ", length(terms_shape))

Table3_tests <- bind_rows(
  wald_joint(m_het, terms_all,
             "Joint test: all domain-by-exposure interactions = 0"),
  wald_joint(m_het, terms_level,
             "Joint test: domain-specific level shifts at release = 0"),
  wald_joint(m_het, terms_shape,
             "Joint test: domain-specific trajectory shapes identical")
)
say("领域异质性联合检验 F = ",
    ifelse(is.na(Table3_tests$F_statistic[1]), "NA",
           Table3_tests$F_statistic[1]),
    "，p = ", ifelse(is.na(Table3_tests$P_value[1]), "NA",
                    Table3_tests$P_value[1]))

Table3_domain <- domain_res %>%
  transmute(
    Domain          = domain,
    N_articles      = n_papers,
    N_journals      = n_journals,
    Level_shift_pct = round(level_shift_pct, 2),
    Month6_pct      = round(m6_pct, 2),
    Month6_CI       = sprintf("[%.2f, %.2f]", m6_lo, m6_hi),
    Month6_sig      = star(m6_p),
    Month_end_pct   = round(end_pct, 2),
    Month_end_CI    = sprintf("[%.2f, %.2f]", end_lo, end_hi),
    Month_end_p     = signif(end_p, 3),
    Month_end_sig   = star(end_p)
  )
names(Table3_domain)[names(Table3_domain) == "Month_end_pct"] <-
  paste0("Month", max_reliable_month, "_pct")
names(Table3_domain)[names(Table3_domain) == "Month_end_CI"] <-
  paste0("Month", max_reliable_month, "_CI")
names(Table3_domain)[names(Table3_domain) == "Month_end_p"] <-
  paste0("Month", max_reliable_month, "_p")
names(Table3_domain)[names(Table3_domain) == "Month_end_sig"] <-
  paste0("Month", max_reliable_month, "_sig")

# --- Figure 2：森林图 ---
plot_forest <- domain_res %>%
  mutate(domain_lab = fct_reorder(domain, end_pct),
         sig = ifelse(end_p < 0.05, "p < 0.05", "n.s."))

fig2 <- ggplot(plot_forest, aes(x = end_pct, y = domain_lab)) +
  geom_vline(xintercept = 0, colour = "grey40", linewidth = 0.4) +
  geom_vline(xintercept = end_eff$pct, linetype = "dashed",
             colour = "#C0392B", linewidth = 0.5) +
  geom_errorbarh(aes(xmin = end_lo, xmax = end_hi), height = 0,
                 linewidth = 0.6, colour = "#34495E") +
  geom_point(aes(fill = sig), shape = 21, size = 2.9,
             colour = "#2C3E50", stroke = 0.5) +
  scale_fill_manual(values = c("p < 0.05" = "#2874A6", "n.s." = "white"),
                    name = NULL) +
  labs(x = paste0("Deviation in peer review duration at month ",
                  max_reliable_month, " (%)"), y = NULL) +
  theme_classic(base_size = 11) +
  theme(legend.position = "top", axis.text.y = element_text(size = 9),
        panel.grid.major.y = element_line(colour = "grey93", linewidth = 0.3))

ggsave(file.path(OUTPUT_DIR, "Figure2_domain_forest.png"), fig2,
       width = 7.5, height = max(4, 0.34 * nrow(plot_forest) + 1.5), dpi = 300)
ggsave(file.path(OUTPUT_DIR, "Figure2_domain_forest.pdf"), fig2,
       width = 7.5, height = max(4, 0.34 * nrow(plot_forest) + 1.5))

# --- Figure 3：分领域轨迹小多图 ---
traj_plot <- domain_traj %>%
  mutate(domain = factor(domain, levels = domain_res$domain))

fig3 <- ggplot(traj_plot, aes(x = month, y = pct)) +
  geom_hline(yintercept = 0, colour = "grey45", linewidth = 0.35) +
  geom_ribbon(aes(ymin = lo_pct, ymax = hi_pct), fill = "#2874A6",
              alpha = 0.16) +
  geom_line(colour = "#1B4F72", linewidth = 0.7) +
  facet_wrap(~ domain, ncol = 3, scales = "free_y") +
  scale_x_continuous(breaks = seq(0, max_reliable_month, by = 3)) +
  labs(x = "Months since public release of ChatGPT",
       y = "Deviation from extrapolated pre-release trend (%)") +
  theme_bw(base_size = 10) +
  theme(panel.grid.minor = element_blank(),
        strip.background = element_rect(fill = "grey95", colour = NA),
        strip.text = element_text(face = "bold", size = 9))

n_rows_fig3 <- ceiling(nrow(domain_res) / 3)
ggsave(file.path(OUTPUT_DIR, "Figure3_domain_panel.png"), fig3,
       width = 9, height = 2.3 * n_rows_fig3 + 0.8, dpi = 300)
ggsave(file.path(OUTPUT_DIR, "Figure3_domain_panel.pdf"), fig3,
       width = 9, height = 2.3 * n_rows_fig3 + 0.8)

# ========================= 8. 写出 Excel =====================================
readme <- tibble(
  Output = c("Table1", "Table2", "Table3_domain", "Table3_tests",
             "Figure1_pooled_trajectory.png", "Figure2_domain_forest.png",
             "Figure3_domain_panel.png"),
  Manuscript_label = c("Table 1", "Table 2", "Table 3 (domain estimates)",
                       "Table 3 (joint tests, reported in the notes)",
                       "Figure 1", "Figure 2", "Figure 3"),
  Description = c(
    "Functional-form comparison and main regression coefficients across four specifications.",
    "Month-by-month cumulative deviation from the extrapolated pre-release trend, with month-over-month increments.",
    "Domain-specific level shift, month-6 and terminal-month deviations with 95% CIs.",
    "Wald tests of domain heterogeneity: overall, level shifts, and trajectory shapes.",
    "Pooled monthly trajectory with 95% confidence band.",
    "Forest plot of domain estimates at the terminal month; dashed line is the pooled estimate.",
    "Small-multiple domain trajectories, free vertical scales."
  )
)

key_numbers <- tibble(
  Item = c("Analytic sample (articles)", "Stable journals",
           "Reliable follow-up window (months)",
           "Unadjusted pre-release mean days", "Unadjusted post-release mean days",
           "Deviation at month 12 (%)",
           paste0("Deviation at month ", max_reliable_month, " (%)"),
           paste0("Month-over-month increment at month ", max_reliable_month, " (%)"),
           "First month with significant positive deviation",
           "Spline vs linear Wald F", "Spline vs linear p",
           "Domain heterogeneity F", "Domain heterogeneity p",
           "Descriptive statistics (for Table 1 notes)"),
  Value = c(format(nrow(dat), big.mark = ","),
            length(stable_journals),
            max_reliable_month,
            sprintf("%.1f", desc$Mean[desc$Period == "Pre-release"]),
            sprintf("%.1f", desc$Mean[desc$Period == "Post-release"]),
            sprintf("%.2f", m12_eff$pct),
            sprintf("%.2f", end_eff$pct),
            sprintf("%.2f (p = %s)", end_inc$inc_pct, signif(end_inc$inc_p, 3)),
            ifelse(nrow(first_sig) == 1, first_sig$Month_since_release, "none"),
            sprintf("%.3f", wald_shape$F_statistic),
            ifelse(wald_shape$P_value < 0.001, "< 0.001", wald_shape$P_value),
            ifelse(is.na(Table3_tests$F_statistic[1]), "NA",
                   Table3_tests$F_statistic[1]),
            ifelse(is.na(Table3_tests$P_value[1]), "NA",
                   Table3_tests$P_value[1]),
            desc_note)
)

hs <- createStyle(textDecoration = "bold", halign = "center",
                  fgFill = "#F2F2F2", border = "TopBottom")
add_sheet <- function(wb, name, df) {
  addWorksheet(wb, name)
  writeData(wb, name, df, headerStyle = hs)
  freezePane(wb, name, firstRow = TRUE)
  setColWidths(wb, name, cols = seq_len(ncol(df)), widths = "auto")
}

wb <- createWorkbook()
add_sheet(wb, "README", readme)
add_sheet(wb, "Key_numbers", key_numbers)
add_sheet(wb, "Table1", Table1)
add_sheet(wb, "Table2", Table2_out)
add_sheet(wb, "Table3_domain", Table3_domain)
add_sheet(wb, "Table3_tests", Table3_tests)
saveWorkbook(wb, file.path(OUTPUT_DIR, "LLM_results.xlsx"), overwrite = TRUE)

# ========================= 9. 图题草稿 =======================================
captions <- c(
  paste0("Figure 1. Peer review duration relative to the extrapolated ",
         "pre-release trend, by month since the public release of ChatGPT."),
  paste0("Notes. The line plots the estimated cumulative deviation in ",
         "submission-to-acceptance time from the trend the pre-release period ",
         "was already following; the shaded band is the 95% confidence ",
         "interval from journal-clustered standard errors. Estimates come from ",
         "a model with journal and calendar-month fixed effects, a pre-release ",
         "linear trend frozen at the release date, an indicator for submission ",
         "after 30 November 2022, and post-release exposure time entered as a ",
         "restricted cubic spline. The window ends at month ",
         max_reliable_month, " because the underlying data extend only to ",
         "31 December 2024 and a 365-day follow-up buffer is retained."),
  "",
  paste0("Figure 2. Domain-specific change in peer review duration at month ",
         max_reliable_month, " after the public release of ChatGPT."),
  paste0("Notes. Points are domain-specific estimates from separate models ",
         "sharing a common spline basis; horizontal bars are 95% confidence ",
         "intervals from journal-clustered standard errors. Filled points are ",
         "significant at the 5% level, hollow points are not. The dashed line ",
         "marks the pooled all-domain estimate (",
         sprintf("%.2f", end_eff$pct), "%); the solid line marks no change. ",
         "Domains with fewer than ",
         format(MIN_DOMAIN_PAPERS, big.mark = ","),
         " eligible articles were excluded. A non-significant estimate ",
         "indicates limited precision rather than an established absence of ",
         "effect."),
  "",
  "Figure 3. Domain-specific trajectories of peer review duration.",
  paste0("Notes. Each panel plots one domain's estimated cumulative deviation ",
         "from its own extrapolated pre-release trend, with 95% confidence ",
         "bands. Vertical scales differ across panels to make the shape of ",
         "each trajectory legible; magnitudes should therefore be compared ",
         "using Figure 2 or Table 3 rather than by eye across panels.")
)
writeLines(captions, file.path(OUTPUT_DIR, "Figure_captions.txt"))
writeLines(log_lines, file.path(OUTPUT_DIR, "run_log.txt"))

# ========================= 10. 输出清单 ======================================
cat("\n================ 运行完成，请回传以下文件 ================\n")
cat("目录：", normalizePath(OUTPUT_DIR), "\n\n")
cat("1) LLM_results.xlsx                （含 Table1 / Table2 / Table3_domain / Table3_tests / Key_numbers）\n")
cat("2) Figure1_pooled_trajectory.png\n")
cat("3) Figure2_domain_forest.png\n")
cat("4) Figure3_domain_panel.png\n")
cat("5) Figure_captions.txt\n")
cat("6) run_log.txt\n")
cat("\n（.pdf 版与 .qs 模型对象也已保存，撰写 Results 时不需要回传。）\n")
print(Table2_out)
print(Table3_domain)
