#for 5000 iterations
Uniform <- readRDS("path/Results_uniform.rds")
Skewed <- readRDS("path/Results_uniform.rds")


add_dist <- function(df, dist_label) df %>% mutate(dist = dist_label, log_RI_true = log(RI_true))


all_res <- bind_rows(
  add_dist(Skewed,   "Skewed"),
  add_dist(Uniform, "Uniform")
)

# --- 2) Performance measures on the LOG scale ---
perf1 <- all_res %>%
  group_by(dist, scenario, years, RI_true, variation, approach, dose) %>%
  summarise(
    n                = n(),
    
    # --- point estimates ---
    mean_log         = mean(log_RI_est),
    median_log       = median(log_RI_est),
    emp_se_log       = sd(log_RI_est),
    avg_model_se_log = mean(SE_est),
    bias_mean_log    = mean_log   - first(log_RI_true),
    bias_median_log  = median_log - first(log_RI_true),
    mse_log          = mean((log_RI_est - first(log_RI_true))^2),
    cover95          = mean({
      lo <- log_RI_est - 1.96 * SE_est
      hi <- log_RI_est + 1.96 * SE_est
      (first(log_RI_true) >= lo) & (first(log_RI_true) <= hi)
    }),
    
    # --- Monte Carlo standard errors ---
    mcse_bias_mean    = emp_se_log / sqrt(n),
    
    # approximate under normality
    mcse_bias_median  = 1.2533 * emp_se_log / sqrt(n),
    
    mcse_emp_se       = emp_se_log / sqrt(2 * (n - 1)),
    
    mcse_avg_model_se = sd(SE_est) / sqrt(n),
    
    mcse_mse          = sqrt(
      var((log_RI_est - first(log_RI_true))^2) / n
    ),
    
    mcse_rmse         = ifelse(
      mse_log > 0,
      mcse_mse / (2 * sqrt(mse_log)),
      NA_real_
    ),
    
    mcse_cover95      = sqrt(cover95 * (1 - cover95) / n),
    
    .groups = "drop"
  )


# --- 3) Bias plots 
y_rng <- range(perf1$bias_median_log, na.rm = TRUE)
y_pad <- ifelse(diff(y_rng) == 0, 0.1, 0.08 * diff(y_rng))
ylim_use <- c(y_rng[1] - y_pad, y_rng[2] + y_pad)

brks_x <- sort(unique(all_res$variation))


# Skewed
p_bias_skew <- perf1 %>%
  filter(dist == "Skewed") %>%
  ggplot(aes(x = variation, y = bias_median_log, colour = approach, group = approach)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_line() +
  geom_point(size = 1.6) +
  scale_x_continuous(breaks = brks_x) +
  coord_cartesian(ylim = ylim_use) +
  labs(title = "Bias of the estimate: Left-skewed vaccine distribution",
       x = "Variation ID (0 = base)", y = "Median(log RI) − true(log RI)", colour = "Method") +
  facet_grid(dose ~ scenario) +
  theme_minimal(base_size = 13) +
  theme(panel.grid.minor = element_blank(),
        legend.position = "bottom")

# Uniform 
p_bias_unif <- perf1 %>%
  filter(dist == "Uniform") %>%
  ggplot(aes(x = variation, y = bias_median_log, colour = approach, group = approach)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_line() +
  geom_point(size = 1.6, shape = 17) +
  scale_x_continuous(breaks = brks_x) +
  coord_cartesian(ylim = ylim_use) +
  labs(title = "Bias of the estimate: Uniform vaccine distribution",
       x = "Variation ID (0 = base)", y = "Median(log RI) − true(logRI)", colour = "Method") +
  facet_grid(dose ~ scenario) +
  theme_minimal(base_size = 13) +
  theme(panel.grid.minor = element_blank(),
        legend.position = "bottom")

(p_bias_unif | p_bias_skew) + plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

