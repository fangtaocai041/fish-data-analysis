# ============================================================
# G — 时间序列与动态分析 R 代码模板
# 方陶文库 · 数量生态学
# ============================================================

library(tidyverse)
library(forecast)     # ARIMA, STL
library(mgcv)         # GAM 时间平滑
library(changepoint)  # 变化点检测
library(earlywarnings) # 早期预警信号
library(ggplot2)

# ============================================================
# 0. 模拟数据
# ============================================================
set.seed(42)
t <- 1:100

# 模拟种群动态：趋势 + 季节 + 噪声 + 临界转换
trend <- 50 + 0.33 * t
seasonal <- 10 * sin(2 * pi * t / 12)
noise <- rnorm(100, 0, 3)

# 在 t=70 处加入状态转换
shift <- c(rep(0, 70), seq(0, -30, length.out = 30))
pop <- trend + seasonal + noise + shift

ts_data <- data.frame(
  time = t,
  abundance = pop,
  date = seq.Date(from = as.Date("2010-01-01"), by = "month", length.out = 100)
)

ggplot(ts_data, aes(x = date, y = abundance)) +
  geom_line(color = "steelblue", size = 1) +
  geom_smooth(method = "loess", se = TRUE, color = "coral") +
  labs(title = "模拟种群动态", y = "丰度") + theme_bw()

# ============================================================
# G1 — 时间序列分解
# ============================================================

# --- 1a. 转为 ts 对象 ---
pop_ts <- ts(ts_data$abundance, start = c(2010, 1), frequency = 12)

# --- 1b. STL 分解 ---
stl_decomp <- stl(pop_ts, s.window = "periodic")
plot(stl_decomp, main = "STL 分解：趋势 + 季节 + 残差")

# 提取分量
trend_component <- as.numeric(stl_decomp$time.series[, "trend"])
seasonal_component <- as.numeric(stl_decomp$time.series[, "seasonal"])
residual_component <- as.numeric(stl_decomp$time.series[, "remainder"])

# --- 1c. 变化点检测 ---
cpt_mean <- cpt.mean(pop_ts, method = "PELT")
plot(cpt_mean, type = "l", cpt.col = "red",
     main = "均值变化点检测")
cpts(cpt_mean)  # 变化点位置

cpt_var <- cpt.var(pop_ts, method = "PELT")
plot(cpt_var, main = "方差变化点检测")

# ============================================================
# G2 — 时间序列建模
# ============================================================

# --- 2a. ARIMA ---
fit_arima <- auto.arima(pop_ts, seasonal = TRUE)
summary(fit_arima)

# 诊断残差
checkresiduals(fit_arima)

# 预测 12 个月
forecast_arima <- forecast(fit_arima, h = 12)
autoplot(forecast_arima) + labs(title = "ARIMA 预测")

# --- 2b. GAM 时间平滑 ---
gam_time <- gam(abundance ~ s(time, k = 20) + s(cos(2*pi*time/12), k = 5),
                data = ts_data)
ts_data$gam_fitted <- fitted(gam_time)

ggplot(ts_data, aes(x = time)) +
  geom_line(aes(y = abundance), alpha = 0.5) +
  geom_line(aes(y = gam_fitted), color = "red", size = 1) +
  labs(title = "GAM 拟合（非线性趋势 + 季节）") + theme_bw()

# ============================================================
# G3 — 早期预警信号
# ============================================================

# --- 3a. 计算 EWS ---
# 在最后 50 个时间点滑动窗口计算
ews_results <- generic_ews(
  pop_ts,
  winsize = 50,
  detrending = "loess",
  bandwidth = 10,
  logtransform = FALSE,
  interpolate = FALSE,
  AR_n = FALSE,
  powerspectrum = FALSE
)
ews_results  # 包含方差、自相关(AR1)、偏度、峰度等指标

# --- 3b. 可视化 EWS ---
# AR1 和方差的变化趋势
ggplot(ews_results, aes(x = timeindex)) +
  geom_line(aes(y = ar1, color = "AR1")) +
  geom_line(aes(y = variance, color = "Variance")) +
  labs(title = "早期预警信号：AR1 和方差", y = "值") +
  scale_color_manual(values = c("AR1" = "coral", "Variance" = "steelblue")) +
  theme_bw()

# --- 3c. Kendall tau 趋势检验 ---
# 检验 EWS 是否随时间显著增加
cor.test(ews_results$timeindex, ews_results$ar1, method = "kendall")
cor.test(ews_results$timeindex, ews_results$variance, method = "kendall")
