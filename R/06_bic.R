# Compare model performances for the SPD of each dating method using BIC

library(ggplot2)
library(patchwork)

# Load models fit to radiocarbon SPD
load(file = here::here("analysis/data/derived_data/rcarbon_models.RData"))

# Find BIC values
BIC.exp <- 1*log(134) - 2*(-exp$value)
BIC.log <- 2*log(134) - 2*(-log$value)
BIC.uniform <- 0 - 2*(-uniform)

# Assemble data frame for plotting
bics <- data.frame(model = c("Exponential", "Logistic", "Uniform"),
                   bic = c(BIC.exp, BIC.log, BIC.uniform))

# Create plot
bicplt_rcarbon <- ggplot(data = bics, aes(x = bic, y = model)) +
  geom_point(size = 2.5) +
  labs(x = "BIC", y = "", subtitle = "Models fit to RSPD") +
  theme_bw()

# Repeat for SPD of shoreline dates
load(file = here::here("analysis/data/derived_data/shore_models.RData"))

BIC.exp <- 1*log(921) - 2*(-exp$value)
BIC.log <- 2*log(921) - 2*(-logi$value)
BIC.uniform <- 0 - 2*(unif)

bics <- data.frame(model = c("Exponential", "Logistic", "Uniform"),
                   bic = c(BIC.exp, BIC.log, BIC.uniform))

bicplt_shore <- ggplot(data = bics, aes(x = bic, y = model)) +
  geom_point(size = 2.5) +
  labs(x = "BIC", y = "", subtitle = "Models fit to SSPD") +
  theme_bw()

# Combine plots
bicplt_shore + bicplt_rcarbon

# Save plot
ggsave(filename = here::here("analysis/figures/bic.png"),
       units = "px", width = 2000, height = 900)
