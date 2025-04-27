
# These codes are used to generate the figures for the manuscript:
# "Atmospheric deposition enhances marine methane production and emissions from global oceans"
# Zhuang and Mao et al., 2025 (under review in Nature Communications).
#
# The raw data for figure generation are available in the manuscript itself or from Figshare:
# https://doi.org/10.6084/m9.figshare.23912640.v1
#
# The main software and tools used include:
# - R (data analysis and plotting)
# - Illustrator (figure panel arrangement and final composition)
# - Ocean Data View (ODV) (for mapping and oceanographic data visualization)
#
# All individual plots were created in R and combined into composite figures using Illustrator.

# Author: [Shihai Mao]
# Date: [2024-11-22]

#######################################################################################################################################

# Figure 1 is composed of six panels that were combined in Illustrator. The individual plots were created using R. 
# You can find the raw data for these figures in Figshare at [https://doi.org/10.6084/m9.figshare.23912640.v1].

#######################################################################################################################################

# Set the working directory
setwd("D:/D/R Languese")

# Read data
mydata <- read.csv(file.choose(), header = TRUE)

colnames(mydata) <- c("Time", "Group", "Mean", "SD")

# Extract the day number from "Time" (remove " T " and " day")
mydata$Day <- gsub("T | day", "", mydata$Time)

# Convert "Day" to factor to control the order on x-axis
mydata$Day <- factor(mydata$Day, levels = c("0", "0.25", "1", "3", "6", "10"))

library(ggplot2)

# Plot
ggplot(mydata, aes(x = Day, y = Mean, group = Group)) +
  geom_errorbar(aes(ymin = Mean - SD, ymax = Mean + SD), colour = "black", width = 0.1, size = 0.6) +
  geom_line(aes(color = "black"), size = 0.5) +
  geom_point(aes(fill = Group), shape = 21, size = 5, color = "black", stroke = 1) +
  scale_fill_manual(values = c("Aerosol" = "yellow", 
                               "MPn" = "lightblue1", 
                               "Aerosol+MPn" = "orangered", 
                               "MPn+Pi" = "gray88", 
                               "Control" = "white")) +
  scale_color_manual(values = c("Aerosol" = "yellow", 
                                "MPn" = "lightblue1", 
                                "Aerosol+MPn" = "orangered", 
                                "MPn+Pi" = "gray88", 
                                "Control" = "black")) +
  xlab("Incubation time (day)") +
  ylab(expression("BP (nmol C L"^-1 * " d"^-1 * ")")) +
  theme_bw() +
  theme(axis.title.x = element_text(size = rel(1.3)),
        axis.title.y = element_text(size = rel(1.3)),
        axis.text.x = element_text(size = rel(1.6)),
        axis.text.y = element_text(size = rel(1.6)))

ggsave("figure1.png", width = 6, height = 4, dpi = 600)

#######################################################################################################################################

# Figure 2 is composed of five panels that were combined in Illustrator. The individual plots were created using R. 
# You can find the raw data for these figures in Figshare at [https://doi.org/10.6084/m9.figshare.23912640.v1].

#######################################################################################################################################
# 1. Panels a-d: Experimental results

# Create dataframe for panels (a-d)
mydata <- data.frame(
  Concentration = rep(c("0.1 ¦ÌM", "0.5 ¦ÌM", "1 ¦ÌM", "5 ¦ÌM"), each = 7),
  NP = rep(c(1, 5, 10, 16, 20, 25, 32), times = 4),
  Mean = c(10.08, 14.81, 24.81, 50.72, 57.32, 63.40, 64.85, 
           3.34, 11.59, 23.76, 51.27, 55.43, 59.07, 64.83,
           2.85, 10.44, 28.90, 53.04, 49.13, 56.99, 60.95,
           0.97, 14.70, 35.22, 54.35, 51.51, 56.20, 62.54),
  sd = c(3.93, 3.64, 2.67, 3.62, 3.76, 1.76, 1.62,
         0.21, 0.78, 0.94, 2.06, 2.73, 2.80, 2.96,
         0.82, 0.62, 1.94, 4.77, 2.00, 1.51, 4.89,
         0.30, 1.59, 0.95, 3.26, 2.73, 13.78, 2.66)
)

library(ggplot2)

# Set NP as factor
mydata$NP <- factor(mydata$NP, levels = c(1, 5, 10, 16, 20, 25, 32))

# Plot for 0.1 ¦ÌM
ggplot(subset(mydata, Concentration == "0.1 ¦ÌM"), aes(x = NP, y = Mean, group = NP)) +
  geom_errorbar(aes(ymax = Mean + sd, ymin = Mean - sd), colour = "black", width = 0.1, size = 0.6) +
  geom_line(size = 0.5, color = "black") +
  geom_point(aes(fill = NP), shape = 23, size = 5, color = "black", stroke = 1) +
  scale_fill_manual(values = c("1" = "white", "5" = "#D5E8F3", "10" = "skyblue", 
                               "16" = "#9183E6", "20" = "yellow", "25" = "orange", 
                               "32" = "red1")) +
  xlab("The ratio of N:P") +
  ylab("MPn Conversion Ratio (%)") +
  scale_y_continuous(limits = c(0, 72), breaks = seq(0, 72, by = 18)) +
  theme_bw() +
  theme(
    axis.title.x = element_text(size = rel(1.3)),
    axis.title.y = element_text(size = rel(1.3)),
    axis.text.x = element_text(size = rel(1.6)),
    axis.text.y = element_text(size = rel(1.6)),
    plot.title = element_text(size = 18),
    legend.position = "none"
  ) +
  ggtitle("a")

ggsave("MPn_conversion_ratio_0.1uM_fixed.png", width = 5, height = 3, dpi = 300)

############################################################################################################################
# 2. Panel e: Sigmoid curve fitting

library(minpack.lm)

# Experimental data
x <- c(1, 5, 10, 16, 20, 25, 32)
y_obs <- c(4.31, 12.89, 28.17, 52.35, 53.35, 58.92, 63.29)

# Initial values
a_init <- min(y_obs)
b_init <- max(y_obs) - min(y_obs)
c_init <- 1
d_init <- median(x)

# Nonlinear fitting
fit <- nlsLM(y_obs ~ a + b / (1 + exp(-c * (x - d))),
             start = list(a = a_init, b = b_init, c = c_init, d = d_init))

# Extract fitted parameters
params <- coef(fit)
a_fit <- params['a']
b_fit <- params['b']
c_fit <- params['c']
d_fit <- params['d']

# Calculate R2
y_fit <- predict(fit, newdata = data.frame(x = x))
y_mean <- mean(y_obs)
ss_total <- sum((y_obs - y_mean)^2)
ss_res <- sum((y_obs - y_fit)^2)
r_squared <- 1 - (ss_res / ss_total)

# Build dataframe
data <- data.frame(x = x, y_obs = y_obs)
data_fit <- data.frame(x = seq(min(x), max(x), length.out = 400))
data_fit$y_fit <- predict(fit, newdata = data_fit)

# Plot
ggplot() +
  geom_point(data = data, aes(x = x, y = y_obs, fill = factor(x)), 
             colour = "black", stroke = 1, size = 5, shape = 23) +
  scale_fill_manual(values = c("1" = "white", "5" = "#D5E8F3", "10" = "skyblue", 
                               "16" = "#9183E6", "20" = "yellow", "25" = "orange", "32" = "red1")) +
  geom_line(data = data_fit, aes(x = x, y = y_fit), color = 'red', size = 1.2) +
  xlab('The ratio of N:P') + 
  ylab('MPn conversion ratio (%)') +
  scale_y_continuous(limits = c(0, 72), breaks = seq(0, 72, by = 18)) +
  scale_x_continuous(limits = c(0, 35)) +
  ggtitle('e') +
  
  # Manual legend
  geom_segment(aes(x = 22, xend = 24, y = 10, yend = 10), color = "red", size = 1) +
  geom_text(aes(x = 24.5, y = 10, label = "Fitted Curve"), hjust = 0, vjust = 0.5, size = 4) +
  geom_point(aes(x = 22.8, y = 3), shape = 23, size = 4, fill = "white", colour = "black", stroke = 1) +
  geom_text(aes(x = 24.5, y = 3, label = "Observed Data"), hjust = 0, vjust = 0.5, size = 4) +
  
  # Formula and R2
  geom_text(aes(x = 1, y = 70, 
                label = paste0(
                  "f(N:P) = ", round(a_fit,2)," + ", 
                  round(b_fit,2)," / (1 + EXP(-", round(c_fit,2)," ¡Á (N:P - ", round(d_fit,2),")))\n",
                  "R2 = ", round(r_squared, 4)
                )), 
            hjust = 0, vjust = 1, size = 3.5, family = "serif") +
  theme_bw() +
  theme(
    plot.title = element_text(size = 18),
    legend.position = "none",
    axis.title.x = element_text(size = rel(1.3)),
    axis.title.y = element_text(size = rel(1.3)),
    axis.text.x = element_text(size = rel(1.6)),
    axis.text.y = element_text(size = rel(1.6))
  )

ggsave("MPn_conversion_curve_e.png", width = 5, height = 4, dpi = 300)


#######################################################################################################################################

# Figure 3 is composed of two panels that were combined in Illustrator. The panel a was created using Ocean Data View (V5.6.3) and the panel b was generated using R.

# You can find the raw data for these figures in Figshare at [https://doi.org/10.6084/m9.figshare.23912640.v1].

#######################################################################################################################################

library(ggplot2)
library(reshape2)

# Import data
data <- read.csv("D:/C/Desktop/Done.csv")

# Reshape data to long format for ggplot
data_long <- reshape2::melt(data, id.vars = c("Lon", "Lat"), measure.vars = c("phnJ", "Dep.MLD"))

# Create boxplot
ggplot(data_long, aes(x = value, y = variable, fill = variable)) +
  stat_boxplot(geom = "errorbar", width = 0.15) +
  geom_boxplot(outlier.shape = NA, coef = 1.0, fatten = 0) +
  labs(title = "",
       x = expression("Modeled " * italic(phnJ) * " (%)"),
       y = "Group") +
  stat_summary(fun = mean, geom = "point", shape = 23, size = 3, color = "red", fill=NA, stroke = 1.5) +
  scale_y_discrete(limits = c("Dep.MLD","phnJ"),
                   labels = c("Dep-MLD","phnJ")) +
  scale_x_continuous(limits = c(0, 30), breaks = c(0, 5, 10, 20)) +
  theme_bw() + 
  coord_trans(x = "log1p") +
  scale_fill_manual(values = c("phnJ" = "lightblue", "Dep.MLD" = "peachpuff")) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(size = rel(1.3)),
    axis.text.y = element_text(size = rel(1.3), face= "bold"),
    axis.title.x = element_text(size = 13),
    axis.title.y = element_text(size = 13)
  ) +
  guides(fill = FALSE)  

ggsave("phnJ_box.pdf", width = 5, height = 3, dpi = 300)


#######################################################################################################################################

# Figure 4 is composed of two panels that were created using Ocean Data View (V5.6.3) and combined in Illustrator. 

# You can find the raw data for these figures in Figshare at [https://doi.org/10.6084/m9.figshare.23912640.v1].

#######################################################################################################################################


#######################################################################################################################################

# Figure 5 was generated using Illustrator. 

#######################################################################################################################################

 