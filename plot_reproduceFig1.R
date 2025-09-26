library(dplyr)
library(ggplot2)
library(stringr)
library(patchwork)
library(DescTools)
library(ggfortify)  # for autoplot
library(FactoMineR) # for more advanced PCA
library(factoextra)
library(ggrepel)
library(pls) # for PLS


data = read.csv("~/Downloads/results_SexAllocation_InfiniteSeeds.csv")
theo = read.csv("~/SexAllocationEvolution/Nt_K_Value_Table.csv")

data %>% 
  filter(GEN > 10000, SELF == 0.0, K<10) %>% 
  group_by(E, K) %>% 
  summarise(Z = mean(MEAN.Z), 
            dZ = sqrt(mean(VARB.Z+VARW.Z))) %>% 
  ggplot(aes(x = E, y = Z, col=factor(K))) +
  geom_hline(yintercept = 0.5, col = "gray", lty="dashed") +
  geom_errorbar(aes(ymin = Z-dZ, ymax = Z+dZ), width = 0.05, alpha = 0.5) +
  geom_point(alpha = 0.75) +
  geom_line(data=filter(theo, SELF==0), aes(x = Nt, y=Value)) + 
  theme_classic() +
  scale_x_log10(breaks = c(1, 2, 5, 10, 50, 100)) +
  labs(x = "Average deme lifetime",
       y = "Stationary sex allocation",
       col = "# colonizers") + ylim(c(0.4, 1.0)) +
  scale_color_grey()
ggsave("~/SexAllocationEvolution/figs/fig1A_v2_updatedSimulations.svg", height = 3, width = 4)


data %>% 
  filter(GEN > 40000, K == 2) %>% 
  group_by(E, SELF) %>% 
  summarise(Z = mean(MEAN.Z), 
            dZ = sqrt(mean(VARB.Z+VARW.Z))) %>% 
  ggplot(aes(x = E, y = Z, col=factor(SELF))) +
  geom_hline(yintercept = 0.5, col = "gray", lty="dashed") +
  geom_errorbar(aes(ymin = Z-dZ, ymax = Z+dZ), width = 0.05, alpha = 0.5) +
  geom_point(alpha = 0.75) + 
  geom_line(data=filter(theo, K == 2), aes(x = Nt, y=Value)) + 
  theme_classic() +
  scale_x_log10(breaks = c(1, 2, 5, 10, 50, 100)) +
  labs(x = "Average deme lifetime",
       y = "Stationary sex allocation",
       col = "Selfing") + ylim(c(0.4, 1.0)) +
  scale_color_grey()
ggsave("~/SexAllocationEvolution/figs/fig1B.svg", height = 3, width = 4)


data %>% 
  filter(K<10) %>% 
  group_by(E, K, floor(GEN/5000) ) %>% 
  summarize(MEAN.Z = mean(MEAN.Z),
            VARB.Z = mean(VARB.Z),
            VARW.Z = mean(VARW.Z),
            GEN = mean(GEN)) %>% 
  ggplot(aes(x = GEN, y = MEAN.Z)) +
  geom_ribbon(aes(ymin = MEAN.Z-sqrt(VARB.Z+VARW.Z), 
                    ymax = MEAN.Z+sqrt(VARB.Z+VARW.Z)), 
                alpha = 0.5, fill = "red") +
  geom_ribbon(aes(ymin = MEAN.Z-sqrt(VARB.Z), 
                    ymax = MEAN.Z+sqrt(VARB.Z)), 
                alpha = 0.5, fill = "darkgreen") +
  geom_line() +
  geom_hline(yintercept = 0.5, col = "gray", lty = "dashed") +
  facet_grid(rows = vars(E), cols = vars(K)) +
  theme_classic() + scale_x_log10() 

data %>% 
  filter(GEN > 100000) %>% 
  group_by(E, K) %>% 
  summarize(across(everything(), mean, na.rm = TRUE)) %>% View


###########################
### Linear regressions ####
###########################

dataMP = read.csv("~/Downloads/results_SexAllocation_Rnd_Migrant.csv")
dataP = read.csv("~/Downloads/results_SexAllocation_Rnd_Propagule.csv")
data = read.csv("~/Downloads/results_SexAllocation_Rnd_FiniteSeeds.csv")
dataTH = read.csv("~/Downloads/results_SexAllocation_Rnd_Theory.csv", 
                  header = F, col.names = c("FIS", "FIT", "FST", "MEAN.Z", "MODEL") )

# Prepare data 
data_long = data %>%
  pivot_longer(cols = c(FST, FIS, FIT),
               names_to = "Metric", values_to = "Value") %>%
  drop_na(MEAN.Z, Value)

# Compute R² per Metric
r2_labels = data_long %>%
  group_by(Metric, MODEL, DISP) %>%
  summarise(
    r2 = summary(lm(Value ~ MEAN.Z))$r.squared,
    .groups = "drop"
  ) %>%
  mutate(
    label = paste0("R² = ", round(r2, 2)),
    x = 0.05,
    y = 0.95
  )

# Plot with regression lines and R²
data_long %>% 
  filter(Metric %in% c("FST", "FIS", "FIT") ) %>% 
ggplot(aes(x = MEAN.Z, y = Value)) +
  geom_point(color = "gray", alpha = 0.4, size = 0.5) +
  geom_smooth(method = "lm", se = FALSE, color = "black", size = 0.5) +
  geom_text(data = filter(r2_labels, Metric %in% c("FST", "FIS", "FIT") ), aes(x = x, y = y, label = label),
            inherit.aes = FALSE, hjust = 0, vjust = 1, size = 3) +
  geom_vline(xintercept = 0.5, col = "gray", lty = "dashed") +
  geom_hline(yintercept = 0.0, col = "gray", lty = "dashed") +
  facet_grid(Metric~MODEL) +
  scale_color_viridis_d() +
  theme_classic() +
  theme(legend.position = "none") +
  lims(x = c(0, 1), y = c(-0.02, 1))


#######################################################
### LOESS: locally estimated scatterplot smoothing ####
#######################################################

data_long %>%
  #filter(Metric %in% c("FST", "FIS", "FIT")) %>%
  ggplot(aes(x = MEAN.Z, y = Value, col = factor(TYPE), pch = factor(TYPE))) +
  geom_point(alpha = 0.5, size = 0.55) +       # Observed points
  geom_smooth(method = "loess", se = TRUE, size = 0.75, span = 0.5) +  # LOESS curve
  facet_grid(Metric ~ MODEL) +
  #theme_bw() +
  theme_classic() + #theme(legend.position = "none") +
  scale_color_manual(values = c("black", "darkgray")) +
  scale_shape_manual(values = c(16, 15)) +
  lims(x = c(0.5, 1), y = c(-0.02, 1)) +
  labs(x = "MEAN.Z", y = "Value")
ggsave("~/SexAllocationEvolution/figs/fig3_withTHEORY_v2.svg", height = 6, width = 6)


############
### PCA ####
############

# Prepare data
pca_data <- data %>%
  select(FST, FIS, FIT, MEAN.Z) %>%
  drop_na()

# Run PCA including MEAN.Z
pca_res <- prcomp(pca_data, scale. = TRUE)

# Extract scores (individual positions)
scores <- as.data.frame(pca_res$x)
scores$MEAN.Z_color <- pca_data$MEAN.Z  # keep color variable

# Extract loadings (variable contributions)
loadings <- as.data.frame(pca_res$rotation)
loadings$var <- rownames(loadings)

# Scale loadings for plotting
scale_factor <- max(abs(scores$PC1))  # crude scaling
loadings_scaled <- loadings %>%
  mutate(across(PC1:PC2, ~ . * scale_factor))

# Plot
ggplot(scores, aes(x = PC1, y = PC2, color = MEAN.Z_color)) +
  geom_point(alpha = 0.6, size = 0.5) +
  geom_segment(data = filter(loadings_scaled, var != "MEAN.Z"),
             aes(x = 0, y = 0, xend = PC1, yend = PC2),
             arrow = arrow(length = unit(0.2, "cm")),
             color = "black", size = 0.7) +
  geom_segment(data = filter(loadings_scaled, var == "MEAN.Z"),
               aes(x = 0, y = 0, xend = PC1, yend = PC2),
               arrow = arrow(length = unit(0.5, "cm")),
               color = "red", size = 1.2) +
  geom_text_repel(data = loadings_scaled,
                  aes(x = PC1, y = PC2, label = var),
                  color = "black", size = 3) +
  annotate("point", x = 0, y = 0, color = "red", size = 2) +
  scale_color_viridis_c() +
  theme_classic() +
  labs(title = "Principal Component Analysis")

#####################
##### Histograms ####
#####################
data %>%
  filter(K>1) %>% 
  ggplot(aes(x = MEAN.Z, fill = factor(DISP))) +
  geom_density(aes(y = after_stat(density)), position = "identity", 
               alpha = 0.6, col = "white") +
  facet_grid(cols = vars(MODEL)) +
  scale_color_viridis_c() +
  theme_classic()
ggsave("~/SexAllocationEvolution/figs/fig_PropVSMigr.svg", height = 3, width = 6)

data %>%
  filter(MODEL == "randomMating") %>% 
  ggplot(aes(col = MEAN.Z, pch = factor(MEAN.Z<0.5), x = MP-2*MS, y = e)) +
  geom_point(alpha = 0.6, size = 2.0) +
  geom_segment(x = 0, y = 0, xend = 0.1, yend = 0.1, col = "black", lty = "dashed") +
  facet_grid(cols = vars(DISP)) +
  #scale_size_identity() +
  scale_color_gradient2(
    midpoint = 0.5,      # value where colours diverge
    low = "blue",         # colour for low end
    mid = "gray",        # colour for midpoint
    high = "red",         # colour for high end
    limits = c(0, 1)      # min and max of the scale
  ) +  theme_classic() + 
  coord_cartesian(ylim = c(0.0, 0.5), xlim = c(-0.2, 0.1), clip = "on", expand = F)

ggsave("~/SexAllocationEvolution/figs/fig_Male-biased.svg", height = 3, width = 6)

