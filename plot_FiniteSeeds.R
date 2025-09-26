library(dplyr)
library(ggplot2)
library(stringr)
library(patchwork)
library(DescTools)

data = read.csv("~/Downloads/results_SexAllocation_FiniteSeeds_DD.csv")

## Population size
data %>% 
  group_by(E, K, FMAX, floor(GEN/500), DENS.DENP) %>% 
  summarize(POPSIZE = mean(POPSIZE),
            MEAN.Z = mean(MEAN.Z),
            VART.Z = mean(VART.Z),
            GEN = mean(GEN)) %>% 
  ggplot(aes(x = GEN, y = POPSIZE / ((1 - 1/(1+E))*500*100 + (1/(1+E))*500*K), color = factor(K))) +
  geom_line() +
  geom_hline(yintercept = 1.0, col = "gray", lty = "dashed") +
  facet_grid(rows = vars(E), cols = vars(FMAX, DENS.DENP)) +
  theme_classic() + labs(y = "Population (fraction of max.)")


## Trait
data %>% 
  group_by(E, K, FMAX, floor(GEN/1000), DENS.DENP) %>% 
  summarize(POPSIZE = mean(POPSIZE),
            MEAN.Z = mean(MEAN.Z),
            VART.Z = mean(VART.Z),
            GEN = mean(GEN)) %>% 
  ggplot(aes(x = GEN, y = MEAN.Z, color = factor(K))) +
  geom_line() +
  geom_hline(yintercept = 0.5, col = "gray", lty = "dashed") +
  facet_grid(rows = vars(E), cols = vars(FMAX, DENS.DENP)) +
  theme_classic() 


## Fraction of migrants
data %>% 
  filter(K<10) %>% 
  group_by(E, K, FMAX, floor(GEN/1000), DENS.DENP) %>% 
  summarize(MIG = mean(FRAC.MIG),
            EXP = 0.9,
            GEN = mean(GEN)) %>% 
  ggplot(aes(x = GEN, y = MIG/EXP, color = factor(K))) +
  geom_point(size = 0.2) +
  geom_hline(yintercept = 1.0, col = "gray", lty = "dashed") +
  facet_grid(rows = vars(E), cols = vars(FMAX, DENS.DENP)) +
  theme_classic() 


data %>% 
  filter(GEN > 5000, K<10, DENS.DENP == 1) %>% 
  group_by(E, K, FMAX, DENS.DENP) %>% 
  summarize(POPSIZE = mean(POPSIZE),
            MEAN.Z = mean(MEAN.Z),
            VART.Z = mean(VART.Z), 
            DZ = sqrt(VART.Z)) %>% 
  ggplot(aes(x = E, y = MEAN.Z, pch = factor(K), color = factor(FMAX)
            )) +
  #geom_errorbar(aes(ymin = MEAN.Z-DZ, ymax = MEAN.Z+DZ), width = 0.1, alpha = 0.3) +
  #geom_line() + 
  geom_point(alpha = 0.75, position = position_jitter(width = 0.02, height = 0)) +
  geom_line(data=filter(theo, SELF==0), aes(x = Nt, y=Value), 
            lwd = 1.5, alpha = 0.3, lty = 'solid', col = "gray") + 
  geom_hline(yintercept = 0.5, col = "gray", lty = "dashed") +
  #facet_grid(~FMAX) +
  theme_classic() + 
  #scale_shape_manual(values = c(16, 17, 15, 8, 0)) +
  scale_x_log10(breaks = c(1, 2, 5, 10, 50, 100)) +
  scale_color_viridis_d()
ggsave("~/SexAllocationEvolution/figs/fig2A.svg", height = 3, width = 6)

theoStab = expand.grid(FMAX = c(3, 10, 30, 300, 1000, 3000),
                       E = c(1, 2, 5, 10, 50)) %>% 
          mutate(e = 1/(1 + E),
                 minZ = 12*(1.0)/FMAX/(1-e) )

data %>% 
  filter(GEN > 40000, K<10, DENS.DENP == 1, E %in% c(2, 5, 10)) %>% 
  group_by(E, K, FMAX, DENS.DENP) %>% 
  summarize(POPSIZE = mean(POPSIZE),
            MEAN.Z = mean(MEAN.Z),
            VART.Z = mean(VART.Z), 
            DZ = sqrt(VART.Z)) %>% 
  ggplot(aes(x = 100*FMAX, y = MEAN.Z, color = factor(K),
  )) +
  geom_errorbar(aes(ymin = MEAN.Z-DZ, ymax = MEAN.Z+DZ), width = 0.2, alpha = 0.25) +
  geom_line(alpha = 0.7) + 
  geom_point(alpha = 0.7) +
  geom_point(data=filter(theo, SELF==0, E %in% c( 2, 5, 10)), aes(x = 5000, y=Value, col = factor(K)), 
             alpha = 0.7, pch = 15) + 
  geom_hline(yintercept = 0.5, col = "gray", lty = "dashed") +
  facet_grid(~E) +
  theme_classic() + 
  scale_x_log10(
    breaks = c(3, 10, 30, 100, 300, 1000, 3000, 5000),
    labels = c("3", "10", "30", "100", "300", "1k", "3k", "inf")) +
  scale_y_continuous(limits = c(0.47, 1.0), breaks = c(0.5, 0.6, 0.7, 0.8, 0.9, 1.0)) +
  scale_color_grey()
ggsave("~/SexAllocationEvolution/figs/fig2A_v2.svg", height = 3, width = 7)


