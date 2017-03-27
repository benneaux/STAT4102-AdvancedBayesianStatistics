# # Agresti \& Coull
#source("Assignment1/Question1/1functions.R")
q1bchart <- ggplot(data = Q1bdata, aes(x = p, y = dens)) +
  geom_hline(yintercept = 1-a, linetype = 5, colour = "red") +
  geom_line() +
  scale_y_continuous(labels = scales::percent,limits = c(0.75,1)) +
  facet_wrap(~covinterval, nrow = 2) +
  theme_tufte(base_size = 14) + 
  scale_x_continuous(breaks = c(0,0.5,1)) +
  ylab("Coverage (%)") +
  xlab("Parameter (p)") +
  theme(panel.margin.x=unit(1.5, "lines")) +
  theme(axis.title.y=element_text(margin = margin(0,15,0,0))) +
  theme(axis.title.x = element_text(margin = margin(15,0,0,0)))