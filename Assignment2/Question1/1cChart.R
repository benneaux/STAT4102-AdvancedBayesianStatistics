# Agresti \& Coull 2

q1cchart <- ggplot(data = Q1cdata, aes(x = p, y = dens)) +
  geom_hline(yintercept = 1-a, linetype = 5, colour = "red") +
  geom_line(size = 0.25) +
  scale_y_continuous(labels = scales::percent, limits = c(0.8,1)) +
  facet_wrap(~covinterval, nrow = 2) +
  theme_tufte(base_size = 14) +
  scale_x_continuous(breaks = c(0,.5,1)) +
  ylab("Coverage (%)") +
  xlab("Parameter (p)") +
  theme(panel.margin.x=unit(1.5, "lines")) +
  theme(axis.title.y=element_text(margin = margin(0,15,0,0))) +
  theme(axis.title.x = element_text(margin = margin(15,0,0,0)))
