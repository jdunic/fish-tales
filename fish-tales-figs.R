### Fish tales quantile regression
library(tidyverse)
library(quantreg)
library(cowplot)
library(gridExtra)
library(kableExtra)
library(grid)

set.seed(1)

data <- 
  read_csv('ft_Dunic.csv') %>% 
  select(-1) %>% 
  mutate(iucn_status = case_when(iucn.2016 %in% c('data.deficient', 'not.evaluated') ~ 'low', 
                                 iucn.2016 %in% c('least.concern', 'near.threatened') ~ 'med', 
                                 iucn.2016 %in% c('vulnerable', 'endangered', 'critically.endangered') ~ 'high')) %>% 
  mutate(iucn_status = factor(iucn_status, levels = c('low', 'med', 'high')))

# Get sample sizes (i.e., low, min, conservative quantile sample sizes)
get_sample_sizes <- function(df) {
  df %>% 
    mutate(`10` = sample_size * 0.10, 
           `15` = sample_size * 0.15, 
           `20` = sample_size * 0.20, 
           `25` = sample_size * 0.25, 
           `50` = sample_size * 0.50, 
           `75` = sample_size * (1-0.75),
           `80` = sample_size * (1-0.80), 
           `85` = sample_size * (1-0.85), 
           `90` = sample_size * (1-0.90) 
           ) %>% 
    mutate(n_text = paste0('n = ', sample_size)) %>% 
    gather(key = tau, value = ratio, `10`:`90`) %>% 
    mutate(tau = as.numeric(tau) / 100) %>% 
    mutate(ratio_cat = case_when(ratio >= 10 ~ "conservative", 
                                 (ratio < 10 & ratio >= 5) ~ "minimum", 
                                 ratio < 5 ~ "low")) %>% 
    mutate(tau_cat = case_when(tau == 0.10 | tau == 0.90 ~ '10th', 
                           tau == 0.15 | tau == 0.85 ~ '15th', 
                           tau == 0.20 | tau == 0.80 ~ '20th', 
                           tau == 0.25 | tau == 0.75 ~ '25th', 
                           tau == 0.50 ~ '50th')) %>% 
    mutate(tau_cat = ifelse(tau == 0.5, '50th', 'other'))
}

get_qr_tl <- function(df) {
  df %>% 
  do(fit = rq(prop.tl ~ year, tau = c(0.10, 0.15, 0.20, 0.25, 0.5, 0.90, 0.85, 0.80, 0.75), data = ., na.action = na.omit) %>% broom::tidy(., se.type = "boot")) %>% 
  unnest(fit) %>% 
  ungroup() %>% 
    mutate(sig = ifelse(p.value <= 0.05, 'Sig', 'Not sig'))  
}

get_qr_wt <- function(df) {
  df %>% 
    do(fit = rq(prop.weight ~ year, tau = c(0.10, 0.15, 0.20, 0.25, 0.5, 0.90, 0.85, 0.80, 0.75), data = ., na.action = na.omit) %>% broom::tidy(., se.type = "boot")) %>% 
    unnest(fit) %>% 
    ungroup() %>% 
    mutate(sig = ifelse(p.value <= 0.05, 'Sig', 'Not sig'))
}

# Add categories for plotting colours and linetypes
mk_cats <- function(df) {
  df %>% 
    mutate(ratio_cat = factor(ratio_cat, levels = c("conservative", "minimum", "low"))) %>% 
    mutate(sig_slope = factor(sig_slope, levels = c('Sig', 'Not sig'))) %>% 
    mutate(tau_cat   = factor(tau_cat, levels = c('50th', 'other'))) %>%
    mutate(linetype  = case_when(sig_slope == 'Sig' ~ 'solid', 
                                 sig_slope != 'Sig' & tau_cat == '50th' ~ 'longdashed', 
                                 sig_slope != 'Sig' & tau_cat == 'other' ~ 'dotted')) %>% 
    mutate(line_col  = case_when(ratio_cat %in% c('conservative', 'minimum') & sig_slope == 'Sig' ~ 'blue', 
                                 sig_slope == 'Not sig' ~ 'grey')) %>%
    mutate(linetype  = factor(linetype, levels = c('solid', 'longdashed', 'dotted')), 
           line_col  = factor(line_col, levels = c('blue', 'grey')))
}

# Plot generating function
mk_qr_figs <- function(df, x, y, quantile_df, text_x = NULL, text_y = NULL, n_text = TRUE) {
  # Setting the positions of the 'n = ' text on the plot. By default calculates 
  # the min values in the data.
  if(is.null(text_x) | !is.numeric(text_x)) {
    text_x = min(df[x], na.rm = TRUE) 
    message(paste0("text_x = ", text_x))
  }
  if(is.null(text_y) | !is.numeric(text_y)) {
    text_y = min(df[y], na.rm = TRUE) 
    message(paste0("text_y = ", text_y))
  }
  qr_plot <- 
    ggplot(data = df, aes_string(x = x, y = y)) + 
      geom_point(size = 2) + 
      theme_classic() + 
      theme(axis.line = element_line(), 
            legend.title = element_blank(), 
            axis.text = element_text(size = 12)) + 
      geom_abline(data = quantile_df, aes(slope = estimate_slope, intercept = estimate_int, colour = line_col, linetype = linetype)) +
      # These can be changed. E.g., http://www.sthda.com/english/wiki/line-types-in-r-lty
      scale_linetype_manual(values = c(1, 2, 3), drop = FALSE) + 
      scale_colour_manual(values = c("blue", "grey"), drop = FALSE)
  if(n_text == TRUE) {
    qr_plot <- 
      qr_plot + 
      geom_text(data = distinct(quantile_df, n_text), aes(x = text_x + 5, y = text_y, label = n_text), size = 5)
  }
  if(n_text == FALSE) {
    qr_plot <- qr_plot
  }
  return(qr_plot)
}

get_summ <- function(df, ...) {
  df %>% 
    select_('tau', 'estimate_slope', 'p.value_slope', 'conf.low_slope', 'conf.high_slope', 'ratio', ...) %>% 
    arrange(tau)
}


# ----------------------------------------------------
# Figure 1
# ----------------------------------------------------
# Get sample sizes
# --------------------------
all_sample_sizes_tl <- 
  data %>% 
    filter(!is.na(prop.tl)) %>% 
    summarise(sample_size = n()) %>%
    ungroup() %>% 
    get_sample_sizes()
#
all_sample_sizes_wt <- 
  data %>% 
    filter(!is.na(prop.weight)) %>% 
    summarise(sample_size = n()) %>%
    ungroup() %>% 
    get_sample_sizes()

# Run quantile regressions
# --------------------------
# All data - (prop.tl)
all_qrs_tl <-
  data %>% 
    get_qr_tl(.)
#
# All data - (prop.weight)
all_qrs_wt <-
  data %>% 
    get_qr_wt(.)

# Make master dataframes for plotting
# ------------------------------------------------------------------------------
all_qdf_tl <- 
  # Needed to identify slope and intercept model values across columns
  left_join(filter(all_qrs_tl, term == '(Intercept)'), filter(all_qrs_tl, term == 'year'), by = c('tau' = 'tau')) %>% 
    rename_at(vars(contains('.x')), funs(sub('.x', '_int', .))) %>% 
    rename_at(vars(contains('.y')), funs(sub('.y', '_slope', .))) %>% 
    left_join(., all_sample_sizes_tl, by = c('tau' = 'tau')) %>% 
    mk_cats()
#
all_qdf_wt <- 
  # Needed to identify slope and intercept model values across columns
  left_join(filter(all_qrs_wt, term == '(Intercept)'), filter(all_qrs_wt, term == 'year'), by = c('tau' = 'tau')) %>% 
    rename_at(vars(contains('.x')), funs(sub('.x', '_int', .))) %>% 
    rename_at(vars(contains('.y')), funs(sub('.y', '_slope', .))) %>% 
    left_join(., all_sample_sizes_wt, by = c('tau' = 'tau')) %>% 
    mk_cats()

# Plot Figure 1
dev.new(width = 9, height = 3.8)

A_length <- 
  mk_qr_figs(df = data, x = 'year', y = 'prop.tl', quantile_df = all_qdf_tl, n_text = TRUE, text_x = 1876, text_y = 0.02) +
    theme(axis.title = element_blank()) +
    ylim(c(0, 1.0)) + 
    xlim(c(1869, 2017)) + 
    guides(colour = FALSE, linetype = FALSE) + 
    ggtitle('(A) Length')
#
B_weight <- 
  mk_qr_figs(df = data, x = 'year', y = 'prop.weight', quantile_df = all_qdf_wt, n_text = TRUE, text_x = 1876, text_y = 0.02) +
    #annotate(geom = 'text', x = 1868, y = 0.02, label = '(B) Weight', size = 5, colour = 'grey30', hjust = 0) + 
    theme(axis.title = element_blank(), 
          axis.text.y = element_blank(), 
          axis.ticks.y = element_blank()) +
    ylim(c(0, 1.0)) + 
    xlim(c(1869, 2017)) + 
    guides(colour = FALSE, linetype = FALSE) + 
    ggtitle('(B) Weight')
#
ylab <- grid::textGrob(label = 'Relative size', rot = 90, gp = gpar(cex = 1.4), vjust = 0.5)
xlab <- grid::textGrob(label = 'Year', hjust = 1.2, gp = gpar(cex = 1.4))
#
gridExtra::grid.arrange(
  ylab,
  cowplot::plot_grid(A_length, B_weight, ncol = 2, rel_widths = c(1, 0.90)), 
  xlab, 
  layout_matrix = rbind(c(1, 2, NA), c(NA, 4, NA)), 
  widths = c(0.05, 1, 0.025), heights = c(1, 0.1)
  )

dev.copy2pdf(file = 'Figure1.pdf')
dev.off()

# --------------------------
# Figures 2 & 3 - Trophic/IUCN prop.tl ~ year - sig and sample size
# -------------------------- 
#
# Get sample sizes
# --------------------------
trophic_sample_sizes_tl <- 
  data %>% 
    filter(!is.na(trophic.category), !is.na(prop.tl)) %>% 
    group_by(trophic.category) %>% 
    summarise(sample_size = n()) %>%
    ungroup() %>% 
    get_sample_sizes()
#
iucn_sample_sizes_tl <- 
  data %>% 
    filter(!is.na(iucn_status), !is.na(prop.tl)) %>% 
    group_by(iucn_status) %>% 
    summarise(sample_size = n()) %>%
    ungroup() %>% 
    get_sample_sizes()


# Run quantile regressions
# --------------------------
# Trophic category - (prop.tl)
trophic_qrs_tl <-
  data %>% 
    filter(!is.na(trophic.category)) %>% 
    group_by(trophic.category) %>% 
    get_qr_tl(.)
# IUCN category - (prop.tl)
iucn_qrs_tl <-
  data %>% 
    filter(!is.na(iucn_status)) %>% 
    group_by(iucn_status) %>% 
    get_qr_tl(.)

# Make master dataframes for plotting
# ------------------------------------------------------------------------------
trophic_qdf_tl <- 
  # Needed to identify slope and intercept model values across columns
  left_join(filter(trophic_qrs_tl, term == '(Intercept)'), filter(trophic_qrs_tl, term == 'year'), by = c('trophic.category' = 'trophic.category', 'tau' = 'tau')) %>% 
    rename_at(vars(contains('.x')), funs(sub('.x', '_int', .))) %>% 
    rename_at(vars(contains('.y')), funs(sub('.y', '_slope', .))) %>% 
    left_join(., trophic_sample_sizes_tl, by = c('trophic.category' = 'trophic.category', 'tau' = 'tau')) %>% 
    mk_cats()
#
iucn_qdf_tl <- 
  left_join(filter(iucn_qrs_tl, term == '(Intercept)'), filter(iucn_qrs_tl, term == 'year'), by = c('iucn_status' = 'iucn_status', 'tau' = 'tau')) %>% 
    rename_at(vars(contains('.x')), funs(sub('.x', '_int', .))) %>% 
    rename_at(vars(contains('.y')), funs(sub('.y', '_slope', .))) %>% 
    left_join(., iucn_sample_sizes_tl) %>% 
    mk_cats()

get_summ(iucn_qdf_tl, 'iucn_status', 'ratio_cat') %>%
  filter(iucn_status == 'high')

# Trophic groups relative length plots
# ------------------------------------------------------------------------------
dev.new(width = 10.5, height = 3.8)

game_qr_plot_tl <- 
  data %>%
    filter(trophic.category == 'gamefish') %>% 
    mk_qr_figs(df = ., x = 'year', y = 'prop.tl', quantile_df = filter(trophic_qdf_tl, trophic.category == 'gamefish', ratio_cat %in% c('conservative', 'minimum')), text_x = 1878, text_y = 0.2) + 
    ggtitle("(A) Pelagic gamefish") + 
    theme(axis.title = element_blank()) + 
    guides(colour = FALSE, linetype = FALSE) + 
    xlim(c(1869, 2020)) + 
    ylim(c(0.2, 1))
#
oshark_qr_plot_tl <- 
  data %>%
    filter(trophic.category == 'oceanic.shark') %>% 
    mk_qr_figs(df = ., x = 'year', y = 'prop.tl', quantile_df = filter(trophic_qdf_tl, trophic.category == 'oceanic.shark', ratio_cat %in% c('conservative', 'minimum')), text_x = 1878, text_y = 0.2) + 
    ggtitle("(B) Oceanic sharks") + 
    theme(axis.title  = element_blank(), 
          axis.text.y = element_blank()) + 
    xlim(c(1869, 2020)) + 
    guides(colour = FALSE, linetype = FALSE) + 
    ylim(c(0.2, 1))
#
mega_qr_plot_tl <- 
    data %>%
    filter(trophic.category == 'megafauna') %>% 
    mk_qr_figs(df = ., x = 'year', y = 'prop.tl', quantile_df = filter(trophic_qdf_tl, trophic.category == 'megafauna', ratio_cat %in% c('conservative', 'minimum')), text_x = 1878, text_y = 0.2) + 
    ggtitle("(C) Charismatic megafish") + 
    theme(axis.title  = element_blank(), 
          axis.text.y = element_blank()) + 
    guides(colour = FALSE, linetype = FALSE) + 
    xlim(c(1869, 2020)) + 
    ylim(c(0.2, 1))


#
ylab <- grid::textGrob(label = 'Relative length', rot = 90, gp = gpar(cex = 1.2), vjust = 0.5)
xlab <- grid::textGrob(label = 'Year', hjust = 0, gp = gpar(cex = 1.2))
#
gridExtra::grid.arrange(
  ylab,
  cowplot::plot_grid(game_qr_plot_tl, oshark_qr_plot_tl, mega_qr_plot_tl, ncol = 3, rel_widths = c(1, 0.93, 0.93)), 
  xlab, 
  layout_matrix = rbind(c(1, 2, NA), c(NA, 4, NA)), 
  widths = c(0.05, 1, 0.025), heights = c(1, 0.1)
  )

dev.copy2pdf(file = 'Figure2.pdf')
dev.off()

# Show the overlap in two of the quantile regression fits
ggplot(data = filter(data, trophic.category == 'oceanic.shark'), aes(x = year, y = prop.tl)) + 
  geom_point() + 
  geom_abline(data = filter(trophic_qdf_tl, trophic.category == 'oceanic.shark', ratio_cat %in% c('conservative', 'minimum')), aes(slope = estimate_slope, intercept = estimate_int, colour = tau))

# IUCN groups relative length plots
# ------------------------------------------------------------------------------
dev.new(width = 10.5, height = 3.8)

low_qr_plot_tl <- 
  data %>%
    filter(iucn_status == 'low') %>% 
    mk_qr_figs(df = ., x = 'year', y = 'prop.tl', quantile_df = filter(iucn_qdf_tl, iucn_status == 'low', ratio_cat %in% c('conservative', 'minimum')), text_x = 1878, text_y = 0.2) + 
    ggtitle("(A) Unknown") + 
    theme(axis.title = element_blank()) + 
    guides(colour = FALSE, linetype = FALSE) + 
    xlim(c(1869, 2020))
#
med_qr_plot_tl <- 
    data %>%
    filter(iucn_status == 'med') %>% 
    mk_qr_figs(df = ., x = 'year', y = 'prop.tl', quantile_df = filter(iucn_qdf_tl, iucn_status == 'med', ratio_cat %in% c('conservative', 'minimum')), text_x = 1878, text_y = 0.2) + 
    ggtitle("(B) Low risk") + 
    theme(axis.title  = element_blank(), 
          axis.text.y = element_blank()) + 
    guides(colour = FALSE, linetype = FALSE) + 
    xlim(c(1869, 2020))
#
high_qr_plot_tl <- 
    data %>% 
    filter(iucn_status == 'high') %>% 
    mk_qr_figs(df = ., x = 'year', y = 'prop.tl', quantile_df = filter(iucn_qdf_tl, iucn_status == 'high', ratio_cat %in% c('conservative', 'minimum')), text_x = 1878, text_y = 0.2) + 
    ggtitle("(C) High risk") + 
    theme(axis.title  = element_blank(), 
          axis.text.y = element_blank()) + 
    xlim(c(1869, 2020)) + 
    guides(colour = FALSE, linetype = FALSE)

#
ylab <- grid::textGrob(label = 'Relative length', rot = 90, gp = gpar(cex = 1.2), vjust = 0.5)
xlab <- grid::textGrob(label = 'Year', hjust = 0, gp = gpar(cex = 1.2))
#
gridExtra::grid.arrange(
  ylab,
  cowplot::plot_grid(low_qr_plot_tl, med_qr_plot_tl, high_qr_plot_tl, ncol = 3, rel_widths = c(1, 0.93, 0.93)), 
  xlab, 
  layout_matrix = rbind(c(1, 2, NA), c(NA, 4, NA)), 
  widths = c(0.05, 1, 0.025), heights = c(1, 0.1)
  )

dev.copy2pdf(file = 'Figure3.pdf')
dev.off()

