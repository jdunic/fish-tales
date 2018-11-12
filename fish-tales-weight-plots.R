# ------------------------------------------------------------------------------
# Supplemental plots (prop.weight)
# ------------------------------------------------------------------------------
#
# Get sample sizes
# --------------------------
trophic_sample_sizes_wt <- 
  data %>% 
    filter(!is.na(trophic.category), !is.na(prop.weight)) %>% 
    group_by(trophic.category) %>% 
    summarise(sample_size = n()) %>%
    ungroup() %>% 
    get_sample_sizes()
#
iucn_sample_sizes_wt <- 
  data %>% 
    filter(!is.na(iucn_status), !is.na(prop.weight)) %>% 
    group_by(iucn_status) %>% 
    summarise(sample_size = n()) %>%
    ungroup() %>% 
    get_sample_sizes()


# Run quantile regressions
# --------------------------
# Trophic category - (prop.weight)
trophic_qrs_wt <-
  data %>% 
    filter(!is.na(trophic.category)) %>% 
    group_by(trophic.category) %>% 
    get_qr_wt(.)
#
# IUCN category - (prop.weight)
iucn_qrs_wt <-
  data %>% 
    filter(!is.na(iucn_status)) %>% 
    group_by(iucn_status) %>% 
    get_qr_wt(.)

# Make master dataframes for plotting
# ------------------------------------------------------------------------------
trophic_qdf_wt <- 
  # Needed to identify slope and intercept model values across columns
  left_join(filter(trophic_qrs_wt, term == '(Intercept)'), filter(trophic_qrs_wt, term == 'year'), by = c('trophic.category' = 'trophic.category', 'tau' = 'tau')) %>% 
    rename_at(vars(contains('.x')), funs(sub('.x', '_int', .))) %>% 
    rename_at(vars(contains('.y')), funs(sub('.y', '_slope', .))) %>% 
    left_join(., trophic_sample_sizes_wt, by = c('trophic.category' = 'trophic.category', 'tau' = 'tau')) %>% 
    mk_cats()
#
iucn_qdf_wt <- 
  left_join(filter(iucn_qrs_wt, term == '(Intercept)'), filter(iucn_qrs_wt, term == 'year'), by = c('iucn_status' = 'iucn_status', 'tau' = 'tau')) %>% 
    rename_at(vars(contains('.x')), funs(sub('.x', '_int', .))) %>% 
    rename_at(vars(contains('.y')), funs(sub('.y', '_slope', .))) %>% 
    left_join(., iucn_sample_sizes_wt) %>% 
    mk_cats()


get_summ(trophic_qdf_wt, 'trophic.category', 'ratio_cat') %>%
  filter(trophic.category == 'megafauna')


get_summ(iucn_qdf_tl, 'iucn_status', 'ratio_cat') %>%
  filter(iucn_status == 'low')


# ------------------------------------------------------------------------------
# Trophic groups relative weight plots
# ------------------------------------------------------------------------------
dev.new(width = 10.5, height = 3.8)

game_qr_plot_wt <- 
  data %>%
    filter(trophic.category == 'gamefish') %>% 
    mk_qr_figs(df = ., x = 'year', y = 'prop.weight', quantile_df = filter(trophic_qdf_wt, trophic.category == 'gamefish', ratio_cat %in% c('conservative', 'minimum')), text_x = 2000, text_y = 0.02) + 
    ggtitle("(A) Pelagic gamefish") + 
    theme(axis.title = element_blank()) + 
    guides(colour = FALSE, linetype = FALSE) + 
    xlim(c(1869, 2020))
#
oshark_qr_plot_wt <- 
  data %>%
    filter(trophic.category == 'oceanic.shark') %>% 
    mk_qr_figs(df = ., x = 'year', y = 'prop.weight', quantile_df = filter(trophic_qdf_wt, trophic.category == 'oceanic.shark', ratio_cat %in% c('conservative', 'minimum')), text_x = 1874, text_y = 0.02) + 
    ggtitle("(B) Oceanic sharks") + 
    theme(axis.title  = element_blank(), 
          axis.text.y = element_blank()) + 
    xlim(c(1869, 2020)) + 
    guides(colour = FALSE, linetype = FALSE)
mega_qr_plot_wt <- 
    data %>%
    filter(trophic.category == 'megafauna') %>% 
    mk_qr_figs(df = ., x = 'year', y = 'prop.weight', quantile_df = filter(trophic_qdf_wt, trophic.category == 'megafauna', ratio_cat %in% c('conservative', 'minimum')), text_x = 1874, text_y = 0.02) + 
    ggtitle("(C) Charismatic megafish") + 
    theme(axis.title  = element_blank(), 
          axis.text.y = element_blank()) + 
    guides(colour = FALSE, linetype = FALSE) + 
    xlim(c(1869, 2020))
#
ylab <- grid::textGrob(label = 'Relative weight', rot = 90, gp = gpar(cex = 1.2), vjust = 0.5)
xlab <- grid::textGrob(label = 'Year', hjust = 0, gp = gpar(cex = 1.2))

# Where the plotting magic happens
gridExtra::grid.arrange(
  ylab,
  cowplot::plot_grid(game_qr_plot_wt, oshark_qr_plot_wt, mega_qr_plot_wt, ncol = 3, rel_widths = c(1, 0.93, 0.93)), 
  xlab, 
  layout_matrix = rbind(c(1, 2, NA), c(NA, 4, NA)), 
  widths = c(0.05, 1, 0.025), heights = c(1, 0.1)
  )

dev.copy2pdf(file = 'FigureS2-trophic-weight.pdf')
dev.off()

# IUCN groups relative weight plots
# ------------------------------------------------------------------------------
dev.new(width = 10.5, height = 3.8)

low_qr_plot_wt <- 
  data %>%
    filter(iucn_status == 'low') %>% 
    mk_qr_figs(df = ., x = 'year', y = 'prop.weight', quantile_df = filter(iucn_qdf_wt, iucn_status == 'low', ratio_cat %in% c('conservative', 'minimum')), text_x = 2000, text_y = 0.02) + 
    ggtitle("(A) Unknown") + 
    theme(axis.title = element_blank()) + 
    guides(colour = FALSE, linetype = FALSE) + 
    xlim(c(1869, 2020))
#
med_qr_plot_wt <- 
    data %>%
    filter(iucn_status == 'med') %>% 
    mk_qr_figs(df = ., x = 'year', y = 'prop.weight', quantile_df = filter(iucn_qdf_wt, iucn_status == 'med', ratio_cat %in% c('conservative', 'minimum')), text_x = 1876, text_y = 0.02) + 
    ggtitle("(B) Low risk") + 
    theme(axis.title  = element_blank(), 
          axis.text.y = element_blank()) + 
    guides(colour = FALSE, linetype = FALSE) + 
    xlim(c(1869, 2020))
#
high_qr_plot_wt <- 
  data %>%
    filter(iucn_status == 'high') %>% 
  mk_qr_figs(df = ., x = 'year', y = 'prop.weight', quantile_df = filter(iucn_qdf_wt, iucn_status == 'high', ratio_cat %in% c('conservative', 'minimum')), text_x = 1878, text_y = 0.02) + 
    ggtitle("(C) High risk") + 
    theme(axis.title  = element_blank(), 
          axis.text.y = element_blank()) + 
    guides(colour = FALSE, linetype = FALSE) + 
    xlim(c(1869, 2020))

#
ylab <- grid::textGrob(label = 'Relative weight', rot = 90, gp = gpar(cex = 1.2), vjust = 0.5)
xlab <- grid::textGrob(label = 'Year', hjust = 0, gp = gpar(cex = 1.2))
#
# Where the plotting magic happens
gridExtra::grid.arrange(
  ylab,
  cowplot::plot_grid(low_qr_plot_wt, med_qr_plot_wt, high_qr_plot_wt, ncol = 3, rel_widths = c(1, 0.93, 0.93)), 
  xlab, 
  layout_matrix = rbind(c(1, 2, NA), c(NA, 4, NA)), 
  widths = c(0.05, 1, 0.025), heights = c(1, 0.1)
  )

dev.copy2pdf(file = 'FigureS3-iucn-weight.pdf')
dev.off()