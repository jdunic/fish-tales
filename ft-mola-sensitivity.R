# Mola mola sensitivity

mola <- 
  read_csv('ft_mola_sensitivity.csv') %>% 
  select(-1) %>% 
  mutate(iucn_status = case_when(iucn.2016 %in% c('data.deficient', 'not.evaluated') ~ 'low', 
                                 iucn.2016 %in% c('least.concern', 'near.threatened') ~ 'med', 
                                 iucn.2016 %in% c('vulnerable', 'endangered', 'critically.endangered') ~ 'high')) %>% 
  mutate(iucn_status = factor(iucn_status, levels = c('low', 'med', 'high'))) %>% 
  mutate(prop.tl = reported.tl.cm / max.length) %>% 
  mutate(prop.tl = replace(prop.tl, prop.tl > 1, 1))

# Get sample sizes
# --------------------------
mola_trophic_sample_sizes_tl <- 
  mola %>% 
    filter(!is.na(trophic.category), !is.na(prop.tl)) %>% 
    group_by(trophic.category) %>% 
    summarise(sample_size = n()) %>%
    ungroup() %>% 
    get_sample_sizes()

# Mola: Trophic category - (prop.tl)
mola_trophic_qrs_tl <-
  mola %>% 
    filter(!is.na(trophic.category)) %>% 
    group_by(trophic.category) %>% 
    get_qr_tl(.)

mola_trophic_qdf_tl <- 
  # Needed to identify slope and intercept model values across columns
  left_join(filter(mola_trophic_qrs_tl, term == '(Intercept)'), filter(mola_trophic_qrs_tl, term == 'year'), by = c('trophic.category' = 'trophic.category', 'tau' = 'tau')) %>% 
    rename_at(vars(contains('.x')), funs(sub('.x', '_int', .))) %>% 
    rename_at(vars(contains('.y')), funs(sub('.y', '_slope', .))) %>% 
    left_join(., mola_trophic_sample_sizes_tl, by = c('trophic.category' = 'trophic.category', 'tau' = 'tau')) %>% 
    mk_cats()


orig_mega_qr_plot_tl <- 
  data %>%
    filter(trophic.category == 'megafauna') %>% 
    mk_qr_figs(df = ., x = 'year', y = 'prop.tl', quantile_df = filter(trophic_qdf_tl, trophic.category == 'megafauna', ratio_cat %in% c('conservative', 'minimum')), text_x = 1875, text_y = 0.2) + 
    theme(axis.title  = element_blank()) + 
    guides(colour = FALSE, linetype = FALSE) + 
    xlim(c(1869, 2020)) + 
    ylim(c(0.2, 1)) + 
    #geom_point(data = filter(data, genus.species == 'Molamola' & article.id %in% c(182, 183, 236, 287)), aes(x = year, y = prop.tl), colour = 'red') + 
    ggtitle(expression("(A) "~italic("Mola mola")))

mola_mega_qr_plot_tl <- 
    mola %>%
    filter(trophic.category == 'megafauna') %>% 
    mk_qr_figs(df = ., x = 'year', y = 'prop.tl', quantile_df = filter(mola_trophic_qdf_tl, trophic.category == 'megafauna', ratio_cat %in% c('conservative', 'minimum')), text_x = 1875, text_y = 0.2) + 
    theme(axis.title  = element_blank(), 
          axis.text.y = element_blank()) + 
    guides(colour = FALSE, linetype = FALSE) + 
    xlim(c(1869, 2020)) + 
    ylim(c(0.2, 1)) + 
    #geom_point(data = filter(mola, genus.species == 'Molamola' & article.id %in% c(182, 183, 236, 287)), aes(x = year, y = prop.tl), colour = 'red') + 
    ggtitle(expression("(B) "~italic("Mola tecta")))

#filter(mola, genus.species == 'Molamola' & article.id %in% c(182, 183, 236, 287)) %>% select(prop.tl)
#filter(data, genus.species == 'Molamola' & article.id %in% c(182, 183, 236, 287)) %>% select(prop.tl)

dev.new(width = 9, height = 3.8)

ylab <- grid::textGrob(label = 'Relative length', rot = 90, gp = gpar(cex = 1.4), vjust = 0.5)
xlab <- grid::textGrob(label = 'Year', hjust = 1.2, gp = gpar(cex = 1.4))
#
gridExtra::grid.arrange(
  ylab,
  cowplot::plot_grid(orig_mega_qr_plot_tl, mola_mega_qr_plot_tl, ncol = 2, rel_widths = c(1, 0.90)), 
  xlab, 
  layout_matrix = rbind(c(1, 2, NA), c(NA, 4, NA)), 
  widths = c(0.05, 1, 0.025), heights = c(1, 0.1)
  )

dev.copy2pdf(file = 'FigureS1-mola-sensitivity.pdf')

dev.copy2pdf(file = 'FigureS1-mola-sensitivity2.pdf')