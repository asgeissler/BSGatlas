source('analysis/00_load.R')
library(venn)

load('analysis/05_bsg_boundaries.rda')

load('data/01_bsub_raw.rda')
load('data/01_nicolas.rda')
load('data/03_dbtbs.rda')

###############################################################################
# count of annotation origins

bsg.boundaries %>%
  map2(c('TSS', 'TTS'), ~ mutate(.x, type = .y)) %>%
  map(select, id, type, src) %>%
  bind_rows -> dat

dat %>%
  separate_rows(src, sep = ';') %>%
  mutate(
    # src = fct_recode(src,
    #  "Nicolas et al.\nDownshift" = "Nicolas et al. downshift",
    #  "Nicolas et al.\n3'UTR end" = "Nicolas et al 3' UTR end",
    #  "Nicolas et al.\n5'UTR start" = "Nicolas et al 5' UTR start",
    #  'Nicolas et al.\nupshift' = 'Nicolas et al upshift'
    # )
    src = ifelse(
      str_detect(src, 'Nicolas'),
      'Nicolas et al.',
      src
    )
  ) %>%
  unique %>%
  bind_rows(
    dat %>%
      mutate(src = 'BSGatlas')
  ) %>%
  arrange(id) -> dat.full

dat.full %>%
  count(type, src) -> dat.count


dat.count %>%
  mutate(src = fct_reorder(src, n, max) %>%
           fct_rev) %>%
  mutate(lab = str_replace(n, '(\\d)(\\d{3})',  '\\1,\\2')) %>%
  mutate(type = ifelse(type == 'TSS',
                       'Transcription Start Sites (TSSs)',
                       'Transcription Termination Sites (TTSs)'
                       )) %>%
  ggplot(aes(x = src, y = n, fill = src, group = src)) +
  geom_bar(stat = 'identity',
           position = 'dodge') +
  facet_wrap(~ type) +
  geom_text(aes(label = lab),
            # size = 6,
            position = position_dodge(width=0.9),
            vjust=-0.25) +
  scale_fill_manual(values = ggsci::pal_jama()(5)[-2]) +
  theme_minimal(14) +
  scale_y_continuous(breaks = NULL) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank()
  ) +
  xlab(NULL) + ylab(NULL) +
  theme(legend.text = element_text(face = "italic"),
        legend.position = 'none',
        legend.justification = c(1, 1)) -> p1

###############################################################################
# comparison via venn diagrams
library(gridGraphics)

helper <- function(i) {
  dat.full %>%
    filter(src != 'BSGatlas') %>%
    filter(type == i) %>%
    mutate_at('src', fct_recode, 'Nicolas\net al.' = 'Nicolas et al.') %>%
    group_by(src) %>%
    do(i = list(.$id)) %>%
    with(set_names(map(i, 1), src)) %>%
    `[`(c('Nicolas\net al.', 'BsubCyc', 'DBTBS')) %>%
    venn(cexil = 1.1,
         cexsn = 1.1,
         opacity = 0.6,
         zcolor =  ggsci::pal_jama()(5)[-c(1, 2)])
  grid.echo()
  # grid.ls()
  grid.remove('graphics-plot-1-lines-1')
  size <- 15
  cowplot::ggdraw() + cowplot::draw_grob(
    grid.grab(),
    scale = 0.8
    # width = (size + 1)/2.54,
    # height = (size + 1)/2.54
  )
}

##############################################################################
# comparison of sigma factor distributions

# getting the raw dist from bsubcyc is annoying

bsub_raw$proteins %>%
  select(name = `COMMON-NAME`, other = SYNONYMS,
         box = `RECOGNIZED-PROMOTERS`) %>%
  drop_na(box) %>%
  separate_rows(box, sep = ';') -> bs.box

bsub_raw$promoters %>%
  select(pro.id = `UNIQUE-ID`, box = `PROMOTER-BOXES`) %>%
  separate_rows(box, sep = ';') %>%
  left_join(bs.box, 'box') %>%
  mutate(sigma = str_sub(other, -1)) %>%
  mutate_at('sigma', replace_na, '?') %>%
  transmute(id = pro.id, src = 'BsubCyc', sigma) -> bs

bind_rows(
  nicolas$upshifts %>%
    transmute(src = 'Nicolas et al.', id, sigma) %>%
    mutate_at('sigma', str_remove, 'Sig'),
  dbtbs$tss %>%
    transmute(src = 'DBTBS', id = 1:n(), sigma) %>%
    mutate_at('id', as.character) %>%
    mutate_at('sigma', str_remove, 'Sig'),
  bs,
  bsg.boundaries$TSS %>%
    transmute(id, sigma, src = 'BSGatlas')
) %>%
  mutate_at('sigma', replace_na, '?') %>%
  mutate_at('sigma', str_replace, '-', '?') %>%
  mutate_at('sigma', str_replace_all, '([A-Z])(?=[A-Z])', '\\1;') %>%
  separate_rows(sigma, sep = ';') %>%
  mutate(sigma = ifelse(sigma == 'C', 'YlaC', sigma)) %>%
  mutate(
    group = ifelse(sigma %in% c('?', 'A', 'E','F', 'G', 'K'),
                   sigma,
                   'other')
  ) -> tss.dat


tss.dat %>%
  count(src, group) %>%
  left_join(
    tss.dat %>%
      select(id, src) %>%
      # unique() %>%
      # nope, one position can have multiple binding sites
      count(src) %>%
      rename(total = n),
    'src'
  ) %>%
  mutate(prop = n / total * 100) -> tss.prop
  
tss.prop %>%
  mutate(group = fct_reorder(group, prop, max) %>% fct_rev) %>%
  mutate(group = fct_recode(group, 'unknown' = '?')) %>%
  ggplot(aes(x = src, fill = group, y = prop)) +
  geom_bar(stat = 'identity') +
  xlab(NULL) + ylab('Proprotion [%]') +
  ggsci::scale_fill_jco(name = 'Sigma Factor') +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank()
  ) -> p2
  

##############################################################################
# put it all togehter

cowplot::plot_grid(
  p1,
  cowplot::plot_grid(
    helper('TSS'), helper('TTS'), p2,
    # function() helper('TSS'),
    # partial(helper, 'TSS'), partial(helper, 'TTS'), p2,
    labels = c('(b) TSS', '(c) TTS', '(d)'),
    rel_widths = c(1, 1, 2),
    ncol = 3, hjust = 0
  ),
  labels = c('(a)', NULL),
  ncol = 1
)

ggsave('analysis/05b_bounds.pdf', 
       width = 30, height = 17, units = 'cm')

