#' Helper for finding all overlapping pairs 
#' of intervals between two sets
#' 
#' The mode describes the specific type of overlap types, and the
#' notation given in the output.
#' Here X, Y are the ids from the input
#' 
#' 1. Positions are equal, but possibly anti-sense
#'    Notation: equal
#' 2. X 3' end overlaps Y 5' end
#'    implies both are on the same strand
#'    Notation: 3-5_overlap
#'    Picto:   5'>XXXXX>3'
#'                5'>YYYYY>3'
#' 3. X 3' end overlaps Y 3' end
#'    implies oppposite strand
#'    Notation: 3-3_overlap
#'    Picto:   5'>XXXXX>3'
#'                3'<YYYYY<5'
#' 4. X 5' end overlaps Y 3' end
#'    implies same strand
#'    Notation: 5-3_overlap
#'    Picto:   3'<XXXXX<5'
#'                3'<YYYYY<5'
#' 5. X 5' end overlaps Y 5' end
#'    implies oppposite strand
#'    Notation: 5-5_overlap
#'    Picto:   3'<XXXXX<5'
#'                5'>YYYYY>3'
#' 6. X fully contained in Y, but possibly antisense
#'    Notation: contianed_by
#'    Picto:    5'>XXX>3'
#'             5'>YYYYYY>3'
#'    or:       5'>XXX>3'
#'             3'<YYYYYY<5'
#' Analogous case:
#' 7. contains
#' Special cases of 6/7 are
#'  a) start positions are (nearly) equal
#'  b) end positions are (nearly) equal
#' These are indicated by additional columns that indicate the distances of
#'  * the X5' to either the Y5' or in the antisense case the Y3'
#'  * the X3' to either the Y3' or in the antisense case the Y5'
#'  
#'  8' Either X or Y has no overlap, in the output the id of the missing partner
#'     is set to NA.
#'     Notation: without_overlap
#'     !Note: Each element X/Y could have two 'without_overlap', for both
#'            sense and antisense case

overlap_matching <- function(x.tbl, y.tbl) {
  x.tbl %<>% mutate(
    seqnames = 'placeholder',
    length = end - start + 1,
    `5prime` = ifelse(strand == '+', start, end),
    `3prime` = ifelse(strand == '+', end, start)
  )
  y.tbl %<>% mutate(
    seqnames = 'placeholder',
    length = end - start + 1,
    `5prime` = ifelse(strand == '+', start, end),
    `3prime` = ifelse(strand == '+', end, start)
  )
  x.range <- plyranges::as_granges(x.tbl)
  y.range <- plyranges::as_granges(y.tbl)
  
  # Find and ocmpute overlaps
  res <- x.range %>%
    plyranges::join_overlap_intersect(y.range) %>%
    as_tibble %>%
    select(x = id.x, y = id.y, overlap = width) %>%
    # add start / end positions
    left_join(
      x.tbl %>%
        set_names(paste0('x.', names(.))),
      c('x' = 'x.id')
    ) %>%
    left_join(
      y.tbl %>%
        set_names(paste0('y.', names(.))),
      c('y' = 'y.id')
    ) %>%
    mutate(
      antisense = x.strand != y.strand
    )
  # Compute distances as described in function explanation
  res %<>%
    mutate(
      # the X5' to either the Y5' or in the antisense case the Y3'
      x5.dist = abs(x.5prime - ifelse(antisense, y.3prime, y.5prime)),
      # the X3' to either the Y3' or in the antisense case the Y5'
      x3.dist = abs(x.5prime - ifelse(antisense, y.5prime, y.3prime))
    )
  # Classify
  res %<>% mutate(
    mode = case_when(
      # 1. Positions are equal, but possibly anti-sense
      (x.start == y.start) & (x.end == y.end) ~ 'equal',
      # 6. X fully contained in Y, but possibly antisense
      (x.start >= y.start) & (x.end <= y.end) ~ 'contained_by',
      # 7. contains
      (x.start <= y.start) & (x.end >= y.end) ~ 'contains',
      # 2. X 3' end overlaps Y 5' end
      # 3. X 3' end overlaps Y 3' end
      # 4. X 5' end overlaps Y 3' end
      # 5. X 5' end overlaps Y 5' end
      (x.start < y.start) & (x.end < y.end) &
        (x.strand == '+') & (y.strand == '+') ~ '3-5_overlap',
      (x.start < y.start) & (x.end < y.end) &
        (x.strand == '-') & (y.strand == '-') ~ '5-3_overlap',
      (x.start < y.start) & (x.end < y.end) &
        (x.strand == '+') & (y.strand == '-') ~ '3-3_overlap',
      (x.start < y.start) & (x.end < y.end) &
        (x.strand == '-') & (y.strand == '+') ~ '5-5_overlap',
      (x.start > y.start) & (x.end > y.end) &
        (x.strand == '+') & (y.strand == '+') ~ '5-3_overlap',
      (x.start > y.start) & (x.end > y.end) &
        (x.strand == '-') & (y.strand == '-') ~ '3-5_overlap',
      (x.start > y.start) & (x.end > y.end) &
        (x.strand == '+') & (y.strand == '-') ~ '5-5_overlap',
      (x.start > y.start) & (x.end > y.end) &
        (x.strand == '-') & (y.strand == '+') ~ '3-3_overlap'
    )
  )
  
  # Add indications of non-overlap
  x.sub <- x.tbl %>%
    set_names(paste0('x.', names(.))) %>%
    rename(x = x.id)
  y.sub <- y.tbl %>%
    set_names(paste0('y.', names(.))) %>%
    rename(y = y.id)
  all.sub <- bind_rows(x.sub, y.sub) %>%
    crossing(antisense = c(TRUE, FALSE)) %>%
    mutate(mode = 'without_overlap') %>%
    anti_join(res, c('x', 'antisense')) %>%
    anti_join(res, c('y', 'antisense'))
  
  res %<>% bind_rows(all.sub)
  
  # Invariant: distances and overlap should add up
  assertthat::are_equal(
    0,
    res %>%
      filter(2 * overlap != X.length + Y.length - X3.dist - X.5dist) %>%
      drop_na %>%
      nrow
  )
  return(res)
}


{
  test_x <- tribble(
    ~id, ~start, ~stop, ~strand,
    'A',  1,  5, '+',
    'B', 20, 30, '+',
    'C', 25, 35, '-',
    'D', 40, 50, '+'
  )
  test_y <- tribble(
    ~id, ~start, ~stop, ~strand,
    'a',  3,  7, '+',
    'b', 27, 28, '+',
    'c', 20, 30, '+',
    'd',  6,  9, '-'
  )
  test_res <- tribble(
    ~x,  ~y,  ~mode,          ~antisense, ~x5.dist, ~x3.dist, ~overlap, ~y.length, ~y.length, ~jaccard,
    'A', 'a', "3-5_overlap",     FALSE,           2,       2,        3,         5,         5,  0.429,
    'B', 'b', "contains",        FALSE,           7,       2,        2,        11,         2,  0.18,
    'B', 'c', "equal",            TRUE,           0,       0,       11,        11,        11,   1,
    'C', 'b', "5-5_overlap",      TRUE,           7,       2,        2,        11,         2,  0.18,
    'C', 'c', "3-3_overlap",      TRUE,           5,       5,        6,        11,        11,  0.375,
    'A',  NA, "without_overlap",  TRUE,          NA,      NA,       NA,         5,        NA,   NA,
    NA,  'a', "without_overlap",  TRUE,          NA,      NA,       NA,        NA,         5,  NA,
    NA,  'c', "without_overlap", FALSE,          NA,      NA,       NA,        NA,        10,  NA,
    NA,  'd', "without_overlap",  TRUE,          NA,      NA,       NA,        NA,         4,  NA,
    NA,  'd', "without_overlap", FALSE,          NA,      NA,       NA,        NA,         4,  NA
  )
}