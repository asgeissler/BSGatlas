#' Helper for finding all pairs of cloest neighbors,
#' including pairs from opposite strands for
#' intervals between two sets.
#' Overlapping pairs are indicated with zeros.
#' 
#' The possible neighbor on the opposite strand is only shown if it
#' is closer than the nearest same strand neighbor in upstream/down-stream
#' direction
#' 
#' The up/downstream indication is relative to th 5' end position,
#' equality is expclitily indicated.
#' For more details on the overlaps refer to `overlap_matching`
#' 
#' If for a gene no such partner exists, then this is indicated with NAs.
#' 
#' For each pair it is indicated if the x is in direction
#' of transcription upstream of y.
#' 
distance_matching <- function(x.tbl, y.tbl) {
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
  
  # Get upstream/downstream pair each with possible opposite strand
  crossing(
    # explicitly include `nearest`, otherwise it would not list all overlapping
    f = c(GenomicRanges::follow, GenomicRanges::precede, GenomicRanges::nearest), 
    strand = c(TRUE, FALSE)
  ) %>%
    # do the work
    rowwise() %>%
    do(res = plyranges:::make_hits(
      x.range, y.range, .$f, select = "all",
      ignore.strand = .$strand
    ) %>%
      plyranges:::expand_by_hits(
        x = x.range, y = y.range,
        suffix = c(".x", ".y")
      )
    ) %>%
    # get the pair
    pull(res) %>%
    map(as_tibble) %>%
    bind_rows %>%
    select(id.x, id.y) %>%
    unique %>%
    # add 5/3' positions
    left_join(x.tbl %>% set_names(paste0(names(.), '.x')), 'id.x') %>%
    left_join(y.tbl %>% set_names(paste0(names(.), '.y')), 'id.y') %>%
    transmute(
      x = id.x, y = id.y,
      mode = case_when(
        `5prime.x` == `5prime.y` ~ '5prime.equal',
        (strand.x == '+') & (`5prime.x` > `5prime.y`) ~ 'x.after.y',
        (strand.x == '-') & (`5prime.x` > `5prime.y`) ~ 'x.before.y',
        (strand.x == '+') & (`5prime.x` < `5prime.y`) ~ 'x.before.y',
        (strand.x == '-') & (`5prime.x` < `5prime.y`) ~ 'x.after.y',
        TRUE ~ 'case.error'
      ),
      antisense = strand.x != strand.y,
      distance = pmax(
       pmax(start.x - end.y, 0),
       pmax(start.y - end.x, 0)
      )
    )
}


{
  test_x <- tribble(
    ~id, ~start, ~end, ~strand,
    'A', 5, 10, '+'
  )
  test_y <- tribble(
    ~id, ~start, ~end, ~strand,
    'a', 15, 20, '+',
    'b',  1,  3, '+',
    'c',  1,  8, '+',
    'd',  1,  8, '-',
    'e',  5,  9, '+'
  )
  test_res <- tribble(
    ~x,  ~mode,          ~y, ~antisense, ~distance,
    'A', 'x.before.y',   'a', FALSE, 5,
    'A', 'x.after.y',    'b', FALSE, 2,
    'A', 'x.after.y',    'c', FALSE, 0,
    'A', 'x.before.y',   'd',  TRUE, 0,
    'A', '5prime.equal', 'e', FALSE, 0
  )
  testthat::expect_equal(
    distance_matching(test_x, test_y) %>%
      arrange(y) %>%
      select(!!! names(test_res)),
    test_res
  )
}
