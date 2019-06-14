#' extract sequence for GRanges and explicitly extract first/last tribblet
#' that might be a codon with respect to strandness
#'
#' Assumption: tibl is one-based coordinate
#'
#' @param seq DNAString
#' @param grange GRange with additional `id` column
extract_seq <- function(seq, grange) {
  seqs <- BSgenome::getSeq(seq, grange)
  return(tibble(
    id = grange$id,
    len = GenomicRanges::width(seqs),
    modulo_three = len %% 3 == 0,
    first_tribblet = GenomicRanges::narrow(seqs, 1, 3) %>% as.vector,
    # -2 not -3 because one based indexing
    last_tribblet = GenomicRanges::narrow(seqs, len - 2, len) %>% as.vector,
    full_sequence = as.vector(seqs)
  ))
}
#' Wrapper for tibbles
#'
#' Assumption: tibble is one-based coordinate
#'
#' @param seqs DNAString as +/- list
#' @param tibble annotation tibble
extract_seqs <- function(seqs, tbl, chr = 'NC_000964.3') {
  pos <- tbl %>%
    mutate(seqnames = 'chr') %>%
    plyranges::as_granges()
  lookup <- extract_seq(seqs, pos)

  return(left_join(tbl, lookup, by = 'id'))
}

# small test case
#1. without strandness, no wrapper
{
  s <- Biostrings::DNAStringSet('AACCTTGGAA')
  names(s) <- 'chr'
  annot <- tribble(
    ~id, ~start, ~end, ~strand,
    'alpha',     4,    9,     '*',
    'beta',      4,    8,     '*'
  ) %>%
    mutate(seqnames = 'chr') %>%
    plyranges::as_granges()
  expect <- tribble(
    ~id, ~len, ~modulo_three, ~first_tribblet, ~last_tribblet, ~full_sequence,
    'alpha',   6L,       TRUE,        'CTT',        'GGA',       'CTTGGA',
    'beta',    5L,      FALSE,        'CTT',        'TGG',       'CTTGG'
  )
  assertthat::are_equal(expect, extract_seq(s, annot))
}

#2. with strands and wrapper
{
  s <- Biostrings::DNAStringSet('AACCTTGGAA')
  names(s) <- 'chr'
  annot <- tribble(
    ~id, ~start, ~end, ~strand,
    'alpha',     4,    9,     '+',
    'beta',      4,    8,     '-'
  )
  expect <- tribble(
    ~id, ~len, ~modulo_three, ~first_tribblet, ~last_tribblet, ~full_sequence,
    'alpha',   6L,       TRUE,        'CTT',        'GGA',       'CTTGGA',
    'beta',    5L,      FALSE,        'CCA',        'AAG',       'CCAAG'
  )
  assertthat::are_equal(
    extract_seqs(s, annot, 'chr'),
    left_join(annot, expect, by = 'id')
  )
}


#' Find pattern
#' Find for a vector of string pattern all occurences in a given genome sequence
#' while also allowing possible mismatches.
#'
#' @param strings vector
#' @param mismatch integer
#' @param seq Biostrings object
#'
#' @return list of tibbles of counts and positions for each strand
#' @export
find_pattern <- function(strings, mismatch = 0, seq = bacillus_genome$seq) {
  # Explicitly check also the complementary part
  comp <- Biostrings::reverseComplement(seq)
  both <- Biostrings::DNAStringSet(list(seq, comp))
  names(both) <- c('+', '-')

  # Helper function
  match <- function(pattern, both, mismatch) {
    # for single string find occurence on both strands
    search <- Biostrings::vmatchPattern(pattern, both, max.mismatch = mismatch)
    # foreach strand
    map2(as.list(search), names(search), function(found, strand) {
      result <- list(
        # count occurences
        counts = tibble(
          sequence = pattern,
          strand = strand,
          count_matches = length(found)
        ),
        # remeber positions
        positions = if(length(found) != 0) {
          as.tibble(found) %>%
            transmute(sequence = pattern, strand = strand, start, end)
        } else {
          NULL
        }
      )
      # set negative strand positions to correctly match the forward
      if (strand == '-' && !is.null(result$positions)) {
        result$positions %<>% mutate(
          tmp = start,
          start = as.integer(length(seq) - end + 1),
          end = as.integer(length(seq) - tmp + 1)
        ) %>%
          select(- tmp)
      }
      return(result)
    })
  }


  cluster <- makeCluster(detectCores() - 1, type = 'FORK')
  result <- strings %>%
    parLapply(cl = cluster, fun = match, both = both, mismatch = mismatch)
  stopCluster(cluster)

  # combine results
  results <- invoke(c, result)
  # put counts and positions together
  list('counts', 'positions') %>%
    purrr::set_names(.) %>%
    map(function(i) {
      map(results, i) %>% bind_rows
    })
}

#find pattern counts and positions, without mismatches
{
  strings <- c('AAA', "ACA", 'GGG', 'TTT')
  seq <- Biostrings::DNAString('AAAGGG')
  
  expected <- list(
    counts = tribble(
      ~ sequence, ~ strand, ~ count_matches,
      'AAA', '+', 1L,
      'ACA', '+', 0L,
      'GGG', '+', 1L,
      'TTT', '+', 0L,
      'AAA', '-', 0L,
      'ACA', '-', 0L,
      'GGG', '-', 0L,
      'TTT', '-', 1L
    ),
    positions = tribble(
      ~ sequence, ~ strand, ~ start, ~ end,
      'AAA', '+', 1L, 3L,
      'GGG', '+', 4L, 6L,
      'TTT', '-', 1L, 3L
    )
  )
  assertthat::are_equal(
    find_pattern(strings, 0, seq),
    expected
  )
}

# find pattern positions with mismatches
{
  strings <- c('AAA')
  seq <- Biostrings::DNAString('AAAGGT')
  
  expected <- list(
    counts = tribble(
      ~ sequence, ~ strand, ~ count_matches,
      'AAA', '+', 5L,
      'AAA', '-', 3L
    ),
    positions = tribble(
      ~ sequence, ~ strand, ~ start, ~ end,
      'AAA', '+', -1L, 1L,
      'AAA', '+', 0L, 2L,
      'AAA', '+', 1L, 3L,
      'AAA', '+', 2L, 4L,
      'AAA', '+', 3L, 5L,
      'AAA', '-', 4L, 6L,
      'AAA', '-', 5L, 7L,
      'AAA', '-', 6L, 8L
    )
  )
  assertthat::are_equal(
    find_pattern(strings, 2, seq),
    expected
  )
}

