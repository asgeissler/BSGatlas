#' Vectorized Helper for concaternating list/vecotr of entries to single string
#' identical(vec_paste(c('A', 'B', 'A'), '_'), 'A_B_A')
vec_paste <- function(v, sep = ';') {
  do.call(partial(paste, sep = sep), as.list(v))
}

#' Wrapper for vec_paste that replaces NA and at least returns an empty string
#' identical(clean_paste(c(NA, 'a', 'a', 'b')), 'a;b')
#' identical(clean_paste(c(NA, NA)), '')
clean_paste <-function(i, sep = ';') {
  res <- replace_na(i, '') %>%
    unique %>%
    discard(. == '') %>%
    unlist %>%
    vec_paste(sep = sep)
  if (identical(res, character(0))) {
    ''
  } else {
    res
  }
}

#' given one tibble merge rows with same id in column `col` with clean.join.vector
#' identical(
#'   merge_variants(
#'     tribble(
#'       ~foo, ~A, ~B,
#'       1, 'a', '',
#'       1, 'b', '',
#'       2, 'c', ''
#'     ),
#'     'foo'
#'   ),
#'   tribble(
#'       ~foo, ~A, ~B,
#'       1, 'a;b', '',
#'       2, 'c', ''
#'   )
#' )
merge_variants <- function(tbl, col, sep = ';') {
  tbl.group <- tbl %>%
    #group_by(eval(parse(text=col)))
    group_by(!! sym(col))
  tbl.group %>% summarise_all(clean_paste, sep = sep)
}



#' Given 2 tibbles `a` and `b` augment `a` with information from `b`.
#' Matching rows are identified by the specified named vector.
#'
#' @param a  tibble
#' @param b  tibble
#' @param by named string vector
#' @param sep separator char, default ';'
#'
#' @return tibble
#' @export
#'
#' @examples
#' library(tidyverse)
#' identical(
#'   augment(
#'     tribble(
#'       ~locus, ~meta,
#'       'a', '',
#'       'b', 'x'
#'     ),
#'     tribble(
#'       ~gene, ~meta, ~meta2,
#'       'a', 'x2', 'sRNA1',
#'       'b', 'x3',   'sRNA2'
#'     ),
#'     c('locus' = 'gene')),
#'   # expected outcome
#'   tribble(
#'     ~locus, ~meta, ~meta2,
#'     'a', 'x2',   'sRNA1',
#'     'b', 'x;x3', 'sRNA2'
#'   )
#' )
augment <- function(a, b, by, sep = ';') {
  # merge on named vector
  large <- left_join(a, b, by, suffix=c('_a', '_b'))
  # idenfify double entries
  select_vars(names(large), ends_with('_a')) %>%
    # extract original name
    str_replace('_a$', '') %>%
    map(1) %>%
    unlist %>%
    # merge for each double entry
    map(function(x) {
      a <- paste0(x, '_a')
      b <- paste0(x, '_b')
      tibble(
        !!x := map2(pull(large, a), pull(large, b), c) %>%
          map(clean_paste, sep = sep) %>%
          unlist
      )
    }) %>%
    # bind all tibbles with the large one
    c(large) %>%
    do.call(what = bind_cols) %>%
    # remove unjoined
    select(-ends_with('_a'), -ends_with('_b'))
}


#' Changed rows
#' How many rows have changed in b when merging over the given named by.
#' It does not state removal/insertion of rows
changed_rows <- function(a, b, by) {
  # similar strucutre to augment
  # shared primary keys
  large<- inner_join(a, b, by = by, suffix=c('_a', '_b'))
  # idenfify double entries and compare
  comparison <- select_vars(names(large), ends_with('_a')) %>%
    # extract original name
    str_replace('_a$', '') %>%
    unlist %>%
    # merge for each double entry
    map(function(x) {
      a <- paste0(x, '_a')
      b <- paste0(x, '_b')
      tibble(
        !!x := map2(pull(large, a), pull(large, b), magrittr::equals) %>%
          unlist
      )
    }) %>%
    do.call(what = bind_cols)
  # compute number
  pmap(comparison, any) %>%
    unlist %>%
    magrittr::not() %>%
    sum
}


#' Update primary key
#'
#' Replace the primary key in `tbl` and also update references in the `target`
#' column.
#'
#' tbl <- tribble(
#'   ~id, ~interaction,
#'   'A', 'B',
#'   'B', 'A;C',
#'   'C', NA
#' )
#' update <- tribble(
#'   ~id2, ~new,
#'   'B', 'flee',
#'   'C', 'glee'
#' )
#' result <- tribble(
#'   ~id, ~interaction,
#'   'A', 'flee',
#'   'flee', 'A;glee',
#'   'glee', ''
#' )
#' identical(
#'   update_key(tbl, 'id', update, c('id' = 'id2'), 'new', 'interaction'),
#'   result
#' )
update_key <- function(tbl, primary, update, by, newkey, target, sep = ';') {
  # newkey column is not ambiguoius
  assertthat::assert_that(! newkey %in% names(tbl))
  # convert to symbols
  primary <- sym(primary)
  newkey <- sym(newkey)
  target <- sym(target)
  # 1. substitute primary key
  tbl %>%
    left_join(update, by) %>%
    mutate(
      !! primary := ifelse(!is.na(!! newkey), !! newkey, !! primary)
    ) %>%
    # 2. replace target
    # 2a. spearate row
    separate_rows(!! target, sep = sep) %>%
    # 2b work
    # remove newkey from id update before joining again
    select(- !! newkey) %>%
    left_join(
      update,
      # name in by is primary key in tbl
      set_names(by, target)
    ) %>%
    mutate(
      !! target := ifelse(!is.na(!! newkey), !! newkey, !! target)
    ) %>%
    # 2c. collect results
    merge_variants(primary, sep) %>%
    # 3. make sure to return only tbl columns
    select(!!! syms(names(tbl)))
}
