# Helper to find in a directed graph all paths between the set
# of nodes with only out nodes and the set of only in nodes

library(tidygraph)

find_paths <- function(grph) {
  print('Finding Paths:')
  grph %>%
    activate(nodes) %>%
    mutate(in.degree = centrality_degree(mode = 'in'),
           out.degree = centrality_degree(mode = 'out')) %>%
    as_tibble -> grph.degrees
  
  starts <- grph.degrees %>%
    filter(out.degree > 0, in.degree == 0) %>%
    pull(name)
  ends <- grph.degrees %>%
    filter(out.degree == 0, in.degree > 0) %>%
    pull(name)
  
  worker <- function(grph, i, ends) {
    igraph::all_simple_paths(grph, i, ends) %>%
      map(names) %>%
      map(as_tibble)
  }
  
  cores <- makeForkCluster(detectCores() - 1)
  res <- parLapply(cores, starts,
                   safely(partial(worker, grph = grph, ends =ends)))
  res %>%
    map('result') %>%
    invoke(.f = c) %>%
    map2(1:length(.), ~ mutate(.x, path = .y)) %>%
    bind_rows -> tbl
  
  stopCluster(cores)
  
  print('Possible errors were:')
  res %>%
    map('error') %>%
    discard(is.null)
  
  # return the combined result
  return(tbl)
}



# # Small toy example
# nodes <- tribble(
#   ~id, ~name,
#   1, 'a',
#   2, 'b',
#   3, 'c',
#   4, 'd',
#   5, 'e',
#   6, 'start.a',
#   7, 'start.c',
#   8, 'start.d',
#   9, 'end.b',
#   10, 'end.c',
#   11, 'end.d'
# )
# edges <- tribble(
#   ~from, ~to,
#   'start.a', 'a',
#   'a', 'b',
#   'b', 'end.b',
#   'start.c', 'c',
#   'c', 'd',
#   'start.d', 'd',
#   'c', 'end.c',
#   'd', 'end.d'
# ) %>%
#   mutate(row = 1:n()) %>%
#   gather('key', 'value', from, to) %>%
#   left_join(nodes, c('value' = 'name')) %>%
#   select(-value) %>%
#   spread(key, id) %>%
#   select(from, to)
# 
# tbl_graph(nodes, edges, directed = TRUE) -> grph
# find_paths(grph)
#
# # Exptected paths:
# # start.a -> a -> b -> end.b
# # start.c -> c -> end.c
# # start.c -> c -> d -> end.d
# # start.d -> d -> end.d
# 
# 
