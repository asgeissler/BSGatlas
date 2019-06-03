# The biocyc database BsubCyc requires a special parser

# parse a list of lines belonging to a single entry
# xs should not contain comments
# expect_equal(
#   parse_entry(c('A - B', '/FOO', 'C - D', 'C - E')),
#   c(A = 'B\nFOO', C = 'D;E')
# )
library(hashmap)
parse_entry <- function(xs) {
  if(length(xs) == 0) return(list())
  res <- hashmap('', '')
  res$clear()
  lastKey <- NA
  lastValue <- NA


  # helper function to populate `res`
  add2res <- function(x, y) {
    if (res$has_key(x)) {
      # 1. that key occured before
      res[[x]] <- paste(res[[x]], y, sep=';')
    } else {
      # 2. create new entry
      res[[x]] <- y
    }
  }

  for (i in xs) {
    if (startsWith(i, '/')) {
      # 1. line extends information from last entry
      lastValue <- paste(lastValue, substring(i, 2), sep='\n')
    } else {
      c(k, v) %<-% str_split(i, ' - ', n=2)[[1]]
      # safe last k, v tupple
      # check na necessary for first iteration
      if(!is.na(lastKey)) {
        add2res(lastKey, lastValue)
      }
      # remember new tupple
      lastKey <- k
      lastValue <- v
    }
  }
  # safe last k, v tupple
  add2res(lastKey, lastValue)
  return(res$data())
}


# helper predicate for comments and meta-information
is_comment <- function(line){
  c('#', '^') %>% map(startsWith, x=line) %>% unlist %>% any
}


# Parser for *.dat files from BioCyc
parse_dat <- function(path) {
  lines <- read_lines(path) %>% discard(.p = is_comment)
  # was helpful for debugging
  #lines <- read_lines(path, n_max = 150) %>% discard(.p=is.comment)

  newpos <- which(lines == '//')
  # remove the entry at the very end
  newpos <- newpos[1:length(newpos) - 1]
  # intervals for each entry
  map2(c(1, newpos + 1), c(newpos - 1, length(lines) - 1), seq) %>%
    # partition
    map(function(i) lines[i]) -> entries

  # parsing
  cluster <- makeCluster(detectCores() - 1, type = 'FORK')
  entries %<>% parLapply(cl = cluster, fun = parse_entry)
  stopCluster(cluster)

  # merge rows with filling NA's
  return(do.call(bind_rows, entries))
}

