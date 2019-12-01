makeLegend <- function(descr) {
  # construct legend (separate long names with newline)
  numLgnd <- length(dscr)
  lgnd <- c(numLgnd)
  for (i in seq(dscr)){
    title <- dscr[i]
    # if title is too long
    if (nchar(title) > 20) {
      # get the location of white space
      space <- unlist(gregexpr(" ", title))
      # if there is just one space in title,
      # split location is the first space
      if (length(space) == 1) {
        split <- space[1]
      } else {
        # otherwise split location is the second space
        split <- space[2]}
      f <- substring(dscr[i], 
                     1, split-1)
      s <- substring(dscr[i],
                     split+1)
      title <- paste(f, s, sep = "\n")
    }
    # store the title into lgnd
    lgnd[i] <- title
  }
}