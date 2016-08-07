# Internal OneR functions

# modified cut function for ensuring consistency of cut points and chosen cut points
# http://stackoverflow.com/questions/37899503/inconsistent-behaviour-of-cut-different-intervals-with-same-number-and-same-d
CUT <- function(x, breaks, ...) {
  if (length(breaks) == 1L) {
    nb <- as.integer(breaks + 1)
    dx <- diff(rx <- range(x, na.rm = TRUE))
    if (dx == 0) {
      dx <- abs(rx[1L])
      breaks <- seq.int(rx[1L] - dx/1000, rx[2L] + dx/1000, length.out = nb)
    } else {
      breaks <- seq.int(rx[1L], rx[2L], length.out = nb)
      breaks[c(1L, nb)] <- c(rx[1L] - dx/1000, rx[2L] + dx/1000)
    }
  }
  breaks.f <- c(breaks[1], as.numeric(formatC(0 + breaks[2:(length(breaks)-1)], digits = 3, width = 1L)), breaks[length(breaks)])
  cut(x, breaks = unique(breaks.f), ...)
}

mode <- function(x) {
  names(sort(-table(x[ , ncol(x)])))[1]
}

ADDNA <- function(x) {
  if (is.factor(x) & !("NA" %in% levels(x))) x <- factor(x, levels = c(levels(x), "NA"))
  x[is.na(x)] <- "NA"
  return(x)
}

add_range <- function(x, midpoints) {
  c(min(x, na.rm = TRUE) - 1/1000 * diff(range(x, na.rm = TRUE)), midpoints, max(x, na.rm = TRUE) + 1/1000 * diff(range(x, na.rm = TRUE)))
}

get_breaks <- function(x) {
  x <- x[x != "NA"]
  lower = as.numeric(sub("\\((.+),.*", "\\1", x))
  upper = as.numeric(sub("[^,]*,([^]]*)\\]", "\\1", x))
  breaks <- unique(na.omit(c(lower, upper)))
  return(breaks)
}

#' @importFrom stats na.omit
#' @importFrom stats filter
naive <- function(x, target) {
  orig <- x
  tmp <- na.omit(cbind(x, target))
  x <- tmp[ , 1]; target <- tmp[ , 2]
  xs <- split(x, target)
  midpoints <- sort(sapply(xs, mean, na.rm = TRUE))
  # The cutpoints are the means of the expected values of the respective target levels.
  breaks <- add_range(x, na.omit(filter(midpoints, c(1/2, 1/2))))
  CUT(orig, breaks = unique(breaks))
}

#' @importFrom stats coef
#' @importFrom stats glm
#' @importFrom stats binomial
logreg_midpoint <- function(data) {
  df <- data.frame(x = unlist(data), target = factor(rep(names(data), sapply(data, length))))
  coefs <-  suppressWarnings(coef(glm(target ~ x, data = df, family = binomial)))
  midpoint <- - coefs[1] / coefs[2]
  # test limits
  range <- sort(sapply(data, mean, na.rm = TRUE))
  if (length(range) == 1) range <- c(range, range)
  if (is.na(midpoint)) return(mean(range, na.rm = TRUE))
  if (midpoint < range[1]) return(range[1])
  if (midpoint > range[2]) return(range[2])
  # ---
  return(midpoint)
}

#' @importFrom stats na.omit
logreg <- function(x, target) {
  orig <- x
  tmp <- na.omit(cbind(x, target))
  x <- tmp[ , 1]; target <- tmp[ , 2]
  xs <- split(x, target)
  midpoints <- sapply(xs, mean)
  nl <- xs[order(midpoints)]
  pairs <- matrix(c(1:(length(nl) - 1), 2:length(nl)), ncol = 2, byrow = TRUE)
  midpoints <- apply(pairs, 1, function(x) logreg_midpoint(c(nl[x[1]], nl[x[2]])))
  breaks <- add_range(x, na.omit(midpoints))
  CUT(orig, breaks = unique(breaks))
}

entropy <- function(x) {
  freqs <- table(x) / length(x)
  - sum(freqs * log2(freqs))
}

#' @importFrom stats na.omit
infogain_midpoint <- function(data) {
  df <- data.frame(numvar = unlist(data), target = factor(rep(names(data), sapply(data, length))))
  data <- na.omit(df[order(df[ , 1]), ])
  numvar <- data$numvar; target <- data$target
  # determine midpoint candidates
  left_thresholds <- which(as.logical(diff(as.numeric(target))))
  midpoints <- (numvar[left_thresholds] + numvar[(left_thresholds + 1)]) / 2
  # calculate average entropies for all midpoint candidates
  belows <- sapply(midpoints, function(x) as.character(data[numvar <= x, 2]))
  aboves <- sapply(midpoints, function(x) as.character(data[numvar > x, 2]))
  below_entropies <- sapply(belows, function(x) length(x)/length(target) * entropy(x))
  above_entropies <- sapply(aboves, function(x) length(x)/length(target) * entropy(x))
  # calculate information gains and chose highest one
  infogains <- entropy(target) - (below_entropies + above_entropies)
  midpoints[which.max(infogains)]
}

#' @importFrom stats na.omit
infogain <- function(x, target) {
  orig <- x
  tmp <- na.omit(cbind(x, target))
  x <- tmp[ , 1]; target <- tmp[ , 2]
  xs <- split(x, target)
  midpoints <- sapply(xs, mean)
  nl <- xs[order(midpoints)]
  pairs <- matrix(c(1:(length(nl) - 1), 2:length(nl)), ncol = 2, byrow = TRUE)
  midpoints <- apply(pairs, 1, function(x) infogain_midpoint(c(nl[x[1]], nl[x[2]])))
  breaks <- add_range(x, na.omit(midpoints))
  CUT(orig, breaks = unique(breaks))
}
