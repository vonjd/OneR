mode <- function(x) {
  names(sort(-table(x[ , ncol(x)])))[1]
}

addNA <- function(x) {
  if (is.factor(x)) x <- factor(x, levels = c(levels(x), "NA"))
  x[is.na(x)] <- "NA"
  return(x)
}

get_breaks <- function(x) {
  lower = as.numeric(sub("\\((.+),.*", "\\1", x))
  upper = as.numeric(sub("[^,]*,([^]]*)\\]", "\\1", x))
  breaks <- unique(c(lower, upper))
  return(breaks)
}

naive <- function(x, target) {
  orig <- x
  tmp <- na.omit(cbind(x, target))
  x <- tmp[ , 1]; target <- tmp[ , 2]
  xs <- split(x, target)
  midpoints <- sort(sapply(xs, mean))
  # The cutpoints are the means of the expected values of the respective target levels.
  breaks <- c(min(x) - 1/1000 * diff(range(x)), na.omit(filter(midpoints, c(1/2, 1/2))), max(x) + 1/1000 * diff(range(x)))
  #breaks <- unique(as.numeric(formatC(0 + breaks, digits = 3, width = 1L)))
  cut(orig, breaks = unique(breaks))
}

logreg_midpoint <- function(data) {
  df <- data.frame(x = unlist(data), target = factor(rep(names(data), sapply(data, length))))
  coefs <-  suppressWarnings(coef(glm(target ~ x, data = df, family = binomial)))
  midpoint <- - coefs[1] / coefs[2]
  # test limits
  range <- sort(sapply(data, mean))
  if (is.na(midpoint)) return(mean(range))
  if (midpoint < range[1]) return(range[1])
  if (midpoint > range[2]) return(range[2])
  # ---
  return(midpoint)
}

logreg <- function(x, target) {
  orig <- x
  tmp <- na.omit(cbind(x, target))
  x <- tmp[ , 1]; target <- tmp[ , 2]
  xs <- split(x, target)
  midpoints <- sapply(xs, mean)
  nl <- xs[order(midpoints)]
  pairs <- matrix(c(1:(length(nl) - 1), 2:length(nl)), ncol = 2, byrow = TRUE)
  midbreaks <- apply(pairs, 1, function(x) logreg_midpoint(c(nl[x[1]], nl[x[2]])))
  breaks <- c(min(x) - 1/1000 * diff(range(x)), midbreaks, max(x) + 1/1000 * diff(range(x)))
  #breaks <- unique(as.numeric(formatC(0 + breaks, digits = 3, width = 1L)))
  cut(orig, breaks = unique(breaks))
}

#' Binning function
#'
#' Discretizes all numerical data in a dataframe into categorical bins of equal length or content.
#' @param data dataframe which contains the data.
#' @param nbins number of bins (= levels).
#' @param labels character vector of labels for the resulting category.
#' @param method a character string specifying the binning method, see 'Details'; can be abbreviated.
#' @param na.omit boolean value whether instances with missing values should be removed.
#' @keywords binning discretization discretize
#' @details Character strings and logical strings are coerced into factors. Matrices are coerced into dataframes. When called with a single vector only the respective factor (and not a dataframe) is returned.
#' Method \code{"length"} gives intervals of equal length, method \code{"content"} gives intervals of equal content (via quantiles).
#'
#' When \code{"na.omit = FALSE"} a new level \code{"NA"} is introduced into each factor.
#' @author Holger von Jouanne-Diedrich, \email{r-project@ephorie.de}
#' @references \url{http://vonjd.github.io/OneR/}
#' @seealso \code{\link{OneR}}, \code{\link{optbin}}
#' @examples
#' data <- iris
#' str(data)
#' str(bin(data))
#' str(bin(data, nbins = 3))
#' str(bin(data, nbins = 3, labels = c("small", "medium", "large")))
#' @export
bin <- function(data, nbins = 5, labels = NULL, method = c("length", "content"), na.omit = TRUE) {
  method <- match.arg(method)
  vec <- FALSE
  if (is.atomic(data) == TRUE & is.null(dim(data)) == TRUE) { vec <- TRUE; data <- data.frame(data) }
  # could be a matrix -> dataframe (even with only one column)
  if (is.list(data) == FALSE) data <- data.frame(data)
  if (na.omit == TRUE) {
    len_rows <- nrow(data)
    data <- na.omit(data)
    if (len_rows > nrow(data)) warning("at least one instance was removed due to missing values")
  }
  if (!is.null(labels)) if (nbins != length(labels)) stop("number of 'nbins' and 'labels' differ")
  if (nbins <= 1) stop("number of 'bins' must be bigger than 1")
  data[] <- lapply(data, function(x) if (is.numeric(x) ) {
    if (length(unique(x)) <= nbins) as.factor(x)
    else {
      if (method == "content") nbins <- c(min(x) - 1/1000 * diff(range(x)), na.omit(quantile(x, (1:(nbins-1)/nbins))), max(x) + 1/1000 * diff(range(x)))
      #if (method == "content") { nbins <- c(min(x) - 1/1000 * diff(range(x)), na.omit(quantile(x, (1:(nbins-1)/nbins))), max(x) + 1/1000 * diff(range(x))); nbins <- unique(as.numeric(formatC(0 + nbins, digits = 3, width = 1L))) }
      cut(x, breaks = nbins, labels = labels)
    }
  } else as.factor(x))
  if (na.omit == FALSE) {
    # data[] <- lapply(data, addNA)
    data[] <- lapply(data, function(x) if(any(is.na(x))) addNA(x) else x)
  }
  if (vec) { data <- unlist(data); names(data) <- NULL }
  return(data)
}

#' Optimal Binning function
#'
#' Discretizes all numerical data in a dataframe into categorical bins where the cut points are optimally aligned with the target categories, thereby a factor is returned.
#' When building a OneR model this could result in fewer rules with enhanced accuracy.
#' @param data dataframe which contains the data. When \code{formula = NULL} (the default) the last column must be the target variable.
#' @param formula formula interface for the \code{optbin} function.
#' @param method a character string specifying the method for optimal binning, see 'Details'; can be abbreviated.
#' @param na.omit boolean value whether instances with missing values should be removed.
#' @keywords binning discretization discretize
#' @details The cutpoints are calculated by pairwise logistic regressions (method \code{"logreg"}) or as the means of the expected values of the respective classes (\code{"naive"}).
#' The function is likely to give unsatisfactory results when the distributions of the respective classes are not (linearly) separable. Method \code{"naive"} should only be used when distributions are (approximately) normal,
#' although in this case \code{"logreg"} should give comparable results, so it is the preferable (and therefore default) method.
#'
#' Character strings and logical strings are coerced into factors. Matrices are coerced into dataframes. If the target is numeric it is turned into a factor with the number of levels equal to the number of values. Additionally a warning is given.
#'
#' When \code{"na.omit = FALSE"} a new level \code{"NA"} is introduced into each factor.
#' @author Holger von Jouanne-Diedrich, \email{r-project@ephorie.de}
#' @references \url{http://vonjd.github.io/OneR/}
#' @seealso \code{\link{OneR}}, \code{\link{bin}}
#' @examples
#' data <- iris # without optimal binning
#' model <- OneR(data, verbose = TRUE)
#' summary(model)
#'
#' data_opt <- optbin(iris) # with optimal binning
#' model_opt <- OneR(data_opt, verbose = TRUE)
#' summary(model_opt)
#'
#' ## The same with the formula interface:
#' data_opt <- optbin(formula = Species ~., data = iris)
#' model_opt <- OneR(data_opt, verbose = TRUE)
#' summary(model_opt)
#'
#' @export
optbin <- function(data, formula = NULL, method = c("logreg", "naive"), na.omit = TRUE) {
  method <- match.arg(method)
  if (class(formula) == "formula") {
    mf <- model.frame(formula = formula, data = data, na.action = NULL)
    data <- mf[c(2:ncol(mf), 1)]
  } else if (is.null(formula) == FALSE) stop("invalid formula")
  if (is.list(data) == FALSE) data <- data.frame(data)
  if (dim(data)[2] < 2) stop("data must have at least two columns")
  if (is.numeric(unlist(data[ncol(data)])) == TRUE) {
    data[ncol(data)] <- as.factor(unlist(data[ncol(data)]))
    warning("target is numeric")
  }
  if (na.omit == TRUE) {
    len_rows <- nrow(data)
    data <- na.omit(data)
    if (len_rows > nrow(data)) warning("at least one instance was removed due to missing values")
  } else {
    # only add NA to target
    if(any(is.na(unlist(data[ncol(data)])))) data[ncol(data)] <- addNA(unlist(data[ncol(data)]))
  }
  target <- data[ncol(data)]
  nbins <- length(unique(unlist(target)))
  if (nbins <= 1) stop("number of target levels must be bigger than 1")
  data[] <- lapply(data, function(x) if (is.numeric(x)) {
    if (length(unique(x)) <= nbins) as.factor(x) else do.call(method, list(x, target))
  } else as.factor(x))
  if (na.omit == FALSE) {
    # data[] <- lapply(data, addNA)
    data[] <- lapply(data, function(x) if(any(is.na(x))) addNA(x) else x)
  }
  return(data)
}

#' Remove factors with too many levels
#'
#' Removes all columns of a dataframe where a factor (or character string) has more than a maximum number of levels.
#' @param data dataframe which contains the data.
#' @param maxlevels number of maximum factor levels.
#' @param na.omit boolean value whether missing values should be treated as a level, defaults to omit missing values before counting.
#' @details Often categories that have very many levels are not useful in modelling OneR rules because they result in too many rules and tend to overfit.
#' Examples are IDs or names.
#'
#' Character strings are treated as factors although they keep their datatype. Numeric data is left untouched.
#' @author Holger von Jouanne-Diedrich, \email{r-project@ephorie.de}
#' @references \url{http://vonjd.github.io/OneR/}
#' @seealso \code{\link{OneR}}
#' @examples
#' df <- data.frame(numeric = c(1:26), alphabet = letters)
#' str(df)
#' str(maxlevels(df))
#' @export
maxlevels <- function(data, maxlevels = 20, na.omit = TRUE) {
  if (is.list(data) == FALSE) stop("data must be a dataframe")
  if (maxlevels <= 2) stop("maxlevels must be bigger than 2")
  tmp <- bin(data, nbins = 2, na.omit = na.omit)
  nlevels <- sapply(tmp, nlevels)
  cols <- nlevels <= maxlevels
  return(data[cols])
}

#' Predict method for OneR models
#'
#' Predict values based on OneR model object.
#' @param object object of class \code{"OneR"}.
#' @param newdata dataframe in which to look for the feature variable with which to predict.
#' @param ... further arguments passed to or from other methods.
#' @details \code{newdata} can have the same format as used for building the model but must at least have the feature variable that is used in the OneR rules.
#' If cases appear that were not present when building the model the predicted value is \code{UNSEEN}.
#' @author Holger von Jouanne-Diedrich, \email{r-project@ephorie.de}
#' @references \url{http://vonjd.github.io/OneR/}
#' @seealso \code{\link{OneR}}
#' @examples
#' model <- OneR(iris)
#' prediction <- predict(model, iris[1:4])
#' eval_model(prediction, iris[5])
#' @export
predict.OneR <- function(object, newdata, ...) {
  if (is.list(newdata) == FALSE) stop("newdata must be a dataframe")
  if (all(names(newdata) != object$feature)) stop("cannot find feature column in newdata")
  model <- object
  data <- newdata
  index <- which(names(data) == model$feature)[1]
  if (is.numeric(data[ , index])) {
    levels <- names(model$rules)
    if (grepl(",", levels[1]) == TRUE) {
      features <- as.character(cut(unlist(data[ , index]), breaks = get_breaks(levels)))
    } else features <- as.character(data[ , index])
  } else features <- as.character(data[ , index])
  features[is.na(features)] <- "NA"
  return(sapply(features, function(x) if (is.null(model$rules[[x]]) == TRUE) "UNSEEN" else model$rules[[x]]))
}

#' Summarize OneR models
#'
#' \code{summary} method for class \code{OneR}.
#' @param object object of class \code{"OneR"}.
#' @param ... further arguments passed to or from other methods.
#' @details Prints the rules of the OneR model, the accuracy, a contingency table of the feature attribute and the target and performs a chi-squared test on this table.
#'
#' In the contingency table the maximum values in each column are highlighted by adding a '*', thereby representing the rules of the OneR model.
#' @author Holger von Jouanne-Diedrich, \email{r-project@ephorie.de}
#' @references \url{http://vonjd.github.io/OneR/}
#' @seealso \code{\link{OneR}}
#' @keywords diagnostics
#' @examples
#' model <- OneR(iris)
#' summary(model)
#' @export
summary.OneR <- function(object, ...) {
  model <- object
  print(model)
  tbl <- model$cont_table
  pos <- cbind(apply(tbl, 2, which.max), 1:dim(tbl)[2])
  tbl <- addmargins(tbl)
  tbl[pos] <- paste("*", tbl[pos])
  cat("Contingency table:\n")
  print(tbl, quote = FALSE, right = TRUE)
  cat("---\nMaximum in each column: '*'\n")
  # Chi-squared test
  digits = getOption("digits")
  x <- suppressWarnings(chisq.test(model$cont_table))
  cat("\nPearson's Chi-squared test:\n")
  out <- character()
  if (!is.null(x$statistic))
    out <- c(out, paste(names(x$statistic), "=", format(signif(x$statistic, max(1L, digits - 2L)))))
  if (!is.null(x$parameter))
    out <- c(out, paste(names(x$parameter), "=", format(signif(x$parameter, max(1L, digits - 2L)))))
  if (!is.null(x$p.value)) {
    fp <- format.pval(x$p.value, digits = max(1L, digits - 3L))
    out <- c(out, paste("p-value", if (substr(fp, 1L, 1L) == "<") fp else paste("=", fp)))
  }
  cat(strwrap(paste(out, collapse = ", ")), sep = "\n")
  cat("\n")
}

#' Print OneR models
#'
#' \code{print} method for class \code{OneR}.
#' @param x object of class \code{"OneR"}.
#' @param ... further arguments passed to or from other methods.
#' @details Prints the rules and the accuracy of an OneR model.
#' @author Holger von Jouanne-Diedrich, \email{r-project@ephorie.de}
#' @references \url{http://vonjd.github.io/OneR/}
#' @seealso \code{\link{OneR}}
#' @examples
#' model <- OneR(iris)
#' print(model)
#' @export
print.OneR <- function(x, ...) {
  model <- x
  longest <- max(nchar(names(model$rules)))
  cat("\nRules:\n")
  for (iter in 1:length(model$rules)) {
    len <- longest - nchar(names(model$rules[iter]))
    cat("If ", model$feature, " = ", names(model$rules[iter]), rep(" ", len)," then ", model$target, " = ", model$rules[[iter]], "\n", sep = "")
  }
  cat("\nAccuracy:\n")
  cat(model$correct_instances, " of ", model$total_instances, " instances classified correctly (", round(100 * model$correct_instances / model$total_instances, 2), "%)\n\n", sep = "")
}

#' Plot Diagnostics for an OneR object
#'
#' Plots a mosaic plot for the feature attribute and the target of the OneR model.
#' @param x object of class \code{"OneR"}.
#' @param ... further arguments passed to or from other methods.
#' @details If more than 20 levels are present for either the feature attribute or the target the function stops with an error.
#' @author Holger von Jouanne-Diedrich, \email{r-project@ephorie.de}
#' @references \url{http://vonjd.github.io/OneR/}
#' @seealso \code{\link{OneR}}
#' @keywords diagnostics
#' @examples
#' model <- OneR(iris)
#' plot(model)
#' @export
plot.OneR <- function(x, ...) {
  model <- x
  if (any(dim(model$cont_table) > 20)) stop("cannot plot: more than 20 levels")
  mosaicplot(t(model$cont_table), color = TRUE, main = "OneR model diagnostic plot")
}

#' Test OneR model objects
#'
#' Test if object is a OneR model.
#' @param x object to be tested.
#' @keywords OneR model
#' @author Holger von Jouanne-Diedrich, \email{r-project@ephorie.de}
#' @references \url{http://vonjd.github.io/OneR/}
#' @examples
#' model <- OneR(iris)
#' is.OneR(model) # evaluates to TRUE
#' @export
is.OneR <- function(x) inherits(x, "OneR")

#' Classification Evaluation function
#'
#' Function for evaluating a OneR classification model. Prints prediction vs. actual in absolute and relative numbers. Additionally it gives the accuracy and error rate.
#' @param prediction vector which contains the predicted values.
#' @param actual dataframe which contains the actual data. When there is more than one column the last last column is taken. A single vector is allowed too.
#' @details Invisibly returns a list with the number of correctly classified and total instances and a contingency table with the absolute numbers.
#' @author Holger von Jouanne-Diedrich, \email{r-project@ephorie.de}
#' @references \url{http://vonjd.github.io/OneR/}
#' @keywords evaluation accuracy
#' @examples
#' data <- iris
#' model <- OneR(data)
#' summary(model)
#' prediction <- predict(model, data)
#' eval_model(prediction, data)
#' @export
eval_model <- function(prediction, actual) {
  data <- actual
  if (is.list(data) == FALSE) data <- data.frame(data)
  if (typeof(as.vector(data[ , ncol(data)])) != typeof(prediction)) warning("data types of prediction and actual are different")
  cont <- table(prediction, data[ , ncol(data)])
  cont.m <- addmargins(cont)
  cat("\n", rep(" ", 11), "actual", sep = "")
  print(cont.m)
  cont.p <- prop.table(table(prediction, data[ , ncol(data)]))
  cont.pm <- round(addmargins(cont.p), 2)
  cat("\n", rep(" ", 11), "actual", sep = "")
  print(cont.pm)
  sum.cont.adj <- sum(cont[colnames(cont)[col(cont)] == rownames(cont)[row(cont)]])
  cat("\nAccuracy:\n", round(sum.cont.adj / sum(cont), 4), " (", sum.cont.adj, "/", sum(cont), ")", sep = "")
  cat("\n\nError rate:\n", round(1 - sum.cont.adj / sum(cont), 4), " (", sum(cont) - sum.cont.adj, "/", sum(cont), ")\n\n", sep = "")
  return(invisible(list(correct_instances = sum.cont.adj, total_instances = sum(cont), cont_table = cont)))
}
