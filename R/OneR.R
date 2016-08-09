# OneR helper functions

#' Binning function
#'
#' Discretizes all numerical data in a dataframe into categorical bins of equal length or content or based on automatically determined clusters.
#' @param data dataframe which contains the data.
#' @param nbins number of bins (= levels).
#' @param labels character vector of labels for the resulting category.
#' @param method character string specifying the binning method, see 'Details'; can be abbreviated.
#' @param na.omit logical value whether instances with missing values should be removed.
#' @return A dataframe or vector.
#' @keywords binning discretization discretize clusters Jenks breaks
#' @details Character strings and logical strings are coerced into factors. Matrices are coerced into dataframes. When called with a single vector only the respective factor (and not a dataframe) is returned.
#' Method \code{"length"} gives intervals of equal length, method \code{"content"} gives intervals of equal content (via quantiles).
#' Method \code{"clusters"} determins \code{"nbins"} clusters via 1D kmeans with deterministic seeding of the initial cluster centres (Jenks natural breaks optimization).
#'
#' When \code{"na.omit = FALSE"} an additional level \code{"NA"} is added to each factor with missing values.
#' @author Holger von Jouanne-Diedrich
#' @references \url{https://github.com/vonjd/OneR}
#' @seealso \code{\link{OneR}}, \code{\link{optbin}}
#' @examples
#' data <- iris
#' str(data)
#' str(bin(data))
#' str(bin(data, nbins = 3))
#' str(bin(data, nbins = 3, labels = c("small", "medium", "large")))
#'
#' ## Difference between methods "length" and "content"
#' set.seed(1); table(bin(rnorm(900), nbins = 3))
#' set.seed(1); table(bin(rnorm(900), nbins = 3, method = "content"))
#'
#' ## Method "clusters"
#' intervals <- paste(levels(bin(faithful$waiting, nbins = 2, method = "cluster")), collapse = " ")
#' hist(faithful$waiting, main = paste("Intervals:", intervals))
#' abline(v = c(42.9, 67.5, 96.1), col = "blue")
#'
#' ## Missing values
#' bin(c(1:10, NA), nbins = 2, na.omit = FALSE) # adds new level "NA"
#' bin(c(1:10, NA), nbins = 2)                  # omits missing values by default (with warning)
#' @importFrom stats quantile
#' @importFrom stats kmeans
#' @export
bin <- function(data, nbins = 5, labels = NULL, method = c("length", "content", "clusters"), na.omit = TRUE) {
  method <- match.arg(method)
  vec <- FALSE
  if (is.atomic(data) == TRUE & is.null(dim(data)) == TRUE) { vec <- TRUE; data <- data.frame(data) }
  # could be a matrix -> dataframe (even with only one column)
  if (is.list(data) == FALSE) data <- data.frame(data)
  if (na.omit == TRUE) {
    len_rows_orig <- nrow(data)
    data <- na.omit(data)
    len_rows_new <- nrow(data)
    no_removed <- len_rows_orig - len_rows_new
    if (no_removed > 0) warning(paste(no_removed, "instance(s) removed due to missing values"))
  }
  if (!is.null(labels)) if (nbins != length(labels)) stop("number of 'nbins' and 'labels' differ")
  if (nbins <= 1) stop("nbins must be bigger than 1")
  data[] <- lapply(data, function(x) if (is.numeric(x)) {
    if (length(unique(x)) <= nbins) as.factor(x)
    else {
      if (method == "content") nbins <- add_range(x, na.omit(quantile(x, (1:(nbins-1)/nbins), na.rm = TRUE)))
      if (method == "clusters") {
        midpoints <- sort(kmeans(na.omit(x), centers = seq(min(x, na.rm = TRUE), max(x, na.rm = TRUE), length = nbins))$centers)
        nbins <- add_range(x, na.omit(filter(midpoints, c(1/2, 1/2))))
      }
      CUT(x, breaks = unique(nbins), labels = labels)
    }
  } else as.factor(x))
  data[] <- lapply(data, function(x) if(any(is.na(as.character(x)))) ADDNA(x) else x)
  if (vec) { data <- unlist(data); names(data) <- NULL }
  data
}

#' Optimal Binning function
#'
#' Discretizes all numerical data in a dataframe into categorical bins where the cut points are optimally aligned with the target categories, thereby a factor is returned.
#' When building a OneR model this could result in fewer rules with enhanced accuracy.
#' @param data dataframe which contains the data. When \code{formula = NULL} (the default) the last column must be the target variable.
#' @param formula formula interface for the \code{optbin} function.
#' @param method character string specifying the method for optimal binning, see 'Details'; can be abbreviated.
#' @param na.omit logical value whether instances with missing values should be removed.
#' @return A dataframe with the target variable being in the last column.
#' @keywords binning discretization discretize
#' @details The cutpoints are calculated by pairwise logistic regressions (method \code{"logreg"}), information gain (method \code{"infogain"}) or as the means of the expected values of the respective classes (\code{"naive"}).
#' The function is likely to give unsatisfactory results when the distributions of the respective classes are not (linearly) separable. Method \code{"naive"} should only be used when distributions are (approximately) normal,
#' although in this case \code{"logreg"} should give comparable results, so it is the preferable (and therefore default) method.
#'
#' Method \code{"infogain"} is an entropy based method which calculates cut points based on information gain. The idea is that uncertainty is minimized by making the resulting bins as pure as possible. This method is the standard method of many decision tree algorithms.
#'
#' Character strings and logical strings are coerced into factors. Matrices are coerced into dataframes. If the target is numeric it is turned into a factor with the number of levels equal to the number of values. Additionally a warning is given.
#'
#' When \code{"na.omit = FALSE"} an additional level \code{"NA"} is added to each factor with missing values.
#' If the target contains unused factor levels (e.g. due to subsetting) these are ignored and a warning is given.
#' @author Holger von Jouanne-Diedrich
#' @references \url{https://github.com/vonjd/OneR}
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
optbin <- function(data, formula = NULL, method = c("logreg", "infogain", "naive"), na.omit = TRUE) {
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
    len_rows_orig <- nrow(data)
    data <- na.omit(data)
    len_rows_new <- nrow(data)
    no_removed <- len_rows_orig - len_rows_new
    if (no_removed > 0) warning(paste(no_removed, "instance(s) removed due to missing values"))
  } else {
    # only add NA to target
    if(any(is.na(unlist(data[ncol(data)])))) data[ncol(data)] <- ADDNA(unlist(data[ncol(data)]))
  }
  target <- data[ , ncol(data)]
  # Test if unused factor levels and drop them for analysis
  nlevels_orig <- nlevels(target)
  target <- droplevels(target)
  nlevels_new <- nlevels(target)
  if (nlevels_new < nlevels_orig) warning("target contains unused factor levels")
  nbins <- length(unique(target))
  if (nbins <= 1) stop("number of target levels must be bigger than 1")
  data[] <- lapply(data, function(x) if (is.numeric(x)) {
    if (length(unique(x)) <= nbins) as.factor(x) else optcut(x, target, method)
  } else as.factor(x))
  data[] <- lapply(data, function(x) if(any(is.na(as.character(x)))) ADDNA(x) else x)
  data
}

#' Remove factors with too many levels
#'
#' Removes all columns of a dataframe where a factor (or character string) has more than a maximum number of levels.
#' @param data dataframe which contains the data.
#' @param maxlevels number of maximum factor levels.
#' @param na.omit logical value whether missing values should be treated as a level, defaults to omit missing values before counting.
#' @return A dataframe.
#' @details Often categories that have very many levels are not useful in modelling OneR rules because they result in too many rules and tend to overfit.
#' Examples are IDs or names.
#'
#' Character strings are treated as factors although they keep their datatype. Numeric data is left untouched.
#' If data contains unused factor levels (e.g. due to subsetting) these are ignored and a warning is given.
#' @author Holger von Jouanne-Diedrich
#' @references \url{https://github.com/vonjd/OneR}
#' @seealso \code{\link{OneR}}
#' @examples
#' df <- data.frame(numeric = c(1:26), alphabet = letters)
#' str(df)
#' str(maxlevels(df))
#' @export
maxlevels <- function(data, maxlevels = 20, na.omit = TRUE) {
  if (is.list(data) == FALSE) stop("data must be a dataframe")
  if (maxlevels <= 2) stop("maxlevels must be bigger than 2")
  tmp <- suppressWarnings(bin(data, nbins = 2, na.omit = na.omit))
  # Test if unused factor levels and drop them for analysis
  nlevels_orig <- sapply(tmp, nlevels)
  tmp <- droplevels(tmp)
  nlevels_new <- sapply(tmp, nlevels)
  if (sum(nlevels_new) < sum(nlevels_orig)) warning("data contains unused factor levels")
  cols <- nlevels_new <= maxlevels
  data[cols]
}

#' Predict method for OneR models
#'
#' Predict cases or probabilities based on OneR model object.
#' @param object object of class \code{"OneR"}.
#' @param newdata dataframe in which to look for the feature variable with which to predict.
#' @param type character string denoting the type of predicted value returned. Default \code{"class"} gives a named vector with the predicted classes, \code{"prob"} gives a matrix whose columns are the probability of the first, second, etc. class.
#' @param ... further arguments passed to or from other methods.
#' @return The default is a named vector with the predicted classes, if \code{"type = prob"} a matrix is returned whose columns are the probability of the first, second, etc. class.
#' @details \code{newdata} can have the same format as used for building the model but must at least have the feature variable that is used in the OneR rules.
#' If cases appear that were not present when building the model the predicted case is \code{UNSEEN} or \code{NA} when \code{"type = prob"}.
#' @author Holger von Jouanne-Diedrich
#' @references \url{https://github.com/vonjd/OneR}
#' @seealso \code{\link{OneR}}
#' @examples
#' model <- OneR(iris)
#' prediction <- predict(model, iris[1:4])
#' eval_model(prediction, iris[5])
#'
#' ## type prob
#' predict(model, data.frame(Petal.Width = seq(0.5, 2.5, 0.5)), type = "prob")
#' @export
predict.OneR <- function(object, newdata, type = c("class", "prob"), ...) {
  type <- match.arg(type)
  if (is.list(newdata) == FALSE) stop("newdata must be a dataframe")
  if (all(names(newdata) != object$feature)) stop("cannot find feature column in newdata")
  model <- object
  data <- newdata
  index <- which(names(data) == model$feature)[1]
  if (is.numeric(data[ , index])) {
    levels <- names(model$rules)
    if (substring(levels[1], 1, 1) == "(" & grepl(",", levels[1]) == TRUE & substring(levels[1], nchar(levels[1]), nchar(levels[1])) == "]") {
      features <- as.character(cut(unlist(data[ , index]), breaks = get_breaks(levels)))
    } else features <- as.character(data[ , index])
  } else features <- as.character(data[ , index])
  features[is.na(features)] <- "NA"
  if (type == "prob") {
    probs <- prop.table(model$cont_table, margin = 2)
    probrules <- lapply(names(model$rules), function(x) probs[ , x])
    names(probrules) <- names(model$rules)
    return(t(sapply(features, function(x) if (is.null(probrules[[x]]) == TRUE) rep(NA, dim(model$cont_table)[1]) else probrules[[x]])))
  }
  sapply(features, function(x) if (is.null(model$rules[[x]]) == TRUE) "UNSEEN" else model$rules[[x]])
}

#' Summarize OneR models
#'
#' \code{summary} method for class \code{OneR}.
#' @param object object of class \code{"OneR"}.
#' @param ... further arguments passed to or from other methods.
#' @details Prints the rules of the OneR model, the accuracy, a contingency table of the feature attribute and the target and performs a chi-squared test on this table.
#'
#' In the contingency table the maximum values in each column are highlighted by adding a '*', thereby representing the rules of the OneR model.
#' @author Holger von Jouanne-Diedrich
#' @references \url{https://github.com/vonjd/OneR}
#' @seealso \code{\link{OneR}}
#' @keywords diagnostics
#' @examples
#' model <- OneR(iris)
#' summary(model)
#' @importFrom stats addmargins
#' @importFrom stats chisq.test
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
#' @author Holger von Jouanne-Diedrich
#' @references \url{https://github.com/vonjd/OneR}
#' @seealso \code{\link{OneR}}
#' @examples
#' model <- OneR(iris)
#' print(model)
#' @export
print.OneR <- function(x, ...) {
  model <- x
  cat("\nCall:\n")
  print(model$call)
  cat("\nRules:\n")
  longest <- max(nchar(names(model$rules)))
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
#' @author Holger von Jouanne-Diedrich
#' @references \url{https://github.com/vonjd/OneR}
#' @seealso \code{\link{OneR}}
#' @keywords diagnostics
#' @examples
#' model <- OneR(iris)
#' plot(model)
#' @importFrom graphics mosaicplot
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
#' @return a logical whether object is of class "OneR".
#' @keywords OneR model
#' @author Holger von Jouanne-Diedrich
#' @references \url{https://github.com/vonjd/OneR}
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
#' @return Invisibly returns a list with the number of correctly classified and total instances and a confusion matrix with the absolute numbers.
#' @author Holger von Jouanne-Diedrich
#' @references \url{https://github.com/vonjd/OneR}
#' @keywords evaluation accuracy
#' @examples
#' data <- iris
#' model <- OneR(data)
#' summary(model)
#' prediction <- predict(model, data)
#' eval_model(prediction, data)
#' @importFrom stats addmargins
#' @export
eval_model <- function(prediction, actual) {
  if (is.list(actual) == FALSE) actual <- data.frame(actual)
  actual <- actual[ , ncol(actual)]
  if (any(is.na(actual))) warning("'actual' contains missing values, results may be false")
  if (typeof(as.vector(actual)) != typeof(prediction)) warning("data types of prediction and actual are different")
  conf <- table(prediction, actual)
  conf.m <- addmargins(conf)
  cat("\nConfusion matrix (absolute):\n")
  print(conf.m)
  conf.p <- prop.table(conf)
  conf.pm <- round(addmargins(conf.p), 2)
  cat("\nConfusion matrix (relative):\n")
  print(conf.pm)
  sum.conf.adj <- sum(conf[colnames(conf)[col(conf)] == rownames(conf)[row(conf)]])
  cat("\nAccuracy:\n", round(sum.conf.adj / sum(conf), 4), " (", sum.conf.adj, "/", sum(conf), ")", sep = "")
  cat("\n\nError rate:\n", round(1 - sum.conf.adj / sum(conf), 4), " (", sum(conf) - sum.conf.adj, "/", sum(conf), ")\n\n", sep = "")
  invisible(list(correct_instances = sum.conf.adj, total_instances = sum(conf), conf_matrix = conf))
}
