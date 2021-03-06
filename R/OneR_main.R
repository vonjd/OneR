# OneR main function

#' One Rule function
#'
#' Builds a model according to the One Rule (OneR) machine learning classification algorithm.
#' @param x data frame with the last column containing the target variable.
#' @param formula formula, additionally the argument \code{data} is needed.
#' @param data data frame which contains the data, only needed when using the formula interface.
#' @param ties.method character string specifying how ties are treated, see 'Details'; can be abbreviated.
#' @param verbose if \code{TRUE} prints rank, names and predictive accuracy of the attributes in decreasing order (with \code{ties.method = "first"}).
#' @param ... arguments passed to or from other methods.
#' @return Returns an object of class "OneR". Internally this is a list consisting of the function call with the specified arguments, the names of the target and feature variables,
#' a list of the rules, the number of correctly classified and total instances and the contingency table of the best predictor vs. the target variable.
#' @keywords 1R OneR One Rule
#' @details All numerical data is automatically converted into five categorical bins of equal length. Instances with missing values are removed.
#' This is done by internally calling the default version of \code{\link{bin}} before starting the OneR algorithm.
#' To finetune this behaviour data preprocessing with the \code{\link{bin}} or \code{\link{optbin}} functions should be performed.
#' If data contains unused factor levels (e.g. due to subsetting) these are ignored and a warning is given.
#'
#' When there is more than one attribute with best performance either the first (from left to right) is being chosen (method \code{"first"}) or
#' the one with the lowest p-value of a chi-squared test (method \code{"chisq"}).
#' @author Holger von Jouanne-Diedrich
#' @references \url{https://github.com/vonjd/OneR}
#' @seealso \code{\link{bin}}, \code{\link{optbin}}, \code{\link{eval_model}}, \code{\link{maxlevels}}
#' @examples
#' data <- optbin(iris)
#' model <- OneR(data, verbose = TRUE)
#' summary(model)
#' plot(model)
#' prediction <- predict(model, data)
#' eval_model(prediction, data)
#'
#' ## The same with the formula interface:
#' data <- optbin(iris)
#' model <- OneR(Species ~., data = data, verbose = TRUE)
#' summary(model)
#' plot(model)
#' prediction <- predict(model, data)
#' eval_model(prediction, data)
#' @importFrom stats model.frame
#' @importFrom stats chisq.test
#' @export
OneR <- function(x, ...) UseMethod("OneR")

#' @export
OneR.default <- function(x, ...) {
  stop("data type not supported")
}

#' @export
#' @describeIn OneR method for formulas.
OneR.formula <- function(formula, data, ties.method = c("first", "chisq"), verbose = FALSE, ...) {
  call <- match.call()
  method <- match.arg(ties.method)
  mf <- model.frame(formula = formula, data = data, na.action = NULL)
  data <- mf[c(2:ncol(mf), 1)]
  OneR.data.frame(x = data, ties.method = ties.method, verbose = verbose, fcall = call)
}

#' @export
#' @describeIn OneR method for data frames.
OneR.data.frame <- function(x, ties.method = c("first", "chisq"), verbose = FALSE, ...) {
  if (!is.null(list(...)$fcall)) call <- list(...)$fcall
  else call <- match.call()
  method <- match.arg(ties.method)
  data <- x
  if (dim(data.frame(data))[2] < 2) stop("data must have at least two columns")
  data <- bin(data)
  if (nrow(data) == 0) stop("no data to analyse")
  # test if unused factor levels and drop them for analysis
  nlevels_orig <- sum(sapply(data, nlevels))
  data <- droplevels(data)
  nlevels_new <- sum(sapply(data, nlevels))
  if (nlevels_new < nlevels_orig) warning("data contains unused factor levels")
  # main routine for finding the best predictor(s)
  tables <- lapply(data[ , 1:(ncol(data)-1), drop = FALSE], table, data[ , ncol(data)])
  errors <- sapply(tables, nerrors)
  perf <- nrow(data) - errors
  target <- names(data[ , ncol(data), drop = FALSE])
  best <- which(perf == max(perf))
  # method "chisq
  if (length(best) > 1) {
    if (method == "chisq") {
      features <- names(data[ , best, drop = FALSE])
      p.values <- sapply(features, function(x) suppressWarnings(chisq.test(table(c(data[target], data[x])))$p.value))
      p.values[is.na(p.values)] <- Inf
      if (all(p.values == Inf)) warning("chi-squared tests failed, first best attribute is chosen instead")
      best <- best[which.min(p.values)]
    } else best <- best[1]
  }
  # preparation and output of results
  groups <- split(data[ , ncol(data), drop = FALSE], data[ , best])
  majority <- lapply(groups, mode)
  feature <- names(data[ , best, drop = FALSE])
  cont_table <- table(c(data[target], data[feature]))
  output <- c(call = call,
              target = target,
              feature = feature,
              rules = list(majority),
              correct_instances = max(perf),
              total_instances = nrow(data),
              cont_table = list(cont_table))
  class(output) <- "OneR"
  # print additional diagnostic information if wanted
  if (verbose == TRUE) {
    newbest <- which(which(perf == max(perf)) == best)
    accs <- round(100 * sort(perf, decreasing = TRUE) / nrow(data), 2)
    attr <- colnames(data[order(perf, decreasing = TRUE)])
    M <- matrix(c(as.character(attr), paste0(accs, "%")), ncol = 2)
    rownames(M) <- rank((100 - sort(perf, decreasing = TRUE)), ties.method = "min")
    rownames(M)[newbest] <- paste0(rownames(M)[newbest], " *")
    colnames(M) <- c("Attribute", "Accuracy")
    cat("\n")
    print(M, quote = FALSE)
    cat("---\nChosen attribute due to accuracy\nand ties method (if applicable): '*'\n\n")
  }
  output
}
