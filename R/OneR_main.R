# OneR main function

#' One Rule function
#'
#' Builds a model according to the One Rule (OneR) machine learning classification algorithm.
#' @param data dataframe, which contains the data. When \code{formula = NULL} (the default) the last column must be the target variable.
#' @param formula formula interface for the \code{OneR} function.
#' @param ties.method a character string specifying how ties are treated, see 'Details'; can be abbreviated.
#' @param verbose If \code{TRUE} prints rank, names and predictive accuracy of the attributes in decreasing order (with \code{ties.method = "first"}).
#' @return Returns an object of class "OneR". Internally this is a list consisting of the names of the target and feature variables, a list of the rules,
#' the number of correctly classified and total instances and the contingency table of the best predictor vs. the target variable.
#' @keywords 1R OneR One Rule
#' @details All numerical data is automatically converted into five categorical bins of equal length. Instances with missing values are removed.
#' This is done by internally calling the default version of \code{\link{bin}} before starting the OneR algorithm.
#' To finetune this behaviour data preprocessing with the \code{\link{bin}} or \code{\link{optbin}} functions should be performed.
#'
#' When there is more than one attribute with best performance either the first (from left to right) is being chosen (method \code{"first"}) or
#' the one with the lowest p-value of a chi-squared test (method \code{"chisq"}).
#' @author Holger von Jouanne-Diedrich
#' @references \url{http://vonjd.github.io/OneR/}
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
#' model <- OneR(formula = Species ~., data = data, verbose = TRUE)
#' summary(model)
#' plot(model)
#' prediction <- predict(model, data)
#' eval_model(prediction, data)
#' @importFrom stats model.frame
#' @importFrom stats chisq.test
#' @export
OneR <- function(data, formula = NULL, ties.method = c("first", "chisq"), verbose = FALSE) {
  method <- match.arg(ties.method)
  if (class(formula) == "formula") {
    mf <- model.frame(formula = formula, data = data, na.action = NULL)
    data <- mf[c(2:ncol(mf), 1)]
  } else if (is.null(formula) == FALSE) stop("invalid formula")
  if (dim(data.frame(data))[2] < 2) stop("data must have at least two columns")
  data <- bin(data)
  if (nrow(data) == 0) stop("no data to analyse")
  perf <- c()
  for (iter in 1:(ncol(data) - 1)) {
    groups <- split(data, data[ , iter])
    majority <- lapply(groups, mode)
    real <- lapply(groups, function(x) x[ , ncol(x)])
    perf <- c(perf, sum(unlist(Map("==", real, majority))))
  }
  target <- names(data[ , ncol(data), drop = FALSE])
  best <- which(perf == max(perf))
  if (length(best) > 1) {
    if (method == "chisq") {
      features <- names(data[ , best, drop = FALSE])
      p.values <- sapply(features, function(x) suppressWarnings(chisq.test(table(c(data[target], data[x])))$p.value))
      p.values[is.na(p.values)] <- Inf
      if (all(p.values == Inf)) warning("chi-squared tests failed, first best attribute is chosen instead")
      best <- best[which.min(p.values)]
    } else best <- best[1]
  }
  groups <- split(data, data[ , best])
  majority <- lapply(groups, mode)
  feature <- names(data[ , best, drop = FALSE])
  cont_table <- table(c(data[target], data[feature]))
  output <- c(target = target, feature = feature, rules = list(majority), correct_instances = max(perf), total_instances = nrow(data), cont_table = list(cont_table))
  class(output) <- "OneR"
  if (verbose == TRUE) {
    newbest <- which(which(perf == max(perf)) == best)
    accs <- round(100 * sort(perf, decreasing = TRUE) / nrow(data), 2)
    attr <- colnames(data[order(perf, decreasing = TRUE)])
    M <- matrix(c(as.character(attr), paste(accs, "%", sep = "")), ncol = 2)
    rownames(M) <- rank((100 - sort(perf, decreasing = TRUE)), ties.method = "min")
    rownames(M)[newbest] <- paste(rownames(M)[newbest], " *", sep = "")
    colnames(M) <- c("Attribute", "Accuracy")
    cat("\n")
    print(M, quote = FALSE)
    cat("---\nChosen attribute due to accuracy\nand ties method (if applicable): '*'\n\n")
  }
  return(output)
}
