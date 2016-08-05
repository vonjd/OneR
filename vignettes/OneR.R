## ------------------------------------------------------------------------
library(OneR)

## ------------------------------------------------------------------------
data <- optbin(iris)

## ------------------------------------------------------------------------
model <- OneR(data, verbose = TRUE)

## ------------------------------------------------------------------------
summary(model)

## ---- fig.width=7.15, fig.height=5---------------------------------------
plot(model)

## ------------------------------------------------------------------------
prediction <- predict(model, data)

## ------------------------------------------------------------------------
eval_model(prediction, data)

## ------------------------------------------------------------------------
data(breastcancer)
data <- breastcancer

## ------------------------------------------------------------------------
set.seed(12) # for reproducibility
random <- sample(1:nrow(data), 0.8 * nrow(data))
data_train <- optbin(data[random, ])
data_test <- data[-random, ]

## ------------------------------------------------------------------------
model_train <- OneR(data_train, verbose = TRUE)

## ------------------------------------------------------------------------
summary(model_train)

## ---- fig.width=7.15, fig.height=5---------------------------------------
plot(model_train)

## ------------------------------------------------------------------------
prediction <- predict(model_train, data_test)

## ------------------------------------------------------------------------
eval_model(prediction, data_test)

## ------------------------------------------------------------------------
data <- iris
str(data)
str(bin(data))
str(bin(data, nbins = 3))
str(bin(data, nbins = 3, labels = c("small", "medium", "large")))

## ------------------------------------------------------------------------
set.seed(1); table(bin(rnorm(900), nbins = 3))
set.seed(1); table(bin(rnorm(900), nbins = 3, method = "content"))

## ---- fig.width=7.15, fig.height=5---------------------------------------
intervals <- paste(levels(bin(faithful$waiting, nbins = 2, method = "cluster")), collapse = " ")
hist(faithful$waiting, main = paste("Intervals:", intervals))
abline(v = c(42.9, 67.5, 96.1), col = "blue")

## ------------------------------------------------------------------------
bin(c(1:10, NA), nbins = 2, na.omit = FALSE) # adds new level "NA"
bin(c(1:10, NA), nbins = 2)

## ------------------------------------------------------------------------
df <- data.frame(numeric = c(1:26), alphabet = letters)
str(df)
str(maxlevels(df))

## ------------------------------------------------------------------------
model <- OneR(iris)
predict(model, data.frame(Petal.Width = seq(0.5, 2.5, 0.5)), type = "prob")

## ---- eval=FALSE---------------------------------------------------------
#  help(package = OneR)

