library(bain)
library(lavaan)
example_data1 <- read.csv('example_data1.csv',header = T)
# Italians
Italians <- example_data1[1:313,]
# Colombians
Colombians <- example_data1[314:628,]
psych::describe(Italians)
psych::describe(Colombians)
# CFA model
model <- 'A  =~ Item1 + Item2 + Item3 + Item4 + Item5 + Item6 + Item7 + Item8 +
                Item9 + Item10 + Item11 + Item12 + Item13 + Item14 + Item16'
# computing the model
x <- cfa(model, data = example_data1, std.lv = TRUE, group = "COUNTRY")
# fit indexes
fit<- fitmeasures(x,c("chisq", "df",
                             "pvalue", "cfi", "rmsea","srmr","BIC"),output = "matrix")
fit
# estimated parameters
summary(x)

# testing metric invariance
set.seed(2020)
BF_metric <- bmi(x,type = "metric",tolerance = .2,minBF = 3,fraction = 2)

# testing Scalar invariance
set.seed(2020)
BF_scalar <- bmi(x,type = "scalar",tolerance = .2,minBF = 3,fraction = 2)
