#' ---
#' title: "Estimating false positive rates with simulation"
#' author: Joel Levin
#' date: "April 12, 2020"
#' output: github_document
#' ---

#' This version: `r format(Sys.time())`.  
#' You can contact Joel at joelmlevin@gmail.com.

#+ echo=FALSE, include=FALSE
library(tidyverse)
options(scipen = 999)

#' ## Description

#' This is a script to test the true false positive rate of any statistical test that one is using to evaluate experimental results. The motivation to conduct this test came from [this working paper](https://psyarxiv.com/cyv6d/) by Ryan, Evers, & Moore (2018). I used this stimulation to estimate the false positive rate of a poisson regression, but it can be used to evaluate any test.

#+
#' ## Basic functions

#' This function generates a single vector containing p values of a test with random experimental conditions. If a statistical test is behaving properly, p values of a given magnitude (e.g., .05) should appear with a corresponding frequency (e.g., 5% of the time). The arguments include:
#+
#'1. `outcome_vector`: the real vector of dependent measures from your study. this should be a numerical vector.
#'2. `num_conditions`: the number of experimental conditions in your study. 
#'3. `test_type`: the statistical test used. this currently supports `"poisson", "ols", "t-test"`.
#'4. `replications`: the number of tests to simulate. this defaults to 10,000 tests, which is usually sufficient. if you're not in a hurry and have a modern-ish machine, make it 100,000.

simulate_fp <- function(outcome_vector, num_conditions, test_type = c("poisson", "ols", "t-test"), replications = 10000) {
  
  # tests
  if(test_type != "poisson" & test_type != "ols" & test_type != "t-test") { # test type must conform to available types
    stop("Please specify a test type using the `test_type` argument. The available tests are 'poisson', 'ols', and 't-test.")
    }
  if(is.numeric(outcome_vector != TRUE)) { # outcome vector must be numeric
    stop("outcome_vector must be a numeric vector.")
    }
  if(all.equal(num_conditions, as.integer(num_conditions)) != TRUE | num_conditions < 1) { # must be an integer greater than 1
    stop("Please specify the number of conditions in your study using the num_conditions argument.")
    }
  
  tests <- rep(NA, replications) # generate an empty vector of appropriate length. this will be overwritten with p values
  
  # looping to fill in each element of the `tests` vector
  for(n in 1:replications) {
    false_condition <- sample(c(1, num_conditions), length(outcome_vector), replace = TRUE) # generating a vector of random condition dummies

    # conducting the tests    
    if(test_type == "poisson") {
    tests[n] <-  coef(summary(glm(outcome_vector ~ false_condition, family = poisson)))[2, 1:4][4]
    } 
    if(test_type == "ols") {
      tests[n] <-  coef(summary(lm(outcome_vector ~ false_condition)))[2, 1:4][4]
    }
    if(test_type == "t-test") {
      tests[n] <-  t.test(outcome_vector ~ false_condition)$p.value
    }
  }
  return(tests)
}
    
#' This function takes the output of the previous function and returns a diagnostic table or plot, which indicates whether p values appear as frequently as expected for a well-behaving test. The arguments include:
#+
#'1. `data`: this should be a vector of simulated p values, such as those provided using the above function.
#'1. `type`: this indicates whether to output a table of values or a plot. see below for examples.
#'1. `quantiles`: a vector of percentiles at which to compare the simulated p values to their observed frequency. this argument is prespecified at recommended values and is therefore optional.

diagnostics <- function(data, type = c("table", "plot"), quantiles = c(.01, .05, .10, .5)) {
  values <- round(quantile(data, quantiles), 4)
  difference <- abs(quantiles - values)
  temp <- cbind(quantiles, values, difference)
  
    if(type == "table") {
      return(temp)
      }

    if(type == "plot") {
    plot(temp)
    abline(0, 1)
    }
}

#' ## Using the functions

#' Generating a "real" set of experimental conditions. In practice, this would be the actual vector of experimental conditions in your study.
real_conditions <- sample(c(1, 2), 100, replace = TRUE)

#' Generating a "real" set of dependent measures (continuous). In practice, this would be the actual vector of dependent measures in your study.
real_dv <- rnorm(100, 20, 10)

#' Conducting a t test on the "real" data. 
t.test(real_dv ~ real_conditions)

#' Now using the function to simulate p values for random experimental conditions
simulated_ps <- simulate_fp(outcome_vector = real_dv, num_conditions = 2, test_type = "t-test", replications = 10000)

simulated_ps[1:10]

#' ### Now using the diagnostic functions.
#' 
#' The diagnostic table gives shows you the p values (values) at various percentiles (quantiles). The closer the two are, the better the test is behaving. Note that this is also affected by the number of replications used to generate the simulated data. 
diagnostics(simulated_ps, type = "table")

#' The diagnostic plot simply combines both values in a plot with a reference line. The closer the points are to the reference line, the better behaved the test.
diagnostics(simulated_ps, type = "plot")

#' We can also conduct a test to more rigorously evaluate whether our statistical test is well behaved. The below function was written by [Eric Archer](mailto:eric.archer@noaa.gov) for the `swfscMisc` package. It is included in the text of this script to minimize dependencies.

uniform.test <- function(hist.output, B = NULL) {
  break.diff <- diff(hist.output$breaks)
  probs <- break.diff / sum(break.diff)
  if (is.null(B)) {
    chisq.test(x = hist.output$counts, p = probs)
  } else {
    chisq.test(x = hist.output$counts, p = probs, simulate.p.value = TRUE, B = B)
  }
}

#' This function returns a chi squared statistic that, if significant, tells us that our data are *not* uniform, meaning that the p values cannot be interpreted conventionally. The large p-value here tells us that we're fine.
uniform.test(hist(simulated_ps))

#' (Later, add in the plotting functions)