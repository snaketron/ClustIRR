# Benchmarking and profiling package comparison

pacman::p_load(ggplot2, microbenchmark, bench, profvis)

###########################################
## microbenchmark
##
## useful to benchmark single functions
## multiple function execution possible
##
##########################################

celsius_to_kelvin <- function(temp_C) {
  temp_K <- temp_C + 273.15
  return(temp_K)
}

fahrenheit_to_celsius <- function(temp_F) {
  temp_C <- (temp_F - 32) * 5 / 9
  return(temp_C)
}

# 1000 calls
result_microbenchmark <- microbenchmark(celsius_to_kelvin(-100:100), 
                                        fahrenheit_to_celsius(-100:100), 
                                        times=1000L)

print(result_microbenchmark)
ggplot2::autoplot(result_microbenchmark)



###########################################
## bench

## useful to run benchmarks against parameter grids 
## -> but official example code not working, maybe avoid
## can only compare similar functions, not different ones
##
##########################################

## NOT WORKING
# create_df <- function(rows, cols) {
#   as.data.frame(setNames(
#     replicate(cols, runif(rows, 1, 1000), simplify = FALSE),
#     rep_len(c("x", letters), cols)))
# }
# 
# # Run 4 data sizes across 3 samples with 2 replicates (24 total benchmarks)
# press(
#   rows = c(1000, 10000),
#   cols = c(10, 100),
#   rep = 1:2,
#   {
#     dat <- create_df(rows, cols)
#     bench::mark(
#       min_time = .05,
#       bracket = dat[dat$x > 500, ],
#       which = dat[which(dat$x > 500), ],
#       subset = subset(dat, x > 500)
#     )
#   }
# )

# compare different computing approaches
set.seed(42)
dat <- data.frame(x = runif(10000, 1, 1000), y=runif(10000, 1, 1000))
results_benchmark <- bench::mark(
  dat[dat$x > 500, ],
  dat[which(dat$x > 500), ],
  subset(dat, x > 500))

results_benchmark
autoplot(results_benchmark)

# won't work for different packages:
## bench::mark(celsius_to_kelvin(-100:100), fahrenheit_to_celsius(-100:100))
## Fehler: Each result must equal the first result:
## `celsius_to_kelvin(-100:100)` does not equal 
## `fahrenheit_to_celsius(-100:100)`


###########################################
## profile / profvis
##
## useful to compare larger code parts
## internal function, mark lines and use "Profile->Profile Selected Line(s)"
## profiles can be saved by wrapping code into function
## can't compare function faster than 5ms
##
##########################################

# Simulate random numbers
n <- 2e4
set.seed(2020)
ygroup <- data.frame(y = rnorm(mean = 100, sd = 20, n = n),
                     group = sample(LETTERS[1:16], size = n, replace = TRUE))

simulation <- data.frame(replicate(100, rnorm(mean = 50, sd = 10, n = n) +
                                     0.1 * ygroup$y))
simulation <- cbind(ygroup, simulation)
rm(ygroup)


# find x-variable that is correlating the most with y
index <- which.max(cor(simulation$y, simulation[, -(1:2)])) + 2
names(simulation[index])
cor(simulation$y, simulation[index])
simulation$X <- simulation[, index]


ggplot(simulation, aes(x = X, y = y)) +
  geom_jitter(alpha = 0.6, size = 1.2) +
  geom_smooth(method = "loess") +
  facet_wrap(~ group, nrow = 4) +
  labs(title = "loess, span = default = 0.75")


ggplot(simulation, aes(x = X, y = y)) +
  geom_jitter(alpha = 0.6, size = 1.2) +
  geom_smooth(method = "loess", span = 0.95) +
  facet_wrap(~ group, nrow = 4) +
  labs(title = "loess, span = default = 0.95")


ggplot(simulation, aes(x = X, y = y)) +
  geom_jitter(alpha = 0.6, size = 1.2) +
  geom_smooth(method = "loess", span = 0.1) +
  facet_wrap(~ group, nrow = 4) +
  labs(title = "loess, span = default = 0.1")


ggplot(simulation, aes(x = X, y = y)) +
  geom_jitter(alpha = 0.6, size = 1.2) +
  geom_smooth(method = "auto") +
  facet_wrap(~ group, nrow = 4) +
  labs(title = "auto, uses gam = generalized additive model")

ggplot(simulation, aes(x = X, y = y)) +
  geom_jitter(alpha = 0.6, size = 1.2) +
  geom_smooth(method = "lm") +
  facet_wrap(~ group, nrow = 4) +
  labs(title = "lm")


## save profile

profile_profvis = profvis({
# Simulate random numbers
n <- 2e4
set.seed(2020)
ygroup <- data.frame(y = rnorm(mean = 100, sd = 20, n = n),
                     group = sample(LETTERS[1:16], size = n, replace = TRUE))

simulation <- data.frame(replicate(100, rnorm(mean = 50, sd = 10, n = n) +
                                     0.1 * ygroup$y))
simulation <- cbind(ygroup, simulation)
rm(ygroup)


# find x-variable that is correlating the most with y
index <- which.max(cor(simulation$y, simulation[, -(1:2)])) + 2
names(simulation[index])
cor(simulation$y, simulation[index])
simulation$X <- simulation[, index]


ggplot(simulation, aes(x = X, y = y)) +
  geom_jitter(alpha = 0.6, size = 1.2) +
  geom_smooth(method = "loess") +
  facet_wrap(~ group, nrow = 4) +
  labs(title = "loess, span = default = 0.75")


ggplot(simulation, aes(x = X, y = y)) +
  geom_jitter(alpha = 0.6, size = 1.2) +
  geom_smooth(method = "loess", span = 0.95) +
  facet_wrap(~ group, nrow = 4) +
  labs(title = "loess, span = default = 0.95")


ggplot(simulation, aes(x = X, y = y)) +
  geom_jitter(alpha = 0.6, size = 1.2) +
  geom_smooth(method = "loess", span = 0.1) +
  facet_wrap(~ group, nrow = 4) +
  labs(title = "loess, span = default = 0.1")


ggplot(simulation, aes(x = X, y = y)) +
  geom_jitter(alpha = 0.6, size = 1.2) +
  geom_smooth(method = "auto") +
  facet_wrap(~ group, nrow = 4) +
  labs(title = "auto, uses gam = generalized additive model")

ggplot(simulation, aes(x = X, y = y)) +
  geom_jitter(alpha = 0.6, size = 1.2) +
  geom_smooth(method = "lm") +
  facet_wrap(~ group, nrow = 4) +
  labs(title = "lm")
}, interval = 0.005)

profile_profvis
