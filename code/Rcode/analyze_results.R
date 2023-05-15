library(tidyverse)

results <- tibble()

for (n in c(3,4,5)){
  for (d in c(2,3,4,5,6)){
    f1 <-  paste0("n_", n, "_d_", d, ".csv")
    f2 <- paste0("n_", n, "_d_", d, "_times.csv")
    if (file.exists(f1)){
      dt <- read_csv(f1)
      n <- dt$n[1]
      d <- dt$d[1]
      obs <- mean(dt$npos)
      expected <- 1/(2^n) * (d-1)^(n / 2)
      exppfeas <- 1 - (1 - 1/(2^n))^((d-1)^(n / 2))
      exppois <- 1- (exp(-2^(-n) * (d-1)^(n/2)))
      pfeas <- mean(dt$npos > 0)
      # now timing
      dt2 <- read_csv(f2)
      tot_time <- sum(dt2$time)
      results <- rbind(results, tibble(
        n = n, 
        d = d, 
        observed = obs,
        expected = expected,
        pfeas = pfeas,
        exppfeas = exppfeas,
        exppois = exppois,
        time_sec = tot_time
      ))
    }
  }
}

print(results)
par(mfcol = c(1,2))
plot(results$observed, results$expected, main = "# equil"); abline(c(0,1))
plot(results$pfeas, results$exppfeas, main = "p feas"); abline(c(0,1))
points(results$pfeas, results$exppois, col = "blue")

# dt <- read_csv("n_3_d_5.csv")
# mu <- mean(dt$nsol)
# vr <- mean(dt$nsol^2) - mu^2
# pp <- mu / vr
# rr <- mu^2 / (vr - mu)
# 
# hist(dt$nsol, prob = TRUE, ylim = c(0, .25), breaks = 15) # may need to tweak the y axis.
# lines(0:max(dt$nsol), dnbinom(0:max(dt$nsol), p = pp, size = rr), col = 'red')
# lines(0:max(dt$nsol), dpois(0:max(dt$nsol), mu), col = 'blue')
