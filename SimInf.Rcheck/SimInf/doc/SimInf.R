### R code from vignette source 'SimInf.Rnw'

###################################################
### code chunk number 1: SimInf.Rnw:106-107
###################################################
options(prompt = "R> ", continue = "+  ", width = 70, useFancyQuotes = FALSE)


###################################################
### code chunk number 2: install-SimInf (eval = FALSE)
###################################################
## install.packages("SimInf")


###################################################
### code chunk number 3: load-SimInf
###################################################
library("SimInf")


###################################################
### code chunk number 4: SIR-u0
###################################################
n <- 1000
u0 <- data.frame(S = rep(999, n), I = rep(1, n), R = rep(0, n))


###################################################
### code chunk number 5: SIR-tspan
###################################################
tspan <- seq(from = 1, to = 180, by = 7)


###################################################
### code chunk number 6: SIR-run
###################################################
model <- SIR(u0 = u0, tspan = tspan, beta = 0.16, gamma = 0.077)
set.seed(123)
set_num_threads(1)
result <- run(model = model)


###################################################
### code chunk number 7: SIR-show
###################################################
result


###################################################
### code chunk number 8: SIR-plot (eval = FALSE)
###################################################
## plot(result)
## plot(result, index = 1:10, range = FALSE)


###################################################
### code chunk number 9: SIR-I
###################################################
pdf("SimInf-SIR-I.pdf", width = 10, height = 5)
plot(result)
dev.off()


###################################################
### code chunk number 10: SIR-II
###################################################
pdf("SimInf-SIR-II.pdf", width = 10, height = 5)
plot(result, index = 1:10, range = FALSE)
dev.off()


###################################################
### code chunk number 11: SIR-head-trajectory
###################################################
head(trajectory(model = result, index = 1))


###################################################
### code chunk number 12: SIR-events
###################################################
events <- data.frame(
event = c("enter", "extTrans", "exit", "enter"),
time = c(2, 3, 4, 4), node = c(3, 1, 2, 1),
dest = c(0, 3, 0, 0), n = c(5, 7, 0, 1),
proportion = c(0, 0, 0.2, 0), select = c(1, 4, 4, 2),
shift = c(0, 0, 0, 0))


###################################################
### code chunk number 13: SIR-events-show
###################################################
events


###################################################
### code chunk number 14: SIR-u0-add
###################################################
u0 <- data.frame(S = rep(0, 5), I = rep(0, 5), R = rep(0, 5))
add <- data.frame(event = "enter", time = rep(1:10, each = 5),
  node = 1:5, dest = 0, n = 1:5, proportion = 0, select = 1, shift = 0)


###################################################
### code chunk number 15: SIR-infect
###################################################
infect <- data.frame(event = "enter", time = 25, node = 5,
  dest = 0, n = 1, proportion = 0, select = 2, shift = 0)


###################################################
### code chunk number 16: SIR-move
###################################################
move <- data.frame(event = "extTrans", time = 35:45, node = c(5, 5, 5,
  5, 4, 4, 4, 3, 3, 2, 1), dest = c(4, 3, 3, 1, 3, 2, 1, 2, 1, 1, 2),
  n = 5, proportion = 0, select = 4, shift = 0)


###################################################
### code chunk number 17: SIR-remove
###################################################
remove <- data.frame(event = "exit", time = c(70, 110),
  node = rep(1:5, each = 2), dest = 0, n = 0, proportion = 0.2,
  select = 4, shift = 0)


###################################################
### code chunk number 18: SimInf.Rnw:1361-1368 (eval = FALSE)
###################################################
## events <- rbind(add, infect, move, remove)
## model <- SIR(u0 = u0, tspan = 1:180, events = events, beta = 0.16,
##   gamma = 0.077)
## set.seed(3)
## set_num_threads(1)
## result <- run(model)
## plot(result, index = 1:5, range = FALSE)


###################################################
### code chunk number 19: SIR-events-I
###################################################
events <- rbind(add, infect, move, remove)
model <- SIR(u0 = u0, tspan = 1:180, events = events, beta = 0.16,
gamma = 0.077)
set.seed(3)
set_num_threads(1)
result <- run(model)
pdf("SimInf-SIR-events-I.pdf", width = 10, height = 5)
plot(result, index = 1:5, range = FALSE)
dev.off()


###################################################
### code chunk number 20: SIR-events-II
###################################################
model_no_infected <- SIR(u0 = u0, tspan = 1:180,
events = rbind(add, move, remove), beta = 0.16,
gamma = 0.077)
set.seed(3)
set_num_threads(1)
result_no_infected <- run(model_no_infected)
pdf("SimInf-SIR-events-II.pdf", width = 10, height = 5)
plot(result_no_infected, index = 1:5, range = FALSE)
dev.off()


###################################################
### code chunk number 21: SIR-replicate
###################################################
set.seed(123)
set_num_threads(1)
mean(replicate(n = 1000, {
  nI <- trajectory(run(model = model), index = 1:4)$I
  sum(nI) > 0
}))


###################################################
### code chunk number 22: load-SISe_sp-data
###################################################
data("nodes", package = "SimInf")
u0 <- u0_SISe()
events <- events_SISe()


###################################################
### code chunk number 23: distance-matrix
###################################################
d_ik <- distance_matrix(x = nodes$x, y = nodes$y, cutoff = 2500)


###################################################
### code chunk number 24: create-SISesp-u0
###################################################
set.seed(123)
i <- sample(x = 1:1600, size = 160)
u0$I[i] <- as.integer(u0$S[i] * 0.05)
u0$S[i] <- u0$S[i] - u0$I[i]


###################################################
### code chunk number 25: create-SISesp-model
###################################################
model <- SISe_sp(u0 = u0, tspan = 1:1460, events = events, phi = 0,
  upsilon = 0.012, gamma = 0.1, alpha = 1, beta_t1 = 0.10,
  beta_t2 = 0.12, beta_t3 = 0.12, beta_t4 = 0.10, end_t1 = 91,
  end_t2 = 182, end_t3 = 273, end_t4 = 365, distance = d_ik,
  coupling = 0.2)


###################################################
### code chunk number 26: SimInf.Rnw:1558-1567 (eval = FALSE)
###################################################
## plot(NULL, xlim = c(0, 1500), ylim = c(0, 0.18), ylab = "Prevalance",
##   xlab = "Time")
## set.seed(123)
## set_num_threads(1)
## replicate(5, {
##   result <- run(model = model)
##   p <- prevalence(model = result, formula = I ~ S + I, level = 2)
##   lines(p)
## })


###################################################
### code chunk number 27: SimInf.Rnw:1576-1582 (eval = FALSE)
###################################################
## gdata(model, "coupling") <- 0.1
## replicate(5, {
##   result <- run(model = model)
##   p <- prevalence(model = result, formula = I ~ S + I, level = 2)
##   lines(p, col = "blue", lty = 2)
## })


###################################################
### code chunk number 28: SIR-mparse-I
###################################################
transitions <- c("S -> b*S*I/(S+I+R) -> I", "I -> g*I -> R")
compartments <- c("S", "I", "R")


###################################################
### code chunk number 29: SIR-mparse-II
###################################################
n <- 1000
u0 <- data.frame(S = rep(99, n), I = rep(5, n), R = rep(0, n))
model <- mparse(transitions = transitions, compartments = compartments,
  gdata = c(b = 0.16, g = 0.077), u0 = u0, tspan = 1:180)


###################################################
### code chunk number 30: SIR-mparse-III (eval = FALSE)
###################################################
## set.seed(123)
## set_num_threads(1)
## result <- run(model = model)
## plot(result)


###################################################
### code chunk number 31: SIR-mparse-IV
###################################################
set.seed(123)
set_num_threads(1)
result <- run(model = model)
plot(result)


###################################################
### code chunk number 32: SIR-mparse-incidence (eval = FALSE)
###################################################
## transitions <- c("S -> b*S*I/(S+I+R) -> I + Icum", "I -> g*I -> R")
## compartments <- c("S", "I", "Icum", "R")


###################################################
### code chunk number 33: SIR-mparse-incidence-run (eval = FALSE)
###################################################
## n <- 1000
## u0 <- data.frame(S = rep(99, n), I = rep(1, n), Icum = rep(0, n),
##   R = rep(0, n))
## model <- mparse(transitions = transitions, compartments = compartments,
##   gdata = c(b = 0.16, g = 0.077), u0 = u0, tspan = 1:150)
## set.seed(123)
## set_num_threads(1)
## result <- run(model = model)


###################################################
### code chunk number 34: SIR-mparse-incidence-trajectory (eval = FALSE)
###################################################
## traj <- trajectory(model = result, compartments = "Icum")
## cases <- stepfun(result@tspan[-1], diff(c(0, traj$Icum[traj$node == 1])))
## avg_cases <- c(0, diff(by(traj, traj$time, function(x) sum(x$Icum))) / n)


###################################################
### code chunk number 35: SIR-mparse-incidence-plot (eval = FALSE)
###################################################
## plot(cases, main = "", xlab = "Time", ylab = "Number of cases",
##   do.points = FALSE)
## lines(avg_cases, col = "blue", lwd = 2, lty = 2)


###################################################
### code chunk number 36: SIR-mparse-incidence-plot
###################################################
transitions <- c("S -> b*S*I/(S+I+R) -> I + Icum", "I -> g*I -> R")
compartments <- c("S", "I", "Icum", "R")
n <- 1000
u0 <- data.frame(S = rep(99, n), I = rep(1, n), Icum = rep(0, n),
R = rep(0, n))
model <- mparse(transitions = transitions, compartments = compartments,
gdata = c(b = 0.16, g = 0.077), u0 = u0, tspan = 1:150)
set.seed(123)
set_num_threads(1)
result <- run(model = model)
traj <- trajectory(model = result, compartments = "Icum")
cases <- stepfun(result@tspan[-1], diff(c(0, traj$Icum[traj$node == 1])))
avg_cases <- c(0, diff(by(traj, traj$time, function(x) sum(x$Icum))) / n)
plot(cases, main = "", xlab = "Time", ylab = "Number of cases",
do.points = FALSE)
lines(avg_cases, col = "blue", lwd = 2, lty = 2)


###################################################
### code chunk number 37: mparse-scheduled-events
###################################################
transitions <- c("S -> b*S*I/(S+I+R+V) -> I + Icum", "I -> g*I -> R")
compartments <- c("S", "I", "Icum", "R", "V")


###################################################
### code chunk number 38: mparse-scheduled-events-data
###################################################
u0 <- u0_SIR()
u0$Icum <- 0
u0$V <- 0
events <- events_SIR()


###################################################
### code chunk number 39: mparse-scheduled-events-vaccination
###################################################
vaccination <- data.frame(event = "intTrans", time = rep(21:52,
  each = 50), node = 1:1600, dest = 0, n = 0, proportion = 0.8,
  select = 3, shift = 1)


###################################################
### code chunk number 40: mparse-scheduled-events-E-N
###################################################
E <- matrix(c(1, 0, 0, 0, 0, 1, 1, 0, 1, 1, 1, 1, 0, 1, 0), nrow = 5,
  ncol = 3, dimnames = list(c("S", "I", "Icum", "R", "V"),
  c("1", "2", "3")))
N <- matrix(c(4, 3, 0, 1, 0), nrow = 5, ncol = 1,
  dimnames = list(c("S", "I", "Icum", "R", "V"), "1"))


###################################################
### code chunk number 41: mparse-recode-events
###################################################
events$select[events$select == 4] <- 2


###################################################
### code chunk number 42: mparse-scheduled-events-epicurve
###################################################
epicurve <- function(model, n = 1000) {
  Icum <- numeric(length(model@tspan))
  for (i in seq_len(n)) {
    model@u0["S", ] <- model@u0["S", ] + model@u0["I", ]
    model@u0["I", ] <- 0L
    j <- sample(seq_len(n_nodes(model)), 1)
    model@u0["I", j] <- 1L
    model@u0["S", j] <- model@u0["S", j] - 1L
    result <- run(model = model)
    traj <- trajectory(model = result, compartments = "Icum",
      format = "matrix")
    Icum <- Icum + colSums(traj)
  }
  stepfun(model@tspan[-1], diff(c(0, Icum / n)))
}


###################################################
### code chunk number 43: mparse-scheduled-events-model-no-vac (eval = FALSE)
###################################################
## model_no_vac <- mparse(transitions = transitions,
##   compartments = compartments, gdata = c(b = 0.16, g = 0.077),
##   u0 = u0, tspan = 1:300, events = events, E = E, N = N)
## cases_no_vac <- epicurve(model_no_vac)


###################################################
### code chunk number 44: mparse-scheduled-events-model-vac (eval = FALSE)
###################################################
## model_vac <- mparse(transitions = transitions,
##   compartments = compartments, gdata = c(b = 0.16, g = 0.077),
##   u0 = u0, tspan = 1:300, events = rbind(events, vaccination),
##   E = E, N = N)
## cases_vac <- epicurve(model_vac)


###################################################
### code chunk number 45: mparse-scheduled-events-plot (eval = FALSE)
###################################################
## plot(cases_no_vac, main = "", xlim = c(0, 300), xlab = "Time",
##   ylab = "Number of cases", do.points = FALSE)
## lines(cases_vac, col = "blue", do.points = FALSE, lty = 2)
## abline(v = 21, col = "red", lty = 3)
## legend("topright", c("No vaccination", "Vaccination"),
##   col = c("black", "blue"), lty = 1:2)


###################################################
### code chunk number 46: mparse-predator-prey-I
###################################################
transitions <- c("@ -> bR*R -> R", "R -> (dR+(bR-dR)*R/K)*R -> @",
  "R -> alpha/(1+w*R)*R*F -> @", "@ -> bF*alpha/(1+w*R)*R*F -> F",
  "F -> dF*F -> @")
compartments <- c("R", "F")
parameters <- c(bR = 2, bF = 2, dR = 1, K = 1000, alpha = 0.007,
  w = 0.0035, dF = 2)


###################################################
### code chunk number 47: mparse-predator-prey-II
###################################################
n <- 1000
u0 <- data.frame(R = rep(1000, n), F = rep(100, n))
model <- mparse(transitions = transitions, compartments = compartments,
  gdata = parameters, u0 = u0, tspan = 1:100)


###################################################
### code chunk number 48: create-package-skeleton (eval = FALSE)
###################################################
## path <- tempdir()
## package_skeleton(model = model, name = "PredatorPrey", path = path)


###################################################
### code chunk number 49: install-package-skeleton (eval = FALSE)
###################################################
## pkg <- file.path(path, "PredatorPrey")
## install.packages(pkg, repos = NULL, type = "source")


###################################################
### code chunk number 50: load-predator-prey (eval = FALSE)
###################################################
## library("PredatorPrey")


###################################################
### code chunk number 51: create-predator-prey-model (eval = FALSE)
###################################################
## model <- PredatorPrey(u0 = u0, tspan = 1:100, gdata = parameters)
## set.seed(123)
## set_num_threads(1)
## result <- run(model)


###################################################
### code chunk number 52: predator-prey-plot (eval = FALSE)
###################################################
## opar <- par(mfrow = c(1, 2))
## plot(R ~ F, trajectory(model = result), cex = 0.3, pch = 20,
##   xlab = "Number of predators", ylab = "Number of prey",
##   col = rgb(0, 0, 0, alpha = 0.1))
## plot(R ~ time, trajectory(model = result, index = 4), type = "l",
##   xlab = "Time", ylab = "N")
## lines(F ~ time, trajectory(model = result, index = 4), type = "l", lty = 2)
## legend("right", c("Prey", "Predator"), lty = 1:2)
## par(opar)


