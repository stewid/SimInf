## -------------------------------------------------------------------
library(SimInf)

set.seed(123)
set_num_threads(1)

u0 <- data.frame(S = c(100, 101, 102, 103, 104, 105),
                 I = c(1, 2, 3, 4, 5, 6),
                 R = c(0, 0, 0, 0, 0, 0))

model  <- SIR(u0 = u0,
              tspan = 1:10,
              beta = 0.16,
              gamma = 0.077)

result <- run(model)

## -------------------------------------------------------------------
trajectory(result)

## -------------------------------------------------------------------
trajectory(result, compartments = "R", index = 1)

## -------------------------------------------------------------------
trajectory(result, compartments = "R", index = c(1, 3))

## -------------------------------------------------------------------
prevalence(result, I ~ S + I + R)

## -------------------------------------------------------------------
prevalence(result, I ~ .)

## -------------------------------------------------------------------
prevalence(result, I ~ S + I + R, level = 2)

## -------------------------------------------------------------------
prevalence(result, I ~ S + I + R, level = 3)

## ----fig.width=7, fig.height=4, fig.align="left"--------------------
plot(result)

## ----fig.width=7, fig.height=4, fig.align="left"--------------------
plot(result, range = 0.95)

## ----fig.width=7, fig.height=4, fig.align="left"--------------------
plot(result, "I")

## ----fig.width=7, fig.height=4, fig.align="left"--------------------
plot(result, ~I)

## ----fig.width=7, fig.height=4, fig.align="left"--------------------
plot(result, index = 1:3, range = FALSE)

## ----fig.width=7, fig.height=4, fig.align="left"--------------------
plot(result, index = 1:3, range = FALSE, type = "l")

## ----fig.width=7, fig.height=4, fig.align="left"--------------------
plot(result, "I", index = 1, range = FALSE)

## ----fig.width=7, fig.height=4, fig.align="left"--------------------
plot(result, I ~ S + I + R)

## ----fig.width=7, fig.height=4, fig.align="left"--------------------
plot(result, I ~ S + I + R, level = 2)

## ----fig.width=7, fig.height=4, fig.align="left"--------------------
plot(result, I ~ S + I + R, level = 3)

## ----fig.width=7, fig.height=4, fig.align="left"--------------------
plot(result, I ~ S + I + R, level = 3, index = 1:3, range = FALSE)

## ----eval=FALSE-----------------------------------------------------
# help("plot,SimInf_model-method", package = "SimInf")

