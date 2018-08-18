## SimInf, a framework for stochastic disease spread simulations
## Copyright (C) 2015 - 2018  Stefan Engblom
## Copyright (C) 2015 - 2018  Stefan Widgren
##
## This program is free software: you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.
##
## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with this program.  If not, see <http://www.gnu.org/licenses/>.

library("SimInf")

## For debugging
sessionInfo()

## Check that invalid arguments to mparse raises error
res <- tools::assertError(
                  mparse(compartments = c("D","W"),
                         gdata = c(c1 = 0.5, c2 = 1, c3 = 0.005, c4 = 0.6),
                         u0 = data.frame(D = 10, W = 10), tspan = 1:5))
stopifnot(length(grep("'transitions' must be specified in a character vector.",
                      res[[1]]$message, fixed = TRUE)) > 0)

res <- tools::assertError(
                  mparse(transitions = 5,
                         compartments = c("D","W"),
                         gdata = c(c1 = 0.5, c2 = 1, c3 = 0.005, c4 = 0.6),
                         u0 = data.frame(D = 10, W = 10), tspan = 1:5))
stopifnot(length(grep("'transitions' must be specified in a character vector.",
                      res[[1]]$message, fixed = TRUE)) > 0)

res <- tools::assertError(
                  mparse(transitions = c("@->c1->D", "D->c2*D->D+D",
                                         "D+W->c3*D*W->W+W","W->c4*W->@"),
                         gdata = c(c1 = 0.5, c2 = 1, c3 = 0.005, c4 = 0.6),
                         u0 = data.frame(D = 10, W = 10), tspan = 1:5))
stopifnot(length(grep("'compartments' must be specified in a character vector.",
                      res[[1]]$message, fixed = TRUE)) > 0)

res <- tools::assertError(
                  mparse(transitions = c("@->c1->D", "D->c2*D->D+D",
                                         "D+W->c3*D*W->W+W","W->c4*W->@"),
                         compartments = 5,
                         gdata = c(c1 = 0.5, c2 = 1, c3 = 0.005, c4 = 0.6),
                         u0 = data.frame(D = 10, W = 10), tspan = 1:5))
stopifnot(length(grep("'compartments' must be specified in a character vector.",
                      res[[1]]$message, fixed = TRUE)) > 0)

res <- tools::assertError(
                  mparse(transitions = c("@->c1->D", "D->c2*D->D+D",
                                         "D+W->c3*D*W->W+W","W->c4*W->@"),
                         compartments = c("D","W"),
                         gdata = letters,
                         u0 = data.frame(D = 10, W = 10), tspan = 1:5))
stopifnot(length(grep("'gdata' must be a named numeric vector with unique names.",
                      res[[1]]$message, fixed = TRUE)) > 0)

res <- tools::assertError(
                  mparse(transitions = c("@->c1->D", "D->c2*D->D+D",
                                         "D+W->c3*D*W->W+W","W->c4*W->@"),
                         compartments = c("D", "W", "D"),
                         gdata = c(c1 = 0.5, c2 = 1, c3 = 0.005, c4 = 0.6),
                         u0 = data.frame(D = 10, W = 10), tspan = 1:5))
stopifnot(length(grep("'compartments' must be specified in a character vector.",
                      res[[1]]$message, fixed = TRUE)) > 0)

res <- tools::assertError(
                  mparse(transitions = c("@->c1->D", "D->c2*D->D+D",
                                         "D+W->c3*D*W->W+W","W->c4*W->@"),
                         compartments = c("D","W"),
                         gdata = c(c1 = 0.5, c2 = 1, c3 = 0.005, c4 = 0.6, c1 = 2),
                         u0 = data.frame(D = 10, W = 10), tspan = 1:5))
stopifnot(length(grep("'gdata' must be a named numeric vector with unique names.",
                      res[[1]]$message, fixed = TRUE)) > 0)

res <- tools::assertError(
                  mparse(transitions = c("@->c1->D", "D->c2*D->D+D",
                                         "D+W->c3*D*W->W+W","W->c4*W->@"),
                         compartments = c("v_new", "u", "v", "ldata", "gdata",
                                          "node", "t", "rng"),
                         gdata = c(c1 = 0.5, c2 = 1, c3 = 0.005, c4 = 0.6),
                         u0 = data.frame(D = 10, W = 10), tspan = 1:5))
stopifnot(length(grep("Invalid compartment names: v_new, u, v, ldata, gdata, node, t, rng",
                      res[[1]]$message)) > 0)

res <- tools::assertError(
                  mparse(transitions = c("@->c1->D", "D->c2*D->D+D",
                                         "D+W->c3*D*W->W+W","W->c4*W->@"),
                         compartments = c("D","W"),
                         gdata = c(v_new = 0.5, u = 1, v = 0.005, ldata = 0.6,
                                   gdata = 2, node = 3, t = 4, rng = 5),
                         u0 = data.frame(D = 10, W = 10), tspan = 1:5))
stopifnot(length(grep("Invalid 'gdata' names: v_new, u, v, ldata, gdata, node, t, rng",
                      res[[1]]$message)) > 0)

res <- tools::assertError(
                  mparse(transitions = c("@->c1->D", "D->c2*D->D+D",
                                         "D+W->c3*D*W->W+W","W->c4*W"),
                         compartments = c("D","W"),
                         gdata = c(c1 = 0.5, c2 = 1, c3 = 0.005, c4 = 0.6),
                         u0 = data.frame(D = 10, W = 10), tspan = 1:5))
stopifnot(length(grep("Invalid transition: 'W->c4[*]W'",
                      res[[1]]$message)) > 0)

res <- tools::assertError(
                  mparse(transitions = c("A->c1->D", "D->c2*D->D+D",
                                         "D+W->c3*D*W->W+W","W->c4*W->@"),
                         compartments = c("D","W"),
                         gdata = c(c1 = 0.5, c2 = 1, c3 = 0.005, c4 = 0.6),
                         u0 = data.frame(D = 10, W = 10), tspan = 1:5))
stopifnot(length(grep("Unknown compartment: 'A'[.]",
                      res[[1]]$message)) > 0)

res <- tools::assertError(
                  mparse(transitions = c("@->c1->D", "D->c2*D->D+D",
                                         "D+W->c3*D*W->W+W","W->c4*W->B"),
                         compartments = c("D","W"),
                         gdata = c(c1 = 0.5, c2 = 1, c3 = 0.005, c4 = 0.6),
                         u0 = data.frame(D = 10, W = 10), tspan = 1:5))
stopifnot(length(grep("Unknown compartment: 'B'[.]",
                      res[[1]]$message)) > 0)

res <- tools::assertError(
                  mparse(transitions = c("@->c1->D", "D->c2*D->D+D",
                                         "D+W->c3*D*W->W+W","W->c4*W->@"),
                         compartments = c("D","W"),
                         gdata = c(c1 = 0.5, c2 = 1, c3 = 0.005, c4 = 0.6),
                         u0 = matrix(c(10, 10), nrow = 1, ncol = 2,
                                     dimnames = list(NULL, c("A", "W"))),
                         tspan = 1:5))
stopifnot(length(grep("Missing columns in u0",
                      res[[1]]$message)) > 0)


res <- tools::assertError(
                  mparse(transitions = c("@->c1->D", "D->c2*D->D+D",
                                         "D+W->c3*D*W->W+W","W->c4*W->@"),
                         compartments = c("D","W"),
                         ldata = 1:5,
                         gdata = c(c1 = 0.5, c2 = 1, c3 = 0.005, c4 = 0.6),
                         u0 = data.frame(D = rep(10, 5), W = 10), tspan = 1:5))
stopifnot(length(grep("'ldata' must be a numeric matrix with non-duplicated rownames.",
                      res[[1]]$message)) > 0)

res <- tools::assertError(
                  mparse(transitions = c("@->c1->D", "D->c2*D->D+D",
                                         "D+W->c3*D*W->W+W","W->c4*W->@"),
                         compartments = c("D","W"),
                         ldata = matrix(1:5,, nrow = 1, dimnames = list("u", NULL)),
                         gdata = c(c1 = 0.5, c2 = 1, c3 = 0.005, c4 = 0.6),
                         u0 = data.frame(D = rep(10, 5), W = 10), tspan = 1:5))
stopifnot(length(grep("Invalid 'ldata' rownames: u",
                      res[[1]]$message)) > 0)

res <- tools::assertError(
                  mparse(transitions = c("@->c1->D", "D->c2*D->D+D",
                                         "D+W->c3*D*W->W+W","W->c4*W->@"),
                         compartments = c("D","W"),
                         ldata = matrix(1:5,, nrow = 1, dimnames = list("c4", NULL)),
                         gdata = c(c1 = 0.5, c2 = 1, c3 = 0.005, c4 = 0.6),
                         u0 = data.frame(D = rep(10, 5), W = 10), tspan = 1:5))
stopifnot(length(grep("'gdata' names and 'ldata' rownames have elements in common.",
                      res[[1]]$message)) > 0)

## Check mparse
m <- mparse(transitions = c("@->c1->D", "D->c2*D->D+D",
                            "D+W->c3*D*W->W+W","W->c4*W->@"),
            compartments = c("D","W"),
            gdata = c(c1 = 0.5, c2 = 1, c3 = 0.005, c4 = 0.6),
            u0 = data.frame(D = 10, W = 10), tspan = 1:5)

G <- new("dgCMatrix",
         i = c(1L, 2L, 1L, 2L, 1L, 2L, 3L, 2L, 3L),
         p = c(0L, 2L, 4L, 7L, 9L),
         Dim = c(4L, 4L),
         Dimnames = list(c("@ -> D", "D -> D + D",
                           "D + W -> W + W", "W -> @"),
                         c("1", "2", "3", "4")),
         x = c(1, 1, 1, 1, 1, 1, 1, 1, 1),
         factors = list())
stopifnot(identical(m@G, G))

S <- new("dgCMatrix",
         i = c(0L, 0L, 0L, 1L, 1L),
         p = c(0L, 1L, 2L, 4L, 5L),
         Dim = c(2L, 4L),
         Dimnames = list(c("D", "W"),
                         c("1", "2", "3", "4")),
         x = c(1, 1, -1, 1, -1),
         factors = list())
stopifnot(identical(m@S, S))

C_code <- c(
    "",
    "#include <R_ext/Rdynload.h>",
    "#include \"SimInf.h\"",
    "",
    "double trFun1(",
    "    const int *u,",
    "    const double *v,",
    "    const double *ldata,",
    "    const double *gdata,",
    "    double t)",
    "{",
    "    return gdata[0];",
    "}",
    "",
    "double trFun2(",
    "    const int *u,",
    "    const double *v,",
    "    const double *ldata,",
    "    const double *gdata,",
    "    double t)",
    "{",
    "    return gdata[1]*u[0];",
    "}",
    "",
    "double trFun3(",
    "    const int *u,",
    "    const double *v,",
    "    const double *ldata,",
    "    const double *gdata,",
    "    double t)",
    "{",
    "    return gdata[2]*u[0]*u[1];",
    "}",
    "",
    "double trFun4(",
    "    const int *u,",
    "    const double *v,",
    "    const double *ldata,",
    "    const double *gdata,",
    "    double t)",
    "{",
    "    return gdata[3]*u[1];",
    "}",
    "",
    "int ptsFun(",
    "    double *v_new,",
    "    const int *u,",
    "    const double *v,",
    "    const double *ldata,",
    "    const double *gdata,",
    "    int node,",
    "    double t)",
    "{",
    "    return 0;",
    "}",
    "",
    "SEXP SimInf_model_run(SEXP model, SEXP threads, SEXP solver)",
    "{",
    "    TRFun tr_fun[] = {&trFun1, &trFun2, &trFun3, &trFun4};",
    "    DL_FUNC SimInf_run = R_GetCCallable(\"SimInf\", \"SimInf_run\");",
    "    return SimInf_run(model, threads, solver, tr_fun, &ptsFun);",
    "}",
    "")
stopifnot(identical(m@C_code[-1], C_code)) ## Skip first line that contains time

## Check mparse with both gdata and ldata
m <- mparse(transitions = c("@->c1->D", "D->c2*D->D+D",
                            "D+W->c3*D*W->W+W","W->c4*W->@"),
            compartments = c("D","W"),
            ldata = matrix(rep(0.6, 5), nrow = 1, dimnames = list("c4", NULL)),
            gdata = c(c1 = 0.5, c2 = 1, c3 = 0.005),
            u0 = data.frame(D = rep(10, 5), W = 10), tspan = 1:5)

C_code <- c(
    "",
    "#include <R_ext/Rdynload.h>",
    "#include \"SimInf.h\"",
    "",
    "double trFun1(",
    "    const int *u,",
    "    const double *v,",
    "    const double *ldata,",
    "    const double *gdata,",
    "    double t)",
    "{",
    "    return gdata[0];",
    "}",
    "",
    "double trFun2(",
    "    const int *u,",
    "    const double *v,",
    "    const double *ldata,",
    "    const double *gdata,",
    "    double t)",
    "{",
    "    return gdata[1]*u[0];",
    "}",
    "",
    "double trFun3(",
    "    const int *u,",
    "    const double *v,",
    "    const double *ldata,",
    "    const double *gdata,",
    "    double t)",
    "{",
    "    return gdata[2]*u[0]*u[1];",
    "}",
    "",
    "double trFun4(",
    "    const int *u,",
    "    const double *v,",
    "    const double *ldata,",
    "    const double *gdata,",
    "    double t)",
    "{",
    "    return ldata[0]*u[1];",
    "}",
    "",
    "int ptsFun(",
    "    double *v_new,",
    "    const int *u,",
    "    const double *v,",
    "    const double *ldata,",
    "    const double *gdata,",
    "    int node,",
    "    double t)",
    "{",
    "    return 0;",
    "}",
    "",
    "SEXP SimInf_model_run(SEXP model, SEXP threads, SEXP solver)",
    "{",
    "    TRFun tr_fun[] = {&trFun1, &trFun2, &trFun3, &trFun4};",
    "    DL_FUNC SimInf_run = R_GetCCallable(\"SimInf\", \"SimInf_run\");",
    "    return SimInf_run(model, threads, solver, tr_fun, &ptsFun);",
    "}",
    "")
stopifnot(identical(m@C_code[-1], C_code)) ## Skip first line that contains time

stopifnot(identical(SimInf:::tokens("beta*S*I/(S+I+R)"),
                    c("beta", "*", "S", "*", "I", "/", "(", "S", "+",
                      "I", "+", "R", ")")))

stopifnot(
    identical(SimInf:::rewriteprop("beta*S*I/(S+I+R)", c("S", "I", "R"), NULL, "beta"),
              structure(list(orig_prop = "beta*S*I/(S+I+R)",
                             propensity = "gdata[0]*u[0]*u[1]/(u[0]+u[1]+u[2])",
                             depends = c(1, 1, 1)),
                        .Names = c("orig_prop", "propensity", "depends"))))

## Check init function
model <- mparse(transitions = c("S -> b*S*I/(S+I+R) -> I",
                                "I -> g*I -> R"),
                compartments = c("S", "I", "R"),
                gdata = c(b = 0.16, g = 0.077),
                u0 = data.frame(S = 100, I = 1, R = 0),
                tspan = 1:10)
C_code <- c(
    "",
    "#include <R_ext/Rdynload.h>",
    "#include \"SimInf.h\"",
    "",
    "double trFun1(",
    "    const int *u,",
    "    const double *v,",
    "    const double *ldata,",
    "    const double *gdata,",
    "    double t)",
    "{",
    "    return gdata[0]*u[0]*u[1]/(u[0]+u[1]+u[2]);",
    "}",
    "",
    "double trFun2(",
    "    const int *u,",
    "    const double *v,",
    "    const double *ldata,",
    "    const double *gdata,",
    "    double t)",
    "{",
    "    return gdata[1]*u[1];",
    "}",
    "",
    "int ptsFun(",
    "    double *v_new,",
    "    const int *u,",
    "    const double *v,",
    "    const double *ldata,",
    "    const double *gdata,",
    "    int node,",
    "    double t)",
    "{",
    "    return 0;",
    "}",
    "",
    "SEXP SimInf_model_run(SEXP model, SEXP threads, SEXP solver)",
    "{",
    "    TRFun tr_fun[] = {&trFun1, &trFun2};",
    "    DL_FUNC SimInf_run = R_GetCCallable(\"SimInf\", \"SimInf_run\");",
    "    return SimInf_run(model, threads, solver, tr_fun, &ptsFun);",
    "}",
    "")
stopifnot(identical(model@C_code[-1], C_code)) ## Skip first line that contains time

u0 <- structure(c(100L, 1L, 0L),
                .Dim = c(3L, 1L),
                .Dimnames = list(c("S", "I", "R"), NULL))
stopifnot(identical(model@u0, u0))

## Check mparse with ldata and gdata as data.frames
m1 <- mparse(transitions = c("@->c1->D", "D->c2*D->D+D",
                             "D+W->c3*D*W->W+W","W->c4*W->@"),
             compartments = c("D","W"),
             ldata = matrix(c(0.2, 0.3, 0.4, 0.5, 0.6), nrow = 1,
                            dimnames = list("c4", NULL)),
             gdata = c(c1 = 0.5, c2 = 1, c3 = 0.005),
             u0 = data.frame(D = rep(10, 5), W = 10), tspan = 1:5)

m2 <- mparse(transitions = c("@->c1->D", "D->c2*D->D+D",
                             "D+W->c3*D*W->W+W","W->c4*W->@"),
             compartments = c("D","W"),
             ldata = data.frame(c4 = c(0.2, 0.3, 0.4, 0.5, 0.6)),
             gdata = data.frame(c1 = 0.5, c2 = 1, c3 = 0.005),
             u0 = data.frame(D = rep(10, 5), W = 10), tspan = 1:5)

stopifnot(identical(m1, m2))

## Check that mparse fails with gdata as a 2-row data.frame
res <- tools::assertError(
                  mparse(transitions = c("@->c1->D", "D->c2*D->D+D",
                                         "D+W->c3*D*W->W+W","W->c4*W->@"),
                         compartments = c("D","W"),
                         ldata = matrix(rep(0.6, 5), nrow = 1, dimnames = list("c4", NULL)),
                         gdata = data.frame(c1 = rep(0.5, 2), c2 = 1, c3 = 0.005),
                         u0 = data.frame(D = rep(10, 5), W = 10), tspan = 1:5))
stopifnot(length(grep("When 'gdata' is a data.frame, it must have one row",
                      res[[1]]$message)) > 0)

## Check mparse fails with ldata as data.frames and nrow(ldata) !=
## nrow(u0)
res <- tools::assertError(
                  mparse(transitions = c("@->c1->D", "D->c2*D->D+D",
                                         "D+W->c3*D*W->W+W","W->c4*W->@"),
                         compartments = c("D","W"),
                         ldata = data.frame(c4 = c(0.2, 0.3, 0.4, 0.5)),
                         gdata = data.frame(c1 = 0.5, c2 = 1, c3 = 0.005),
                         u0 = data.frame(D = rep(10, 5), W = 10), tspan = 1:5))
stopifnot(length(grep("'ldata' and 'u0' must have the same number of rows.",
                      res[[1]]$message)) > 0)
