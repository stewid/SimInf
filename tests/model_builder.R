## SimInf, a framework for stochastic disease spread simulations
## Copyright (C) 2015 - 2017  Stefan Engblom
## Copyright (C) 2015 - 2017  Stefan Widgren
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

## Check mparse
m <- mparse(transitions = c("@->c1->D", "D->c2*D->D+D",
                            "D+W->c3*D*W->W+W","W->c4*W->@"),
            compartments = c("D","W"),
            c1 = 0.5, c2 = 1, c3 = 0.005, c4 = 0.6)

latex <- c("\\begin{align}",
           "  \\left",
           "    \\begin{array}{rcl}",
           "      \\emptyset & \\xrightarrow{c1} & D \\\\",
           "      D & \\xrightarrow{c2*D} & D + D \\\\",
           "      D + W & \\xrightarrow{c3*D*W} & W + W \\\\",
           "      W & \\xrightarrow{c4*W} & \\emptyset \\\\",
           "    \\end{array}",
           "  \\right\\}",
           "\\end{align}")
stopifnot(identical(m@latex, latex))

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
    "/* Rate constants */",
    "const double c1 = 0.5;",
    "const double c2 = 1;",
    "const double c3 = 0.005;",
    "const double c4 = 0.6;",
    "",
    "double trFun1(",
    "    const int *u,",
    "    const double *v,",
    "    const double *ldata,",
    "    const double *gdata,",
    "    double t)",
    "{",
    "    return c1;",
    "}",
    "",
    "double trFun2(",
    "    const int *u,",
    "    const double *v,",
    "    const double *ldata,",
    "    const double *gdata,",
    "    double t)",
    "{",
    "    return c2*u[0];",
    "}",
    "",
    "double trFun3(",
    "    const int *u,",
    "    const double *v,",
    "    const double *ldata,",
    "    const double *gdata,",
    "    double t)",
    "{",
    "    return c3*u[0]*u[1];",
    "}",
    "",
    "double trFun4(",
    "    const int *u,",
    "    const double *v,",
    "    const double *ldata,",
    "    const double *gdata,",
    "    double t)",
    "{",
    "    return c4*u[1];",
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
    "SEXP SimInf_model_run(SEXP model, SEXP threads, SEXP seed, SEXP solver)",
    "{",
    "    TRFun tr_fun[] = {&trFun1, &trFun2, &trFun3, &trFun4};",
    "    DL_FUNC SimInf_run = R_GetCCallable(\"SimInf\", \"SimInf_run\");",
    "    return SimInf_run(model, threads, seed, solver, tr_fun, &ptsFun);",
    "}",
    "")
stopifnot(identical(m@C_code[-1], C_code)) ## Skip first line that contains time

stopifnot(identical(SimInf:::tokens("beta*S*I/(S+I+R)"),
                    c("beta", "*", "S", "*", "I", "/", "(", "S", "+",
                      "I", "+", "R", ")")))

stopifnot(
    identical(SimInf:::rewriteprop("beta*S*I/(S+I+R)", c("S", "I", "R")),
              structure(list(orig_prop = "beta*S*I/(S+I+R)",
                             propensity = "beta*u[0]*u[1]/(u[0]+u[1]+u[2])",
                             depends = c(1, 1, 1)),
                        .Names = c("orig_prop", "propensity", "depends"))))

## Check init function
m <- mparse(c("S -> b*S*I/(S+I+R) -> I", "I -> g*I -> R"),
            c("S", "I", "R"), b = 0.16, g = 0.077)
model <- init(m, u0 = data.frame(S = 100, I = 1, R = 0), tspan = 1:10)
C_code <- c(
    "",
    "#include <R_ext/Rdynload.h>",
    "#include \"SimInf.h\"",
    "",
    "/* Rate constants */",
    "const double b = 0.16;",
    "const double g = 0.077;",
    "",
    "double trFun1(",
    "    const int *u,",
    "    const double *v,",
    "    const double *ldata,",
    "    const double *gdata,",
    "    double t)",
    "{",
    "    return b*u[0]*u[1]/(u[0]+u[1]+u[2]);",
    "}",
    "",
    "double trFun2(",
    "    const int *u,",
    "    const double *v,",
    "    const double *ldata,",
    "    const double *gdata,",
    "    double t)",
    "{",
    "    return g*u[1];",
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
    "SEXP SimInf_model_run(SEXP model, SEXP threads, SEXP seed, SEXP solver)",
    "{",
    "    TRFun tr_fun[] = {&trFun1, &trFun2};",
    "    DL_FUNC SimInf_run = R_GetCCallable(\"SimInf\", \"SimInf_run\");",
    "    return SimInf_run(model, threads, seed, solver, tr_fun, &ptsFun);",
    "}",
    "")
stopifnot(identical(model@C_code[-1], C_code)) ## Skip first line that contains time

u0 <- structure(c(100L, 1L, 0L),
                .Dim = c(3L, 1L),
                .Dimnames = list(c("S", "I", "R"), NULL))
stopifnot(identical(model@u0, u0))
