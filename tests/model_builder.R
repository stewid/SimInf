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

library(SimInf)

## For debugging
sessionInfo()

## Check mparse
m <- mparse(transitions = c("D->c1*D->D+D", "D+W->c2*D*W->W+W","W->c3*W->@"),
            compartments = c("D","W"),
            c1 = 1, c2 = 0.005, c3 = 0.6)

latex <- c("\\begin{align}",
           "  \\left",
           "    \\begin{array}{rcl}",
           "      D & \\xrightarrow{c1*D} & D + D \\\\",
           "      D + W & \\xrightarrow{c2*D*W} & W + W \\\\",
           "      W & \\xrightarrow{c3*W} & \\emptyset \\\\",
           "    \\end{array}",
           "  \\right\\}",
           "\\end{align}")
stopifnot(identical(m@latex, latex))

G <- new("dgCMatrix",
         i = c(0L, 1L, 0L, 1L, 2L, 1L, 2L),
         p = c(0L, 2L, 5L, 7L),
         Dim = c(3L, 3L),
         Dimnames = list(c("D -> D + D", "D + W -> W + W", "W -> @"),
                         c("1", "2", "3")),
         x = c(1, 1, 1, 1, 1, 1, 1),
         factors = list())
stopifnot(identical(m@G, G))

S <- new("dgCMatrix",
         i = c(0L, 0L, 1L, 1L),
         p = c(0L, 1L, 3L, 4L),
         Dim = 2:3,
         Dimnames = list(c("D", "W"),
                         c("1", "2", "3")),
         x = c(1, -1, 1, -1),
         factors = list())
stopifnot(identical(m@S, S))

C_code <- c(
    "",
    "#include <R_ext/Rdynload.h>",
    "#include \"SimInf.h\"",
    "",
    "/* Compartments */",
    "enum{D, W};",
    "",
    "/* Rate constants */",
    "const double c1 = 1;",
    "const double c2 = 0.005;",
    "const double c3 = 0.6;",
    "",
    "double trFun1(",
    "    const int *u,",
    "    const double *v,",
    "    const double *ldata,",
    "    const double *gdata,",
    "    double t)", "{",
    "    return c1*u[D];",
    "}",
    "",
    "double trFun2(",
    "    const int *u,",
    "    const double *v,",
    "    const double *ldata,",
    "    const double *gdata,",
    "    double t)",
    "{",
    "    return c2*u[D]*u[W];",
    "}",
    "",
    "double trFun3(",
    "    const int *u,",
    "    const double *v,",
    "    const double *ldata,",
    "    const double *gdata,",
    "    double t)",
    "{",
    "    return c3*u[W];",
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
    "    return 0;", "}",
    "",
    "SEXP SimInf_model_run(SEXP model, SEXP threads, SEXP seed, SEXP solver)",
    "{",
    "    TRFun tr_fun[] = {&trFun1, &trFun2, &trFun3};",
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
                             propensity = "beta*u[S]*u[I]/(u[S]+u[I]+u[R])",
                             depends = c(1, 1, 1)),
                        .Names = c("orig_prop", "propensity", "depends"))))
