## This file is part of SimInf, a framework for stochastic
## disease spread simulations.
##
## Copyright (C) 2015 Pavol Bauer
## Copyright (C) 2017 -- 2019 Robin Eriksson
## Copyright (C) 2015 -- 2019 Stefan Engblom
## Copyright (C) 2015 -- 2026 Stefan Widgren
##
## SimInf is free software: you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.
##
## SimInf is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with this program.  If not, see <https://www.gnu.org/licenses/>.

library(SimInf)
library(Matrix)
library(tools)
source("util/check.R")

## Specify the number of threads to use.
set_num_threads(1)

## For debugging
sessionInfo()

## Check that invalid arguments to mparse raises error
res <- assertError(
    mparse(compartments = c("D", "W"),
           gdata = c(c1 = 0.5, c2 = 1, c3 = 0.005, c4 = 0.6),
           u0 = data.frame(D = 10, W = 10), tspan = 1:5))
check_error(res, "'transitions' must be specified in a character vector.")

res <- assertError(
    mparse(transitions = 5, compartments = c("D", "W"),
           gdata = c(c1 = 0.5, c2 = 1, c3 = 0.005, c4 = 0.6),
           u0 = data.frame(D = 10, W = 10), tspan = 1:5))
check_error(res, "'transitions' must be specified in a character vector.")

res <- assertError(
    mparse(transitions = c("@->c1->D", "D->c2*D->D+D",
                           "D+W->c3*D*W->W+W", "W->c4*W->@"),
           gdata = c(c1 = 0.5, c2 = 1, c3 = 0.005, c4 = 0.6),
           u0 = data.frame(D = 10, W = 10), tspan = 1:5))
check_error(res, "'compartments' must be specified in a character vector.")

res <- assertError(
    mparse(transitions = c("@->c1->D", "D->c2*D->D+D",
                           "D+W->c3*D*W->W+W", "W->c4*W->@"),
           compartments = 5,
           gdata = c(c1 = 0.5, c2 = 1, c3 = 0.005, c4 = 0.6),
           u0 = data.frame(D = 10, W = 10), tspan = 1:5))
check_error(res, "'compartments' must be specified in a character vector.")

res <- assertError(
    mparse(transitions = c("@->c1->D", "D->c2*D->D+D",
                           "D+W->c3*D*W->W+W", "W->c4*W->@"),
           compartments = c("D", "W"), gdata = letters,
           u0 = data.frame(D = 10, W = 10), tspan = 1:5))
check_error(res, "'gdata' must either be a 'data.frame' or a 'numeric' vector.")

res <- assertError(
    mparse(transitions = c("@->c1->D", "D->c2*D->D+D",
                           "D+W->c3*D*W->W+W", "W->c4*W->@"),
           compartments = c("D", "W", "D"),
           gdata = c(c1 = 0.5, c2 = 1, c3 = 0.005, c4 = 0.6),
           u0 = data.frame(D = 10, W = 10), tspan = 1:5))
check_error(res, "'compartments' must be specified in a character vector.")

res <- assertError(
    mparse(transitions = c("@->c1->D", "D->c2*D->D+D",
                           "D+W->c3*D*W->W+W", "W->c4*W->@"),
           compartments = c("D", "W"),
           gdata = c(c1 = 0.5, c2 = 1, c3 = 0.005, c4 = 0.6,
                     c1 = 2),
           u0 = data.frame(D = 10, W = 10), tspan = 1:5))
check_error(res, "'gdata' must have non-duplicated parameter names.")

res <- assertError(
    mparse(transitions = c("@->c1->D", "D->c2*D->D+D",
                           "D+W->c3*D*W->W+W", "W->c4*W"),
           compartments = c("D", "W"),
           gdata = c(c1 = 0.5, c2 = 1, c3 = 0.005, c4 = 0.6),
           u0 = data.frame(D = 10, W = 10), tspan = 1:5))
check_error(res, "Invalid transition: 'W->c4*W'.")

res <- assertError(
    mparse(transitions = c("A->c1->D", "D->c2*D->D+D",
                           "D+W->c3*D*W->W+W", "W->c4*W->@"),
           compartments = c("D", "W"),
           gdata = c(c1 = 0.5, c2 = 1, c3 = 0.005, c4 = 0.6),
           u0 = data.frame(D = 10, W = 10), tspan = 1:5))
check_error(res, "Unknown compartment: 'A'.")

res <- assertError(
    mparse(transitions = c("@->c1->D", "D->c2*D->D+D",
                           "D+W->c3*D*W->W+W", "W->c4*W->B"),
           compartments = c("D", "W"),
           gdata = c(c1 = 0.5, c2 = 1, c3 = 0.005, c4 = 0.6),
           u0 = data.frame(D = 10, W = 10), tspan = 1:5))
check_error(res, "Unknown compartment: 'B'.")

res <- assertError(
    mparse(transitions = c("@->c1->D", "D->c2*D->D+D",
                           "D+W->c3*D*W->W+W", "W->c4*W->@"),
           compartments = c("D", "W"),
           gdata = c(c1 = 0.5, c2 = 1, c3 = 0.005, c4 = 0.6),
           u0 = matrix(c(10, 10), nrow = 1, ncol = 2,
                       dimnames = list(NULL, c("A", "W"))),
           tspan = 1:5))
check_error(res, "Missing columns in u0.")

res <- assertError(
    mparse(transitions = c("@->c1->D", "D->c2*D->D+D",
                           "D+W->c3*D*W->W+W", "W->c4*W->@"),
           compartments = c("D", "W"), ldata = 1:5,
           gdata = c(c1 = 0.5, c2 = 1, c3 = 0.005, c4 = 0.6),
           u0 = data.frame(D = rep(10, 5), W = 10),
           tspan = 1:5))
check_error(res, "'ldata' must either be a 'data.frame' or a 'matrix'.")

res <- assertError(
    mparse(transitions = c("@->c1->D", "D->c2*D->D+D",
                           "D+W->c3*D*W->W+W", "W->c4*W->@"),
           compartments = c("D", "W"),
           ldata = matrix(rep(0, 10), nrow = 2, ncol = 5,
                          dimnames = list(c("c1", "c1"))),
           gdata = c(c2 = 1, c3 = 0.005, c4 = 0.6),
           u0 = data.frame(D = rep(10, 5), W = 10),
           tspan = 1:5))
check_error(res, "'ldata' must have non-duplicated parameter names.")

res <- assertError(
    mparse(transitions = c("@->c1->D", "D->c2*D->D+D",
                           "D+W->c3*D*W->W+W", "W->c4*W->@"),
           compartments = c("D", "W"), v0 = 1:5,
           gdata = c(c1 = 0.5, c2 = 1, c3 = 0.005, c4 = 0.6),
           u0 = data.frame(D = rep(10, 5), W = 10),
           tspan = 1:5))
check_error(res, "'v0' must either be a 'data.frame' or a 'matrix'.")

res <- assertError(
    mparse(transitions = c("@->c1->D", "D->c2*D->D+D",
                           "D+W->c3*D*W->W+W", "W->c4*W->@"),
           compartments = c("D", "W"),
           v0 = matrix(rep(0, 10), nrow = 2, ncol = 5,
                       dimnames = list(c("c1", "c1"))),
           gdata = c(c2 = 1, c3 = 0.005, c4 = 0.6),
           u0 = data.frame(D = rep(10, 5), W = 10),
           tspan = 1:5))
check_error(res, "'v0' must have non-duplicated parameter names.")

res <- assertError(
    mparse(transitions = c("@->c1->D", "D->c2*D->D+D",
                           "D+W->c3*D*W->W+W", "W->c4*W->@"),
           compartments = c("D", "W"),
           ldata = matrix(1:5, nrow = 1,
                          dimnames = list("c4", NULL)),
           gdata = c(c1 = 0.5, c2 = 1, c3 = 0.005, c4 = 0.6),
           u0 = data.frame(D = rep(10, 5), W = 10),
           tspan = 1:5))
check_error(res, "Duplicated compartment or variable name detected.")

res <- assertError(
    mparse(transitions = c("@->c1->D", "D->c2*D->D+D",
                           "D+W->c3*D*W->W+W", "W->c4*W->@"),
           compartments = c("D", "W"),
           gdata = c(c1 = 0.5, c2 = 1, c3 = 0.005, c4 = 0.6),
           u0 = data.frame(D = 10, W = 10), tspan = 1:5,
           pts_fun = 5))
check_error(res, "'pts_fun' must be a character vector.")

res <- assertError(
    mparse(transitions = c("S -> beta*S*I/N -> I",
                           "I -> gamma*I -> R",
                           "N <- S+I+R"),
           compartments = c("S", "I", "R", "N_COMPARTMENTS_U"),
           gdata = c(beta = 0.16, gamma = 0.077),
           u0 = data.frame(S = 100, I = 1, R = 0,
                           N_COMPARTMENTS_U = 3),
           tspan = 1:100))
check_error(res, "Invalid compartment or variable name.")

res <- assertError(
    mparse(transitions = c("S -> beta*S*I/N -> I",
                           "I -> gamma*I -> R",
                           "N <- S+I+R"),
           compartments = c("S", "I", "R"),
           gdata = c(beta = 0.16, gamma = 0.077, N_COMPARTMENTS_V = 2),
           u0 = data.frame(S = 100, I = 1, R = 0),
           tspan = 1:100))
check_error(res, "Invalid compartment or variable name.")

res <- assertError(
    mparse(transitions = c("S -> beta*S*I/N -> I",
                           "I -> gamma*I -> R",
                           "N <- S+I+R"),
           compartments = c("S", "I", "R", "N_COMPARTMENTS_CELL"),
           gdata = c(beta = 0.16, gamma = 0.077),
           u0 = data.frame(S = 100, I = 1, R = 0,
                           N_COMPARTMENTS_CELL = 3),
           tspan = 1:100))
check_error(res, "Invalid compartment or variable name.")

## Check mparse
m <- mparse(transitions = c("@->c1->D", "D->c2*D->D+D",
                            "D+W->c3*D*W->W+W", "W->c4*W->@"),
            compartments = c("D", "W"),
            gdata = c(c1 = 0.5, c2 = 1, c3 = 0.005, c4 = 0.6),
            u0 = data.frame(D = 10, W = 10), tspan = 1:5)

G <- new("dgCMatrix",
         i = c(1L, 2L, 1L, 2L, 1L, 2L, 3L, 2L, 3L),
         p = c(0L, 2L, 4L, 7L, 9L),
         Dim = c(4L, 4L),
         Dimnames = list(c("@ -> c1 -> D", "D -> c2*D -> 2*D",
                           "D + W -> c3*D*W -> 2*W", "W -> c4*W -> @"),
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
    "/**",
    " * Make sure the necessary macros are defined so that the",
    " * compiler can replace them when compiling the model.",
    " * 'SIMINF_MODEL_RUN' defines the function name of the function",
    " * that will be called from R to run a trajectory of the model.",
    " * 'SIMINF_R_INIT' is the name of the function that R will call",
    " * when this model is loaded into R. 'SIMINF_FORCE_SYMBOLS'",
    " * defines whether R allows the entry point for the run function",
    " * to be searched for as a character string.",
    " * If this file is compiled from SimInf (when calling run), the",
    " * macros are defined by SimInf before calling 'R CMD SHLIB'.",
    " * If this file is compiled as part of a package, then the",
    " * definitions are set in the variable 'PKG_CPPFLAGS' in",
    " * 'src/Makevars' and 'src/Makevars.in'.",
    " */",
    "#if !defined(SIMINF_MODEL_RUN)",
    "#  error Definition for 'SIMINF_MODEL_RUN' is missing.",
    "#endif",
    "#if !defined(SIMINF_R_INIT)",
    "#  error Definition for 'SIMINF_R_INIT' is missing.",
    "#endif",
    "#if !defined(SIMINF_FORCE_SYMBOLS)",
    "#  error Definition for 'SIMINF_FORCE_SYMBOLS' is missing.",
    "#endif",
    "",
    "/**",
    " * @param u The compartment state vector in the node.",
    " * @param v The continuous state vector in the node.",
    " * @param ldata The local data vector in the node.",
    " * @param gdata The global data vector.",
    " * @param t Current time.",
    " * @return propensity.",
    " */",
    "static double trFun1(",
    "    const int *u,",
    "    const double *v,",
    "    const double *ldata,",
    "    const double *gdata,",
    "    double t)",
    "{",
    "    return gdata[0];",
    "}",
    "",
    "/**",
    " * @param u The compartment state vector in the node.",
    " * @param v The continuous state vector in the node.",
    " * @param ldata The local data vector in the node.",
    " * @param gdata The global data vector.",
    " * @param t Current time.",
    " * @return propensity.",
    " */",
    "static double trFun2(",
    "    const int *u,",
    "    const double *v,",
    "    const double *ldata,",
    "    const double *gdata,",
    "    double t)",
    "{",
    "    return gdata[1]*u[0];",
    "}",
    "",
    "/**",
    " * @param u The compartment state vector in the node.",
    " * @param v The continuous state vector in the node.",
    " * @param ldata The local data vector in the node.",
    " * @param gdata The global data vector.",
    " * @param t Current time.",
    " * @return propensity.",
    " */",
    "static double trFun3(",
    "    const int *u,",
    "    const double *v,",
    "    const double *ldata,",
    "    const double *gdata,",
    "    double t)",
    "{",
    "    return gdata[2]*u[0]*u[1];",
    "}",
    "",
    "/**",
    " * @param u The compartment state vector in the node.",
    " * @param v The continuous state vector in the node.",
    " * @param ldata The local data vector in the node.",
    " * @param gdata The global data vector.",
    " * @param t Current time.",
    " * @return propensity.",
    " */",
    "static double trFun4(",
    "    const int *u,",
    "    const double *v,",
    "    const double *ldata,",
    "    const double *gdata,",
    "    double t)",
    "{",
    "    return gdata[3]*u[1];",
    "}",
    "",
    "/**",
    " * Post time step function.",
    " *",
    " * @param v_new If a continuous state vector is used by a model,",
    " *        this is the new continuous state vector in the node after",
    " *        the post time step.",
    " * @param u The compartment state vector in the node.",
    " * @param v The current continuous state vector in the node.",
    " * @param ldata The local data vector in the node.",
    " * @param gdata The global data vector that is common to all nodes.",
    " * @param node The node index. Note the node index is zero-based,",
    " *        i.e., the first node is 0.",
    " * @param t Current time in the simulation.",
    " * @return error code (<0), or 1 if node needs to update the",
    " *         transition rates, or 0 when it doesn't need to update",
    " *         the transition rates.",
    " */",
    "static int ptsFun(",
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
    "/**",
    " * Run a trajectory of the model.",
    " *",
    " * @param model The model.",
    " * @param solver The name of the numerical solver.",
    " * @return A model with a trajectory attached to it.",
    " */",
    "static SEXP SIMINF_MODEL_RUN(SEXP model, SEXP solver)",
    "{",
    "    static SEXP(*SimInf_run)(SEXP, SEXP, TRFun*, PTSFun) = NULL;",
    "    TRFun tr_fun[] = {&trFun1, &trFun2, &trFun3, &trFun4};",
    "",
    "    if (!SimInf_run) {",
    "        SimInf_run = (SEXP(*)(SEXP, SEXP, TRFun*, PTSFun))",
    "            R_GetCCallable(\"SimInf\", \"SimInf_run\");",
    "",
    "        if (!SimInf_run) {",
    "            Rf_error(\"Cannot find function 'SimInf_run'.\");",
    "        }",
    "    }",
    "",
    "    return SimInf_run(model, solver, tr_fun, &ptsFun);",
    "}",
    "",
    "/**",
    " * A NULL-terminated array of routines to register for the .Call",
    " * interface, see section '5.4 Registering native routines' in",
    " * the 'Writing R Extensions' manual.",
    " */",
    "static const R_CallMethodDef callMethods[] =",
    "{",
    "    SIMINF_CALLDEF(SIMINF_MODEL_RUN, 2),",
    "    {NULL, NULL, 0}",
    "};",
    "",
    "/**",
    " * This routine will be invoked when R loads the shared object/DLL,",
    " * see section '5.4 Registering native routines' in the",
    " * 'Writing R Extensions' manual.",
    " */",
    "void SIMINF_R_INIT(DllInfo *info)",
    "{",
    "    R_registerRoutines(info, NULL, callMethods, NULL, NULL);",
    "    R_useDynamicSymbols(info, FALSE);",
    "    R_forceSymbols(info, SIMINF_FORCE_SYMBOLS);",
    "}")

## Skip first line that contains version
stopifnot(identical(m@C_code[-1], C_code))
stopifnot(identical(m@C_code, C_code(m)))

## Check mparse with both gdata and ldata
m <- mparse(transitions = c("@->c1->D", "D->c2*D->D+D",
                            "D+W->c3*D*W->W+W", "W->c4*W->@"),
            compartments = c("D", "W"),
            ldata = matrix(rep(0.6, 5), nrow = 1, dimnames = list("c4", NULL)),
            gdata = c(c1 = 0.5, c2 = 1, c3 = 0.005),
            u0 = data.frame(D = rep(10, 5), W = 10), tspan = 1:5)

C_code <- c(
    "",
    "#include <R_ext/Rdynload.h>",
    "#include \"SimInf.h\"",
    "",
    "/**",
    " * Make sure the necessary macros are defined so that the",
    " * compiler can replace them when compiling the model.",
    " * 'SIMINF_MODEL_RUN' defines the function name of the function",
    " * that will be called from R to run a trajectory of the model.",
    " * 'SIMINF_R_INIT' is the name of the function that R will call",
    " * when this model is loaded into R. 'SIMINF_FORCE_SYMBOLS'",
    " * defines whether R allows the entry point for the run function",
    " * to be searched for as a character string.",
    " * If this file is compiled from SimInf (when calling run), the",
    " * macros are defined by SimInf before calling 'R CMD SHLIB'.",
    " * If this file is compiled as part of a package, then the",
    " * definitions are set in the variable 'PKG_CPPFLAGS' in",
    " * 'src/Makevars' and 'src/Makevars.in'.",
    " */",
    "#if !defined(SIMINF_MODEL_RUN)",
    "#  error Definition for 'SIMINF_MODEL_RUN' is missing.",
    "#endif",
    "#if !defined(SIMINF_R_INIT)",
    "#  error Definition for 'SIMINF_R_INIT' is missing.",
    "#endif",
    "#if !defined(SIMINF_FORCE_SYMBOLS)",
    "#  error Definition for 'SIMINF_FORCE_SYMBOLS' is missing.",
    "#endif",
    "",
    "/**",
    " * @param u The compartment state vector in the node.",
    " * @param v The continuous state vector in the node.",
    " * @param ldata The local data vector in the node.",
    " * @param gdata The global data vector.",
    " * @param t Current time.",
    " * @return propensity.",
    " */",
    "static double trFun1(",
    "    const int *u,",
    "    const double *v,",
    "    const double *ldata,",
    "    const double *gdata,",
    "    double t)",
    "{",
    "    return gdata[0];",
    "}",
    "",
    "/**",
    " * @param u The compartment state vector in the node.",
    " * @param v The continuous state vector in the node.",
    " * @param ldata The local data vector in the node.",
    " * @param gdata The global data vector.",
    " * @param t Current time.",
    " * @return propensity.",
    " */",
    "static double trFun2(",
    "    const int *u,",
    "    const double *v,",
    "    const double *ldata,",
    "    const double *gdata,",
    "    double t)",
    "{",
    "    return gdata[1]*u[0];",
    "}",
    "",
    "/**",
    " * @param u The compartment state vector in the node.",
    " * @param v The continuous state vector in the node.",
    " * @param ldata The local data vector in the node.",
    " * @param gdata The global data vector.",
    " * @param t Current time.",
    " * @return propensity.",
    " */",
    "static double trFun3(",
    "    const int *u,",
    "    const double *v,",
    "    const double *ldata,",
    "    const double *gdata,",
    "    double t)",
    "{",
    "    return gdata[2]*u[0]*u[1];",
    "}",
    "",
    "/**",
    " * @param u The compartment state vector in the node.",
    " * @param v The continuous state vector in the node.",
    " * @param ldata The local data vector in the node.",
    " * @param gdata The global data vector.",
    " * @param t Current time.",
    " * @return propensity.",
    " */",
    "static double trFun4(",
    "    const int *u,",
    "    const double *v,",
    "    const double *ldata,",
    "    const double *gdata,",
    "    double t)",
    "{",
    "    return ldata[0]*u[1];",
    "}",
    "",
    "/**",
    " * Post time step function.",
    " *",
    " * @param v_new If a continuous state vector is used by a model,",
    " *        this is the new continuous state vector in the node after",
    " *        the post time step.",
    " * @param u The compartment state vector in the node.",
    " * @param v The current continuous state vector in the node.",
    " * @param ldata The local data vector in the node.",
    " * @param gdata The global data vector that is common to all nodes.",
    " * @param node The node index. Note the node index is zero-based,",
    " *        i.e., the first node is 0.",
    " * @param t Current time in the simulation.",
    " * @return error code (<0), or 1 if node needs to update the",
    " *         transition rates, or 0 when it doesn't need to update",
    " *         the transition rates.",
    " */",
    "static int ptsFun(",
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
    "/**",
    " * Run a trajectory of the model.",
    " *",
    " * @param model The model.",
    " * @param solver The name of the numerical solver.",
    " * @return A model with a trajectory attached to it.",
    " */",
    "static SEXP SIMINF_MODEL_RUN(SEXP model, SEXP solver)",
    "{",
    "    static SEXP(*SimInf_run)(SEXP, SEXP, TRFun*, PTSFun) = NULL;",
    "    TRFun tr_fun[] = {&trFun1, &trFun2, &trFun3, &trFun4};",
    "",
    "    if (!SimInf_run) {",
    "        SimInf_run = (SEXP(*)(SEXP, SEXP, TRFun*, PTSFun))",
    "            R_GetCCallable(\"SimInf\", \"SimInf_run\");",
    "",
    "        if (!SimInf_run) {",
    "            Rf_error(\"Cannot find function 'SimInf_run'.\");",
    "        }",
    "    }",
    "",
    "    return SimInf_run(model, solver, tr_fun, &ptsFun);",
    "}",
    "",
    "/**",
    " * A NULL-terminated array of routines to register for the .Call",
    " * interface, see section '5.4 Registering native routines' in",
    " * the 'Writing R Extensions' manual.",
    " */",
    "static const R_CallMethodDef callMethods[] =",
    "{",
    "    SIMINF_CALLDEF(SIMINF_MODEL_RUN, 2),",
    "    {NULL, NULL, 0}",
    "};",
    "",
    "/**",
    " * This routine will be invoked when R loads the shared object/DLL,",
    " * see section '5.4 Registering native routines' in the",
    " * 'Writing R Extensions' manual.",
    " */",
    "void SIMINF_R_INIT(DllInfo *info)",
    "{",
    "    R_registerRoutines(info, NULL, callMethods, NULL, NULL);",
    "    R_useDynamicSymbols(info, FALSE);",
    "    R_forceSymbols(info, SIMINF_FORCE_SYMBOLS);",
    "}")

## Skip first line that contains version
stopifnot(identical(m@C_code[-1], C_code))

stopifnot(identical(SimInf:::tokenize("beta*S*I/(S+I+R)"),
                    c("beta", "*", "S", "*", "I", "/", "(", "S", "+",
                      "I", "+", "R", ")")))

stopifnot(
    identical(SimInf:::rewrite_propensity(propensity = "beta*S*I/(S+I+R)",
                                          variables = list(),
                                          compartments = c("S", "I", "R"),
                                          cell_compartments = character(0),
                                          ldata_names  = NULL,
                                          gdata_names = "beta",
                                          v0_names = NULL,
                                          use_enum = FALSE),
              list(code = "gdata[0]*u[0]*u[1]/(u[0]+u[1]+u[2])",
                   depends = c(1, 1, 1),
                   G_rowname = "beta*S*I/(S+I+R)",
                   variables = character(0))))

stopifnot(
    identical(SimInf:::rewrite_propensity(propensity = "beta*S*I/(S+I+R)",
                                          variables = list(),
                                          compartments = c("S", "I", "R"),
                                          cell_compartments = character(0),
                                          ldata_names = NULL,
                                          gdata_names = "beta",
                                          v0_names = NULL,
                                          use_enum = TRUE),
              list(code = "gdata[BETA]*u[S]*u[I]/(u[S]+u[I]+u[R])",
                   depends = c(1, 1, 1),
                   G_rowname = "beta*S*I/(S+I+R)",
                   variables = character(0))))

stopifnot(
    identical(SimInf:::rewrite_propensity(
                           propensity = "beta*S*cell.contamination/(S+I+R)",
                           variables = list(),
                           compartments = c("S", "I", "R"),
                           cell_compartments = "cell.contamination",
                           ldata_names  = NULL,
                           gdata_names = "beta",
                           v0_names = NULL,
                           use_enum = TRUE),
              list(code =
                       "gdata[BETA]*u[S]*cell[CONTAMINATION]/(u[S]+u[I]+u[R])",
                   depends = c(1, 1, 1, 1),
                   G_rowname = "beta*S*cell.contamination/(S+I+R)",
                   variables = character(0))))

stopifnot(
    identical(SimInf:::rewrite_propensity(
                           propensity = "beta*S*cell.contamination/(S+I+R)",
                           variables = list(),
                           compartments = c("S", "I", "R"),
                           cell_compartments = "cell.contamination",
                           ldata_names  = NULL,
                           gdata_names = "beta",
                           v0_names = NULL,
                           use_enum = FALSE),
              list(code = "gdata[0]*u[0]*cell[0]/(u[0]+u[1]+u[2])",
                   depends = c(1, 1, 1, 1),
                   G_rowname = "beta*S*cell.contamination/(S+I+R)",
                   variables = character(0))))

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
    "/**",
    " * Make sure the necessary macros are defined so that the",
    " * compiler can replace them when compiling the model.",
    " * 'SIMINF_MODEL_RUN' defines the function name of the function",
    " * that will be called from R to run a trajectory of the model.",
    " * 'SIMINF_R_INIT' is the name of the function that R will call",
    " * when this model is loaded into R. 'SIMINF_FORCE_SYMBOLS'",
    " * defines whether R allows the entry point for the run function",
    " * to be searched for as a character string.",
    " * If this file is compiled from SimInf (when calling run), the",
    " * macros are defined by SimInf before calling 'R CMD SHLIB'.",
    " * If this file is compiled as part of a package, then the",
    " * definitions are set in the variable 'PKG_CPPFLAGS' in",
    " * 'src/Makevars' and 'src/Makevars.in'.",
    " */",
    "#if !defined(SIMINF_MODEL_RUN)",
    "#  error Definition for 'SIMINF_MODEL_RUN' is missing.",
    "#endif",
    "#if !defined(SIMINF_R_INIT)",
    "#  error Definition for 'SIMINF_R_INIT' is missing.",
    "#endif",
    "#if !defined(SIMINF_FORCE_SYMBOLS)",
    "#  error Definition for 'SIMINF_FORCE_SYMBOLS' is missing.",
    "#endif",
    "",
    "/**",
    " * @param u The compartment state vector in the node.",
    " * @param v The continuous state vector in the node.",
    " * @param ldata The local data vector in the node.",
    " * @param gdata The global data vector.",
    " * @param t Current time.",
    " * @return propensity.",
    " */",
    "static double trFun1(",
    "    const int *u,",
    "    const double *v,",
    "    const double *ldata,",
    "    const double *gdata,",
    "    double t)",
    "{",
    "    return gdata[0]*u[0]*u[1]/(u[0]+u[1]+u[2]);",
    "}",
    "",
    "/**",
    " * @param u The compartment state vector in the node.",
    " * @param v The continuous state vector in the node.",
    " * @param ldata The local data vector in the node.",
    " * @param gdata The global data vector.",
    " * @param t Current time.",
    " * @return propensity.",
    " */",
    "static double trFun2(",
    "    const int *u,",
    "    const double *v,",
    "    const double *ldata,",
    "    const double *gdata,",
    "    double t)",
    "{",
    "    return gdata[1]*u[1];",
    "}",
    "",
    "/**",
    " * Post time step function.",
    " *",
    " * @param v_new If a continuous state vector is used by a model,",
    " *        this is the new continuous state vector in the node after",
    " *        the post time step.",
    " * @param u The compartment state vector in the node.",
    " * @param v The current continuous state vector in the node.",
    " * @param ldata The local data vector in the node.",
    " * @param gdata The global data vector that is common to all nodes.",
    " * @param node The node index. Note the node index is zero-based,",
    " *        i.e., the first node is 0.",
    " * @param t Current time in the simulation.",
    " * @return error code (<0), or 1 if node needs to update the",
    " *         transition rates, or 0 when it doesn't need to update",
    " *         the transition rates.",
    " */",
    "static int ptsFun(",
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
    "/**",
    " * Run a trajectory of the model.",
    " *",
    " * @param model The model.",
    " * @param solver The name of the numerical solver.",
    " * @return A model with a trajectory attached to it.",
    " */",
    "static SEXP SIMINF_MODEL_RUN(SEXP model, SEXP solver)",
    "{",
    "    static SEXP(*SimInf_run)(SEXP, SEXP, TRFun*, PTSFun) = NULL;",
    "    TRFun tr_fun[] = {&trFun1, &trFun2};",
    "",
    "    if (!SimInf_run) {",
    "        SimInf_run = (SEXP(*)(SEXP, SEXP, TRFun*, PTSFun))",
    "            R_GetCCallable(\"SimInf\", \"SimInf_run\");",
    "",
    "        if (!SimInf_run) {",
    "            Rf_error(\"Cannot find function 'SimInf_run'.\");",
    "        }",
    "    }",
    "",
    "    return SimInf_run(model, solver, tr_fun, &ptsFun);",
    "}",
    "",
    "/**",
    " * A NULL-terminated array of routines to register for the .Call",
    " * interface, see section '5.4 Registering native routines' in",
    " * the 'Writing R Extensions' manual.",
    " */",
    "static const R_CallMethodDef callMethods[] =",
    "{",
    "    SIMINF_CALLDEF(SIMINF_MODEL_RUN, 2),",
    "    {NULL, NULL, 0}",
    "};",
    "",
    "/**",
    " * This routine will be invoked when R loads the shared object/DLL,",
    " * see section '5.4 Registering native routines' in the",
    " * 'Writing R Extensions' manual.",
    " */",
    "void SIMINF_R_INIT(DllInfo *info)",
    "{",
    "    R_registerRoutines(info, NULL, callMethods, NULL, NULL);",
    "    R_useDynamicSymbols(info, FALSE);",
    "    R_forceSymbols(info, SIMINF_FORCE_SYMBOLS);",
    "}")

## Skip first line that contains version
stopifnot(identical(model@C_code[-1], C_code))

u0 <- structure(c(100L, 1L, 0L),
                .Dim = c(3L, 1L),
                .Dimnames = list(c("S", "I", "R"), NULL))
stopifnot(identical(model@u0, u0))

## Check mparse with ldata and gdata as data.frames
m1 <- mparse(transitions = c("@->c1->D", "D->c2*D->D+D",
                             "D+W->c3*D*W->W+W", "W->c4*W->@"),
             compartments = c("D", "W"),
             ldata = matrix(c(2L, 3L, 4L, 5L, 6L), nrow = 1,
                            dimnames = list("c4", NULL)),
             gdata = c(c1 = 0.5, c2 = 1, c3 = 0.005),
             u0 = data.frame(D = rep(10, 5), W = 10), tspan = 1:5)

m2 <- mparse(transitions = c("@->c1->D", "D->c2*D->D+D",
                             "D+W->c3*D*W->W+W", "W->c4*W->@"),
             compartments = c("D", "W"),
             ldata = data.frame(c4 = c(2L, 3L, 4L, 5L, 6L)),
             gdata = data.frame(c1 = 0.5, c2 = 1, c3 = 0.005),
             u0 = data.frame(D = rep(10, 5), W = 10), tspan = 1:5)

## Skip first line that contains version
m1@C_code <- m1@C_code[-1]
m2@C_code <- m2@C_code[-1]

stopifnot(identical(m1, m2))

## Check that mparse fails with gdata as a 2-row data.frame
res <- assertError(
    mparse(transitions = c("@->c1->D", "D->c2*D->D+D",
                           "D+W->c3*D*W->W+W", "W->c4*W->@"),
           compartments = c("D", "W"),
           ldata = matrix(rep(0.6, 5), nrow = 1,
                          dimnames = list("c4", NULL)),
           gdata = data.frame(c1 = rep(0.5, 2), c2 = 1,
                              c3 = 0.005),
           u0 = data.frame(D = rep(10, 5), W = 10),
           tspan = 1:5))
check_error(res, "When 'gdata' is a data.frame, it must have one row.")

## Check mparse fails with ldata as data.frames and nrow(ldata) != nrow(u0)
res <- assertError(
    mparse(transitions = c("@->c1->D", "D->c2*D->D+D",
                           "D+W->c3*D*W->W+W", "W->c4*W->@"),
           compartments = c("D", "W"),
           ldata = data.frame(c4 = c(0.2, 0.3, 0.4, 0.5)),
           gdata = data.frame(c1 = 0.5, c2 = 1, c3 = 0.005),
           u0 = data.frame(D = rep(10, 5), W = 10),
           tspan = 1:5))
check_error(res, "The number of nodes in 'u0' and 'ldata' must match.", FALSE)

## Check 'S + S -> mu -> @'
m  <- mparse(transitions = "S + S -> mu -> @",
             compartments = c("S", "I"),
             gdata = c(mu = 1),
             u0 = data.frame(S = 100, I = 100),
             tspan = 1:100)

S <- new("dgCMatrix", i = 0L, p = 0:1, Dim = 2:1,
         Dimnames = list(c("S", "I"), "1"),
         x = -2, factors = list())
stopifnot(identical(m@S, S))

G <- new("dgCMatrix", i = integer(0), p = c(0L, 0L),
         Dim = c(1L, 1L), Dimnames = list("2*S -> mu -> @", "1"),
         x = numeric(0), factors = list())
stopifnot(identical(m@G, G))

## Check 'S + S -> mu -> S + S'
m  <- mparse(transitions = "S + S -> mu -> S + S",
             compartments = c("S", "I"),
             gdata = c(mu = 1),
             u0 = data.frame(S = 100, I = 100),
             tspan = 1:100)

S <- new("dgCMatrix", i = integer(0), p = c(0L, 0L),
         Dim = 2:1, Dimnames = list(c("S", "I"), "1"),
         x = numeric(0), factors = list())
stopifnot(identical(m@S, S))

G <- new("dgCMatrix", i = integer(0), p = c(0L, 0L),
         Dim = c(1L, 1L), Dimnames = list("2*S -> mu -> 2*S", "1"),
         x = numeric(0), factors = list())
stopifnot(identical(m@G, G))

## Check '@ -> mu-> S + S'
m  <- mparse(transitions = "@ -> mu-> S + S",
             compartments = c("S", "I"),
             gdata = c(mu = 1),
             u0 = data.frame(S = 100, I = 100),
             tspan = 1:100)

S <- new("dgCMatrix", i = 0L, p = 0:1, Dim = 2:1,
         Dimnames = list(c("S", "I"), "1"),
         x = 2, factors = list())
stopifnot(identical(m@S, S))

G <- new("dgCMatrix", i = integer(0), p = c(0L, 0L),
         Dim = c(1L, 1L), Dimnames = list("@ -> mu -> 2*S", "1"),
         x = numeric(0), factors = list())
stopifnot(identical(m@G, G))

## Check parsing replicates of compartments
stopifnot(identical(SimInf:::parse_compartments(x = "S + 2*S",
                                                compartments = c("S", "I"),
                                                cell_compartments = character(0)),
                    c(3L, 0L)))

## Check mparse with a compartment name that contains '.', for
## example, '.S.S' (this is a valid column name in a data.frame).
m  <- mparse(transitions = ".S.S -> 1.2*.S.S -> @",
             compartments = c(".S.S"),
             u0 = data.frame(.S.S = 100),
             tspan = 1:100)
stopifnot(identical(m@C_code[46], "    return 1.2*u[0];"))

## Check mparse with a propensity that contains '->' to handle a case
## where a pointer is used in the propensity.
m  <- mparse(transitions = "S -> a->data[2]*1.2*S -> @",
             compartments = c("S"),
             u0 = data.frame(S = 100),
             tspan = 1:100)
stopifnot(identical(m@C_code[46], "    return a->data[2]*1.2*u[0];"))

## Check that an error is raised if the compilation fails. Define a
## model with undeclared identifiers 'betaSI' and 'gammaI'.
model <- mparse(transitions = c("S -> betaSI -> I",
                                "I -> gammaI -> R"),
                compartments = c("S", "I", "R"),
                gdata = c(beta = 0.16, gamma = 0.077),
                u0 = data.frame(S = 100:105, I = 1:6, R = rep(0, 6)),
                tspan = 1:10)
assertError(run(model))

## Test that the environmental variable PKG_CPPFLAGS is restored after
## compiling C code for a model. First set it to a known value.
model <- mparse(transitions = c("S -> beta*S*I/(S+I+R) -> I",
                                "I -> gamma*I -> R"),
                compartments = c("S", "I", "R"),
                gdata = c(beta = 0.16, gamma = 0.077),
                u0 = data.frame(S = 100:105, I = 1:6, R = rep(0, 6)),
                tspan = 1:10)

Sys.setenv(PKG_CPPFLAGS = "test")

set.seed(22)
result <- run(model)

## Then test that PKG_CPPFLAGS is restored.
stopifnot(identical(Sys.getenv("PKG_CPPFLAGS"), "test"))

U_exp <- data.frame(
    node = c(1L, 2L, 3L, 4L, 5L, 6L, 1L, 2L, 3L, 4L, 5L, 6L, 1L, 2L, 3L, 4L,
             5L, 6L, 1L, 2L, 3L, 4L, 5L, 6L, 1L, 2L, 3L, 4L, 5L, 6L, 1L, 2L,
             3L, 4L, 5L, 6L, 1L, 2L, 3L, 4L, 5L, 6L, 1L, 2L, 3L, 4L, 5L, 6L,
             1L, 2L, 3L, 4L, 5L, 6L, 1L, 2L, 3L, 4L,  5L, 6L),
    time = c(1L, 1L, 1L, 1L, 1L, 1L, 2L, 2L, 2L, 2L, 2L, 2L, 3L, 3L, 3L, 3L,
             3L, 3L, 4L, 4L, 4L, 4L, 4L, 4L, 5L, 5L, 5L, 5L, 5L, 5L, 6L, 6L,
             6L, 6L, 6L, 6L, 7L, 7L, 7L, 7L, 7L, 7L, 8L, 8L, 8L, 8L, 8L, 8L,
             9L, 9L, 9L, 9L, 9L, 9L, 10L, 10L, 10L, 10L, 10L, 10L),
    S = c(100L, 101L, 102L, 103L, 103L, 103L, 100L, 101L,
          101L, 103L, 99L, 103L, 100L, 101L, 101L, 103L, 98L, 102L, 100L,
          101L, 101L, 103L, 98L, 100L, 100L, 100L, 100L, 102L, 98L, 99L,
          100L, 100L, 100L, 102L, 98L, 97L, 100L, 100L, 99L, 102L, 96L,
          96L, 100L, 100L, 99L, 102L, 92L, 94L, 100L, 99L, 98L, 102L, 91L,
          94L, 100L, 99L, 98L, 102L, 87L, 94L),
    I = c(1L, 1L, 3L, 4L, 6L,
          7L, 1L, 1L, 4L, 3L, 9L, 4L, 0L, 1L, 4L, 3L, 9L, 5L, 0L, 1L, 4L,
          3L, 9L, 7L, 0L, 2L, 5L, 4L, 8L, 8L, 0L, 2L, 4L, 4L, 8L, 9L, 0L,
          1L, 5L, 3L, 10L, 10L, 0L, 1L, 4L, 3L, 12L, 11L, 0L, 2L, 4L, 3L,
          13L, 11L, 0L, 2L, 4L, 2L, 14L, 10L),
    R = c(0L, 1L, 0L, 0L, 0L,
          1L, 0L, 1L, 0L, 1L, 1L, 4L, 1L, 1L, 0L, 1L, 2L, 4L, 1L, 1L, 0L,
          1L, 2L, 4L, 1L, 1L, 0L, 1L, 3L, 4L, 1L, 1L, 1L, 1L, 3L, 5L, 1L,
          2L, 1L, 2L, 3L, 5L, 1L, 2L, 2L, 2L, 5L, 6L, 1L, 2L, 3L, 2L, 5L,
          6L, 1L, 2L, 3L, 3L, 8L, 7L))

stopifnot(identical(trajectory(result), U_exp))

## Remove the C code and check that an error is raised when calling
## 'run'.
model@C_code <- character(0)
res <- assertError(run(model))
check_error(res, "The model must contain C code.")

## Check that mparse fails with invalid usage of the empty set '@'.
res <- assertError(
    mparse(transitions = c("S -> beta*S*I/(S+I+R) -> I",
                           "I -> gamma*I -> @+R"),
           compartments = c("S", "I", "R"),
           gdata = c(beta = 0.16, gamma = 0.077),
           u0 = data.frame(S = 100:105, I = 1:6, R = rep(0, 6)),
           tspan = 1:10))
check_error(res, "Invalid usage of the empty set '@'.")

res <- assertError(
    mparse(transitions = c("S -> beta*S*I/(S+I+R) -> I",
                           "I -> gamma*I -> @+@"),
           compartments = c("S", "I", "R"),
           gdata = c(beta = 0.16, gamma = 0.077),
           u0 = data.frame(S = 100:105, I = 1:6, R = rep(0, 6)),
           tspan = 1:10))
check_error(res, "Invalid usage of the empty set '@'.")

res <- assertError(
    mparse(transitions = c("S -> beta*S*I/(S+I+R) -> I",
                           "I -> gamma*I -> 2*@"),
           compartments = c("S", "I", "R"),
           gdata = c(beta = 0.16, gamma = 0.077),
           u0 = data.frame(S = 100:105, I = 1:6, R = rep(0, 6)),
           tspan = 1:10))
check_error(res, "Invalid usage of the empty set '@'.")

res <- assertError(
    SimInf:::parse_variable(x = "3N <- S + I + R",
                            compartments = c("S", "I", "R"),
                            cell_compartments = character(0),
                            ldata_names = character(0),
                            gdata_names = character(0),
                            v0_names = character(0),
                            use_enum = FALSE))
check_error(res, "Invalid variable: '3N <- S + I + R'.")

res <- assertError(
    SimInf:::parse_variable(x = "N <- S + I + R",
                            compartments = c("S", "I", "R"),
                            cell_compartments = character(0),
                            ldata_names = "N",
                            gdata_names = character(0),
                            v0_names = character(0),
                            use_enum = FALSE))
check_error(
    res,
    "Invalid variable name.")

res <- assertError(
    SimInf:::parse_variable(x = "contamination <- 2 * S",
                            compartments = c("S", "row", "col"),
                            cell_compartments = "cell.contamination",
                            ldata_names = character(0),
                            gdata_names = character(0),
                            v0_names = character(0),
                            use_enum = TRUE))
check_error(
    res,
    "Invalid variable name.")

res <- assertError(
    SimInf:::parse_variables(variables = c("N <- S + I + R",
                                           "N <- S + I + R"),
                             compartments = c("S", "I", "R"),
                             cell_compartments = character(0),
                             ldata_names = character(0),
                             gdata_names = character(0),
                             v0_names = character(0),
                             use_enum = FALSE))
check_error(
    res,
    "Variables must have non-duplicated names.")

stopifnot(identical(
    SimInf:::parse_variable(x = "N <- S + I + R",
                            compartments = c("S", "I", "R"),
                            cell_compartments = character(0),
                            ldata_names = character(0),
                            gdata_names = character(0),
                            v0_names = character(0),
                            use_enum = FALSE),
    list(variable = "N",
         tokens = c("u[0]", "+", "u[1]", "+", "u[2]"),
         type = "double",
         compartments = c("S", "I", "R"))))

stopifnot(identical(
    SimInf:::parse_variable(x = "(int)N <- S + I + R",
                            compartments = c("S", "I", "R"),
                            cell_compartments = character(0),
                            ldata_names = character(0),
                            gdata_names = character(0),
                            v0_names = character(0),
                            use_enum = FALSE),
    list(variable = "N",
         tokens = c("u[0]", "+", "u[1]", "+", "u[2]"),
         type = "int",
         compartments = c("S", "I", "R"))))

stopifnot(identical(
    SimInf:::parse_variable(x = "N <- S + I + R",
                            compartments = c("S", "I", "R"),
                            cell_compartments = character(0),
                            ldata_names = character(0),
                            gdata_names = character(0),
                            v0_names = character(0),
                            use_enum = TRUE),
    list(variable = "N",
         tokens = c("u[S]", "+", "u[I]", "+", "u[R]"),
         type = "double",
         compartments = c("S", "I", "R"))))

res <- assertError(
    SimInf:::topological_sort(
                 matrix(c(1, 0, 0, 1, 0, 0, 0, 1, 0),
                        nrow = 3, ncol = 3,
                        dimnames = list(c("A", "B", "C"),
                                        c("A", "B", "C")))))
check_error(res, "Invalid dependencies between variables.")

res <- assertError(
    SimInf:::topological_sort(
                 matrix(c(0, 0, 0, 1, 0, 0, 0, 1, 1),
                        nrow = 3, ncol = 3,
                        dimnames = list(c("A", "B", "C"),
                                        c("A", "B", "C")))))
check_error(res, "Invalid dependencies between variables.")

stopifnot(identical(
    SimInf:::topological_sort(
                 matrix(c(0, 1, 0, 1, 0, 0, 0, 0, 0),
                        nrow = 3, ncol = 3,
                        dimnames = list(c("A", "B", "C"),
                                        c("C", "B", "A")))),
    matrix(c(0, 0, 0, 1, 0, 0, 0, 1, 0),
           nrow = 3, ncol = 3,
           dimnames = list(c("A", "B", "C"),
                           c("A", "B", "C")))))

stopifnot(identical(
    SimInf:::topological_sort(
                 matrix(c(0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
                        nrow = 4, ncol = 4,
                        dimnames = list(c("A", "B", "C", "D"),
                                        c("C", "B", "A", "D")))),
    matrix(c(0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0),
           nrow = 4, ncol = 4,
           dimnames = list(c("A", "D", "B", "C"),
                           c("A", "D", "B", "C")))))

## Check to generate C code for variables.
propensity <- list(code = "N>0?beta*u[0]*u[1]/N:0",
                   depends = c(1, 1, 0), S = c(-1L, 1L, 0L),
                   G_rowname = "S -> N>0?beta*S*I/N:0 -> I",
                   variables = "N")

variables <- list(N1 = list(variable = "N1",
                            type = "double",
                            depends = character(0),
                            code = "u[0]"),
                  N2 = list(variable = "N2",
                            type = "double",
                            depends = "N1",
                            code = "N1+u[1]"),
                  N = list(variable = "N",
                            type = "double",
                           depends = "N2",
                           code = "N2+u[2]"))

stopifnot(identical(
    SimInf:::C_variables(propensity, variables),
    c("    const double N1 = u[0];",
      "    const double N2 = N1+u[1];",
      "    const double N = N2+u[2];",
      "")))

stopifnot(identical(
    SimInf:::C_enum(compartment = structure(c("S", "I", "R"),
                                            value = c(0L, 1L, 2L),
                                            n_values = 3L),
                    cell_compartments = character(0),
                    ldata_names = c("beta", "gamma", "delta",
                                    "epsilon", "zeta", "eta", "theta",
                                    "iota", "kappa", "lambda"),
                    gdata_names = structure("beta", value = 4L, n_values = 5L),
                    v0_names = structure(c("mu", "nu"),
                                         value = c(0L, 1L),
                                         n_values = 2L),
                    use_enum = TRUE),
    c("/* Enumeration constants for indicies in the 'u' vector. */",
      "enum {",
      "    S = 0,",
      "    I = 1,",
      "    R = 2,",
      "    N_COMPARTMENTS_U = 3",
      "};",
      "",
      "/* Enumeration constants for indicies in the 'v' vector. */",
      "enum {",
      "    MU = 0,",
      "    NU = 1,",
      "    N_COMPARTMENTS_V = 2",
      "};",
      "",
      "/* Enumeration constants for indicies in the 'ldata' vector. */",
      "enum {",
      "    BETA,",
      "    GAMMA,",
      "    DELTA,",
      "    EPSILON,",
      "    ZETA,",
      "    ETA,",
      "    THETA,",
      "    IOTA,",
      "    KAPPA,",
      "    LAMBDA",
      "};",
      "",
      "/* Enumeration constants for indicies in the 'gdata' vector. */",
      "enum {",
      "    BETA = 4",
      "};",
      "")))

## Check dependecies to compartments via variables.
model <- mparse(transitions = c("S -> beta*NS*NI/N -> E",
                                "E -> gamma*NE -> I",
                                "I -> delta*NI -> R",
                                "NS <- S",
                                "NE <- E",
                                "NI <- I",
                                "N <- NS+NE+NI+NR"),
                compartments = c("S", "E", "I", "R"),
                ldata = c(beta = 1, gamma = 2),
                u0 = c(S = 100, E = 0, I = 10, R = 0),
                v0 = c(delta = 3),
                tspan = 1:100)

G_expected <- new("dgCMatrix",
                  i = c(0L, 1L, 0L, 1L, 2L, 0L, 2L),
                  p = c(0L, 2L, 5L, 7L),
                  Dim = c(3L, 3L),
                  Dimnames = list(c("S -> beta*NS*NI/N -> E",
                                    "E -> gamma*NE -> I",
                                    "I -> delta*NI -> R"),
                                  c("1", "2", "3")),
                  x = c(1, 1, 1, 1, 1, 1, 1),
                  factors = list())

stopifnot(identical(model@G, G_expected))

## Check that mparse works with gdata-parameters without a name.
model  <- mparse(transitions = "S -> beta * S -> I",
                 compartments = c("S", "I"),
                 gdata = c(beta = 1, 1:4),
                 u0 = data.frame(S = 100, I = 100),
                 tspan = 1:100)

show_expected <- c(
    "Model: SimInf_model",
    "Number of nodes: 1",
    "Number of transitions: 1",
    "Number of scheduled events: 0",
    "",
    "Global data",
    "-----------",
    " Parameter Value",
    " beta      1    ",
    " Number of parameters without a name: 4",
    "",
    "Compartments",
    "------------",
    " - Empty, please run the model first")

show_observed <- capture.output(show(model))
stopifnot(identical(show_observed, show_expected))

model  <- mparse(transitions = "S -> S -> I",
                 compartments = c("S", "I"),
                 gdata = c(1:4),
                 u0 = data.frame(S = 100, I = 100),
                 tspan = 1:100)

show_expected <- c(
    "Model: SimInf_model",
    "Number of nodes: 1",
    "Number of transitions: 1",
    "Number of scheduled events: 0",
    "",
    "Global data",
    "-----------",
    " Number of parameters without a name: 4",
    "",
    "Compartments",
    "------------",
    " - Empty, please run the model first")

show_observed <- capture.output(show(model))
stopifnot(identical(show_observed, show_expected))

stopifnot(identical(
    SimInf:::C_enumeration_constants("ldata", character(0)),
    character(0)))

## Check that mparse fails for a raster model with missing row and col.
res <- assertError(
    mparse(transitions = "@->1->S",
           compartments = c("S", "cell.contamination"),
           u0 = data.frame(S = 0),
           tspan = 1:5))
check_error(res, "'row' and 'col' must exist in compartments.")

## Check that mparse fails for a raster model with duplicated
## compartments, i.e., 'cell.contamination' and 'contamination' since
## 'cell.contamination' will become 'contamination'.
res <- assertError(
    mparse(transitions = "@->1->S",
           compartments = c("S", "row", "col", "contamination",
                            "cell.contamination"),
           u0 = data.frame(S = 0, row = 0, col = 0, contamination = 0),
           tspan = 1:5))
check_error(res, "Duplicated compartment or variable name detected.")

## Check that mparse fails for a raster model with duplicated cell
## compartments.
res <- assertError(
    mparse(transitions = "@->1->S",
           compartments = c("S", "row", "col", "cell.contamination",
                            "cell.contamination"),
           u0 = data.frame(S = 0, row = 0, col = 0),
           tspan = 1:5))
check_error(res, "Duplicated compartment or variable name detected.")
