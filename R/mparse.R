## This file is part of SimInf, a framework for stochastic
## disease spread simulations.
##
## Copyright (C) 2015 Pavol Bauer
## Copyright (C) 2017 -- 2019 Robin Eriksson
## Copyright (C) 2015 -- 2019 Stefan Engblom
## Copyright (C) 2015 -- 2024 Stefan Widgren
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

## Split the code in order to separate preprocessor and punctuator
## tokens from identifiers, for example:
##
## > tokenize(" bR * R ")
## [1] "bR" "*"  "R"
tokenize <- function(code) {
    ## List of valid preprocessor operator or punctuator tokens.
    operators <- c("<<=", ">>=", "!=", "%=", "##", "&&", "&=", "*=",
                   "++", "+=", "--", "-=", "->", "/=", "<<", "<=", "==",
                   ">=", ">>", "^=", "|=", "||", "!", "~", "%", "&", "(",
                   ")", "*", "+", ",", "-", "/", ":", ";", "<", "=",
                   ">", "?", "[", "]", "^", "{", "|", "}", "#")

    ## Create a matrix (1 x 2) of the code, where the first column is
    ## the token and the second column indicates if the token is one
    ## of the operators (indicated with 'op').
    code <- cbind(token = code, type = "")

    ## Iterate over each operator and try to split each row in the
    ## code in smaller pieces.
    for (op in operators) {
        code <- lapply(seq_len(nrow(code)), function(i) {
            x <- code[i, seq_len(ncol(code)), drop = FALSE]

            ## Is it a non-operator token that we could split?
            if (nchar(x[1, 2]) == 0) {
                m <- gregexpr(op, x[1, 1], fixed = TRUE)[[1]]
                if (m[1] != -1) {
                    ## The operator exists in the token. Split the
                    ## token in smaller pieces. The cut-points are
                    ## deterimined by the position and length of op
                    ## e.g. "A op B" -> "A", "op", "B".
                    x <- as.character(x[1, 1])
                    j <- 1
                    xx <- NULL
                    for (i in seq_len(length(m))) {
                        if (m[i] > j)
                            xx <- c(xx, substr(x, j, m[i] - 1))
                        j <- m[i] + attr(m, "match.length")[i]
                        xx <- c(xx, substr(x, m[i], j - 1))
                    }

                    ## Make sure last sub-string is copied.
                    if (j <= nchar(x))
                        xx <- c(xx, substr(x, j, nchar(x)))

                    ## Remove leading and trailing whitespace and drop
                    ## empty strings
                    xx <- trimws(xx)
                    xx <- xx[nchar(xx) > 0]

                    ## Create a 2-column matrix from all sub-strings
                    x <- cbind(token = xx, type = ifelse(xx == op, "op", ""))
                }
            }

            x
        })

        code <- do.call("rbind", code)
    }

    code[, 1]
}

remove_spaces <- function(x) {
    gsub(" ", "", x)
}

## Rewrite propensity
##
## Rewrite the propensity by replacing all compartments by
## \code{u[compartments[j]]} where \code{j} is the numbering in
## compartments. On return, 'depends' contains all compartments upon
## which the propensity depends.
rewrite_propensity <- function(propensity, compartments, ldata_names,
                               gdata_names, v0_names) {
    propensity <- tokenize(propensity)
    G_rowname <- paste0(propensity, collapse = "")
    depends <- integer(length(compartments))

    ## Find compartments in propensity
    i <- match(propensity, compartments)
    propensity <- ifelse(is.na(i), propensity, sprintf("u[%i]", i - 1L))
    i <- i[!is.na(i)]
    if (length(i))
        depends[i] <- 1

    ## Find ldata parameters in the propensity
    i <- match(propensity, ldata_names)
    propensity <- ifelse(is.na(i), propensity, sprintf("ldata[%i]", i - 1L))

    ## Find gdata parameters in the propensity
    i <- match(propensity, gdata_names)
    propensity <- ifelse(is.na(i), propensity, sprintf("gdata[%i]", i - 1L))

    ## Find v0 parameters in the propensity
    i <- match(propensity, v0_names)
    propensity <- ifelse(is.na(i), propensity, sprintf("v[%i]", i - 1L))

    list(propensity = paste0(propensity, collapse = ""),
         depends    = depends,
         G_rowname  = G_rowname)
}

## Generate the 'from' or 'dest' labels in the G rownames.
G_label <- function(x) {
    if (length(x) == 0)
        return("@")

    ## Prefix compartments if more than one unit, e.g., '2*S'.
    lbl <- ifelse(abs(x) > 1, paste0(abs(x), "*"), "")
    lbl <- paste0(lbl, names(x))

    ## Combine all compartments, e.g., 'S + I'
    paste0(lbl, collapse = " + ")
}

## Generate rownames from the parsed transitions
G_rownames <- function(transitions) {
    as.character(do.call("rbind", lapply(transitions, "[[", "G_rowname")))
}

parse_compartments <- function(x, compartments) {
    ## Split into 'compartment1 + compartment2 + ..'
    x <- unlist(strsplit(x, "+", fixed = TRUE))

    ## Replace 'n*compartment' with n replicates of 'compartment'
    x <- unlist(sapply(remove_spaces(x), function(xx) {
        m <- regexpr("^[[:digit:]]+[*]", xx)
        if (m != 1)
            return(xx)

        ## Determine number of replicates and remove 'n*'
        n <- regmatches(xx, m)
        xx <- sub(n, "", xx, fixed = TRUE)
        n <- as.integer(substr(n, 1, nchar(n) - 1))

        rep(xx, n)
    }))

    ## Check for valid usage of the empty set.
    if (any(x == "@") && length(x) > 1)
        stop("Invalid usage of the empty set '@'.", call. = FALSE)
    x <- x[x != "@"]

    ## Assign each compartment into its number according to the
    ## ordering in compartments
    i <- match(x, compartments)
    if (anyNA(i))
        stop(sprintf("Unknown compartment: '%s'.", x[is.na(i)]), call. = FALSE)

    tabulate(i, length(compartments))
}

parse_propensity <- function(x, compartments, ldata_names,
                             gdata_names, v0_names) {
    propensity <- remove_spaces(x[c(-1, -length(x))])
    propensity <- paste0(propensity, collapse = "->")

    ## Determine the corresponding column in the state change vector
    ## S.
    from <- parse_compartments(x[1], compartments)
    dest <- parse_compartments(x[length(x)], compartments)
    S <- dest - from

    propensity <- rewrite_propensity(propensity, compartments,
                                     ldata_names, gdata_names,
                                     v0_names)

    ## Determine the G rowname
    names(from) <- compartments
    names(dest) <- compartments
    from <- G_label(from[which(from > 0)])
    dest <- G_label(dest[which(dest > 0)])
    G_rowname <- paste(from, "->", propensity$G_rowname, "->", dest)

    list(propensity = propensity$propensity,
         depends    = propensity$depends,
         S          = S,
         G_rowname  = G_rowname)
}

parse_propensities <- function(propensities, compartments,
                               ldata_names, gdata_names, v0_names) {
    propensities <- strsplit(propensities, "->", fixed = TRUE)

    lapply(propensities, function(x) {
        if (length(x) < 3) {
            stop("Invalid transition: '",
                 paste0(x, collapse = "->"),
                 "'.",
                 call. = FALSE)
        }

        parse_propensity(x, compartments, ldata_names, gdata_names,
                         v0_names)
    })
}

pattern_variable <- function() {
    "^[[:space:]]*[a-zA-Z_][a-zA-Z_0-9]*[[:space:]]*<-"
}

parse_variable <- function(x, compartments, ldata_names, gdata_names,
                           v0_names) {
    m <- regexpr(pattern_variable(), x)
    if (m != 1)
        stop("Invalid variable: '", x, "'.", call. = FALSE)

    variable <- regmatches(x, m)
    variable <- trimws(substr(variable, 1, nchar(variable) - 2))
    if (variable %in% c(compartments, gdata_names, ldata_names,
                        v0_names)) {
        stop("Variable name already exists in 'u0', 'gdata', 'ldata' or 'v0'.",
             call. = FALSE)
    }

    x <- remove_spaces(substr(x, attr(m, "match.length") + 1, nchar(x)))

    list(variable = variable,
         tokens = tokenize(x))
}

parse_variables <- function(variables, compartments, ldata_names,
                            gdata_names, v0_names) {
    if (length(variables) == 0)
        return(list())

    variables <- lapply(variables, function(x) {
        parse_variable(x, compartments, ldata_names, gdata_names,
                       v0_names)
    })

    ## Determine variable names.
    names(variables) <- vapply(variables, "[[", character(1), "variable")
    if (any(duplicated(names(variables))))
        stop("Variables must have non-duplicated names.", call. = FALSE)

    ## Determine dependencies between variables.
    depends <- do.call("cbind", lapply(variables, function(x) {
        i <- match(x$tokens, names(variables))
        d <- integer(length(variables))
        d[i] <- 1L
        matrix(d, ncol = 1, dimnames = list(names(variables), x$variable))
    }))
    depends <- topological_sort(depends)

    lapply(variables[colnames(depends)], function(x) {
        i <- which(depends[, x$variable] > 0)
        x$depends <- colnames(depends)[i]
        x
    })
}

##' Determine if a transition should be parsed as a variable
##' @noRd
is_variable <- function(transition) {
    ## The variable name must be a valid name in C. Which means upper
    ## and lower case letters, digits, and the underscore character
    ## '_'.  Names must not begin with a digit.
    grepl(pattern_variable(), transition)
}

##' Perform a topological search of the variables using Kahn's
##' algorithm (Kahn, 1962). Kahn, A. B. (1962). Topological sorting of
##' large networks. *Communications of the ACM*, 5(11),
##' p. 558-562. \doi{10.1145/368996.369025}.
##' @noRd
topological_sort <- function(x) {
    ## First, sort lexiographically to break potential ties and get a
    ## consistent solution.
    x <- x[sort(colnames(x)), sort(colnames(x)), drop = FALSE]

    ## Find variables which have no dependencies.
    S <- colnames(x)[which(colSums(x) == 0)]

    ## Character vector that will contain the sorted variables.
    L <- character(0)

    if (length(S) == 0)
        stop("Invalid dependencies between variables.", call. = FALSE)

    m <- x
    while (length(S)) {
        var <- S[1]
        S <- S[-1]
        L <- c(L, var)
        m <- m[, -which(colnames(m) == var), drop = FALSE]
        m[var, ] <- 0L
        S <- c(S, colnames(m)[which(colSums(m) == 0)])
    }

    if (ncol(m))
        stop("Invalid dependencies between variables.", call. = FALSE)

    x[L, L, drop = FALSE]
}

parse_transitions <- function(transitions, compartments, ldata_names,
                              gdata_names, v0_names) {
    ## Determine for each transition whether it is a variable or not.
    i <- vapply(transitions, is_variable, logical(1), USE.NAMES = FALSE)

    ## Extract the variables from the transitions.
    variables <- parse_variables(transitions[i], compartments,
                                 ldata_names, gdata_names, v0_names)

    ## Extract the propensites from the transitions.
    parse_propensities(transitions[!i], compartments, ldata_names,
                       gdata_names, v0_names)
}

##' Extract variable names from data
##'
##' @param x data to extract the variable names from.
##' @param is_vector_ok TRUE if x can be a numeric vector, else FALSE.
##' @noRd
variable_names <- function(x, is_vector_ok) {
    if (is.null(x))
        return(NULL)

    if (is.data.frame(x)) {
        lbl <- colnames(x)
    } else if (isTRUE(is_vector_ok)) {
        if (is.atomic(x) && is.numeric(x)) {
            lbl <- names(x)
        } else {
            stop(paste0("'",
                        as.character(substitute(x)),
                        "' must either be a 'data.frame' ",
                        "or a 'numeric' vector."),
                 call. = FALSE)
        }
    } else if (is.matrix(x)) {
        lbl <- rownames(x)
    } else {
        stop(paste0("'",
                    as.character(substitute(x)),
                    "' must either be a 'data.frame' or a 'matrix'."),
             call. = FALSE)
    }

    if (any(duplicated(lbl)) || any(nchar(lbl) == 0)) {
        stop(paste0("'",
                    as.character(substitute(x)),
                    "' must have non-duplicated parameter names."),
             call. = FALSE)
    }

    lbl
}

## Create the state-change matrix S
state_change_matrix <- function(transitions, compartments) {
    S <- do.call("cbind", lapply(transitions, "[[", "S"))
    colnames(S) <- as.character(seq_len(dim(S)[2]))
    rownames(S) <- compartments
    S
}

## Create the dependency graph G
dependency_graph <- function(transitions, S) {
    depends <- do.call("rbind", lapply(transitions, "[[", "depends"))
    G <- ((depends %*% abs(S)) > 0) * 1
    colnames(G) <- as.character(seq_len(dim(G)[2]))
    rownames(G) <- G_rownames(transitions)
    G
}

##' Model parser to define new models to run in \code{SimInf}
##'
##' Describe your model in a logical way in R. \code{mparse} creates a
##' \code{\linkS4class{SimInf_model}} object with your model
##' definition that is ready to \code{\link{run}}.
##' @param transitions character vector containing transitions on the
##'     form \code{"X -> ... -> Y"}. The left (right) side is the
##'     initial (final) state and the propensity is written in between
##'     the \code{->}-signs. The special symbol \code{@} is reserved
##'     for the empty set. For example, \code{transitions =
##'     c("S -> k1*S*I -> I", "I -> k2*I -> R")} expresses a SIR
##'     model.
##' @param compartments contains the names of the involved
##'     compartments, for example, \code{compartments = c("S", "I",
##'     "R")}.
##' @param ldata optional data for the nodes. Can be specified either
##'     as a numeric matrix where column \code{ldata[, j]} contains
##'     the local data vector for the node \code{j} or as a
##'     \code{data.frame} with one row per node. If it's specified as
##'     a matrix, it must have row names to identify the parameters in
##'     the transitions. If it's specified as a data.frame, each
##'     column is one parameter. The local data vector is passed as an
##'     argument to the transition rate functions and the post time
##'     step function.
##' @param gdata optional data that are common to all nodes in the
##'     model. Can be specified either as a named numeric vector or as
##'     as a one-row data.frame. The names are used to identify the
##'     parameters in the transitions. The global data vector is
##'     passed as an argument to the transition rate functions and the
##'     post time step function.
##' @param u0 A \code{data.frame} (or an object that can be coerced to
##'     a \code{data.frame} with \code{as.data.frame}) with the
##'     initial state i.e. the number of individuals in each
##'     compartment in each node when the simulation starts..
##' @param v0 optional data with the initial continuous state in each
##'     node. Can be specified either as a \code{data.frame} with one
##'     row per node or as a numeric matrix where column \code{v0[,
##'     j]} contains the initial state vector for the node
##'     \code{j}. If \code{v0} is specified as a \code{data.frame},
##'     each column is one parameter. If \code{v0} is specified as a
##'     matrix, the row names identify the parameters. The 'v' vector
##'     is passed as an argument to the transition rate functions and
##'     the post time step function. The continuous state can be
##'     updated in the post time step function.
##' @template tspan-param
##' @param events A \code{data.frame} with the scheduled
##'     events. Default is \code{NULL} i.e. no scheduled events in the
##'     model.
##' @param E matrix to handle scheduled events, see
##'     \code{\linkS4class{SimInf_events}}. Default is \code{NULL}
##'     i.e. no scheduled events in the model.
##' @param N matrix to handle scheduled events, see
##'     \code{\linkS4class{SimInf_events}}. Default is \code{NULL}
##'     i.e. no scheduled events in the model.
##' @param pts_fun optional character vector with C code for the post
##'     time step function. The C code should contain only the body of
##'     the function i.e. the code between the opening and closing
##'     curly brackets.
##' @return a \code{\linkS4class{SimInf_model}} object
##' @export
##' @template mparse-example
mparse <- function(transitions = NULL, compartments = NULL, ldata = NULL,
                   gdata = NULL, u0 = NULL, v0 = NULL, tspan = NULL,
                   events = NULL, E = NULL, N = NULL, pts_fun = NULL) {
    ## Check transitions
    if (!is.atomic(transitions) ||
        !is.character(transitions) ||
        any(nchar(transitions) == 0)) {
        stop("'transitions' must be specified in a character vector.",
             call. = FALSE)
    }

    ## Check u0 and compartments
    u0 <- check_u0(u0, compartments)

    ## Extract variable names from data.
    ldata_names <- variable_names(ldata, FALSE)
    gdata_names <- variable_names(gdata, TRUE)
    v0_names <- variable_names(v0, FALSE)

    if (any(duplicated(c(compartments, gdata_names, ldata_names, v0_names)))) {
        stop("'u0', 'gdata', 'ldata' and 'v0' have names in common.",
             call. = FALSE)
    }

    ## Parse transitions
    transitions <- parse_transitions(transitions, compartments,
                                     ldata_names, gdata_names,
                                     v0_names)

    S <- state_change_matrix(transitions, compartments)
    G <- dependency_graph(transitions, S)

    SimInf_model(G      = G,
                 S      = S,
                 E      = E,
                 N      = N,
                 tspan  = tspan,
                 events = events,
                 ldata  = ldata,
                 gdata  = gdata,
                 u0     = u0,
                 v0     = v0,
                 C_code = C_code_mparse(transitions, pts_fun))
}
