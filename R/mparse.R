## This file is part of SimInf, a framework for stochastic
## disease spread simulations.
##
## Copyright (C) 2015 Pavol Bauer
## Copyright (C) 2017 -- 2019 Robin Eriksson
## Copyright (C) 2015 -- 2019 Stefan Engblom
## Copyright (C) 2015 -- 2023 Stefan Widgren
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

## Split the propensity in order to separate preprocessor and
## punctuator tokens from identifiers, for example:
##
## > tokens(" bR * R ")
## [1] "bR" "*"  "R"
tokens <- function(propensity) {
    ## List of valid preprocessor operator or punctuator tokens.
    operators <- c("<<=", ">>=", "!=", "%=", "##", "&&", "&=", "*=",
                   "++", "+=", "--", "-=", "->", "/=", "<<", "<=", "==",
                   ">=", ">>", "^=", "|=", "||", "!", "~", "%", "&", "(",
                   ")", "*", "+", ",", "-", "/", ":", ";", "<", "=",
                   ">", "?", "[", "]", "^", "{", "|", "}", "#")

    ## Create a matrix (1 x 2) of the propensity, where the first
    ## column is the token and the second column indicates if the
    ## token is one of the operators (indicated with 'op').
    propensity <- cbind(token = propensity, type = "")

    ## Iterate over each operator and try to split each row in the
    ## propensity in smaller pieces.
    for (op in operators) {
        propensity <- lapply(seq_len(nrow(propensity)), function(i) {
            x <- propensity[i, seq_len(ncol(propensity)), drop = FALSE]

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

        propensity <- do.call("rbind", propensity)
    }

    propensity[, 1]
}

## Rewrite propensity
##
## Rewrite the propensity by replacing all compartments by
## \code{u[compartments[j]]} where \code{j} is the numbering in
## compartments. On return, 'depends' contains all compartments upon
## which the propensity depends.
rewrite_propensity <- function(propensity, compartments, ldata_names,
                               gdata_names, v0_names) {
    propensity <- tokens(propensity)
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

    ## Remove spaces.
    x <- gsub(" ", "", x)

    ## Replace 'n*compartment' with n replicates of 'compartment'
    x <- unlist(sapply(x, function(xx) {
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

parse_transitions <- function(transitions, compartments, ldata_names,
                              gdata_names, v0_names) {
    lapply(strsplit(transitions, "->", fixed = TRUE), function(x) {
        if (length(x) < 3) {
            stop("Invalid transition: '",
                 paste0(x, collapse = "->"),
                 "'.",
                 call. = FALSE)
        }

        ## Remove spaces
        propensity <- gsub(" ", "", x[c(-1, -length(x))])
        propensity <- paste0(propensity, collapse = "->")

        ## Determine the corresponding column in the state change
        ## vector S.
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
    })
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
        if (is.vector(x = x, mode = "numeric")) {
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
##'     c("S -> beta*S*I/(S+I+R) -> I", "I -> gamma*I -> R")}
##'     expresses the SIR model. It is also possible to define
##'     variables which can then be used in calculations of
##'     propensities or in calculations of other variables. A variable
##'     is defined by the operator \code{<-}. Using a variable for the
##'     size of the population, the SIR model can instead be written
##'     \code{transitions = c("S -> beta*S*I/N -> I",
##'     "I -> gamma*I -> R", "N <- S+I+R")}. By default, the type of a
##'     variable is defined as a double, but it is possible to also
##'     define it as an integer by writing \code{(int)} before the
##'     variable name. For example, for the SIR model, the population
##'     size can be defined as \code{"(int)N <- S+I+R"}. It is also
##'     possible to explicitly use (double) in front of the variable
##'     name, but it is not needed because it is the default. Note
##'     that the order of propensities and variables does not matter.
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
##' @template u0-param
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
