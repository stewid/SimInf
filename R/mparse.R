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
                    for (i in seq_along(m)) {
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

substitute_tokens <- function(tokens, pattern, replacement, use_enum) {
    i <- match(tokens, pattern)
    j <- which(!is.na(i))
    if (length(j)) {
        if (isTRUE(use_enum)) {
            lbl <- toupper(pattern[i[j]])
        } else {
            lbl <- as.character(i[j] - 1L)
        }

        tokens[j] <- sprintf("%s[%s]", replacement, lbl)
    }

    tokens
}

rewrite_tokens <- function(tokens, compartments, ldata_names,
                           gdata_names, v0_names, use_enum) {
    ## Find compartments in tokens
    tokens <- substitute_tokens(tokens, compartments, "u", use_enum)

    ## Find ldata parameters in tokens
    tokens <- substitute_tokens(tokens, ldata_names, "ldata", use_enum)

    ## Find gdata parameters in tokens
    tokens <- substitute_tokens(tokens, gdata_names, "gdata", use_enum)

    ## Find v0 parameters in tokens
    tokens <- substitute_tokens(tokens, v0_names, "v", use_enum)

    tokens
}

propensity_dependencies <- function(tokens, vars, variables,
                                    compartments) {
    depends <- integer(length(compartments))

    ## Find compartments in propensity
    i <- unique(match(tokens, compartments))
    depends[i[!is.na(i)]] <- 1

    ## Determine dependencies to compartments via variables, both
    ## direct and indirect dependencies through other variables.
    i <- unique(match(vars, names(variables)))
    while (length(i)) {
        j <- unique(match(variables[[i[1]]]$compartments, compartments))
        depends[j[!is.na(j)]] <- 1
        i <- c(i, match(variables[[i[1]]]$depends, names(variables)))
        i <- unique(i[-1])
    }

    depends
}

propensity_variables <- function(tokens, variables) {
    ## Find variables in propensity.
    i <- unique(match(tokens, names(variables)))
    variables <- names(variables)[i[!is.na(i)]]
    if (is.null(variables))
        variables <- character(0)
    variables
}

## Rewrite propensity
##
## Rewrite the propensity by replacing all compartments by
## \code{u[compartments[j]]} where \code{j} is the numbering in
## compartments. On return, 'depends' contains all compartments upon
## which the propensity depends.
rewrite_propensity <- function(propensity, variables, compartments,
                               ldata_names, gdata_names, v0_names,
                               use_enum) {
    tokens <- tokenize(propensity)
    G_rowname <- paste0(tokens, collapse = "")
    vars <- propensity_variables(tokens, variables)
    depends <- propensity_dependencies(tokens, vars, variables,
                                       compartments)
    tokens <- rewrite_tokens(tokens, compartments, ldata_names,
                             gdata_names, v0_names, use_enum)

    list(code      = paste0(tokens, collapse = ""),
         depends   = depends,
         G_rowname = G_rowname,
         variables = vars)
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

parse_propensity <- function(x, variables, compartments, ldata_names,
                             gdata_names, v0_names, use_enum) {
    propensity <- remove_spaces(x[c(-1, -length(x))])
    propensity <- paste0(propensity, collapse = "->")

    ## Determine the corresponding column in the state change vector
    ## S.
    from <- parse_compartments(x[1], compartments)
    dest <- parse_compartments(x[length(x)], compartments)
    S <- dest - from

    propensity <- rewrite_propensity(propensity, variables,
                                     compartments, ldata_names,
                                     gdata_names, v0_names, use_enum)

    ## Determine the G rowname
    names(from) <- compartments
    names(dest) <- compartments
    from <- G_label(from[which(from > 0)])
    dest <- G_label(dest[which(dest > 0)])
    G_rowname <- paste(from, "->", propensity$G_rowname, "->", dest)

    list(code       = propensity$code,
         depends    = propensity$depends,
         S          = S,
         G_rowname  = G_rowname,
         variables  = propensity$variables)
}

parse_propensities <- function(propensities, variables, compartments,
                               ldata_names, gdata_names, v0_names,
                               use_enum) {
    propensities <- strsplit(propensities, "->", fixed = TRUE)

    lapply(propensities, function(x) {
        if (length(x) < 3) {
            stop("Invalid transition: '",
                 paste0(x, collapse = "->"),
                 "'.",
                 call. = FALSE)
        }

        parse_propensity(x, variables, compartments, ldata_names,
                         gdata_names, v0_names, use_enum)
    })
}

pattern_variable <- function() {
    ## The variable name must be a valid name in C. Which means upper
    ## and lower case letters, digits, and the underscore character
    ## '_'. Names must not begin with a digit.
    paste0(c(
        "^[[:space:]]*",
        "([(]double[)]|[(]int[)])?",
        "[[:space:]]*",
        "([a-zA-Z_][a-zA-Z_0-9]*)",
        "[[:space:]]*",
        "<-"),
        collapse = "")
}

parse_variable <- function(x, compartments, ldata_names, gdata_names,
                           v0_names, use_enum) {
    m <- regexec(pattern_variable(), x)
    v <- unlist(regmatches(x, m))
    if (length(v) == 0L)
        stop("Invalid variable: '", x, "'.", call. = FALSE)

    if (startsWith(v[2], "(")) {
        type <- substr(v[2], 2, nchar(v[2]) - 1)
    } else {
        type <- "double"
    }

    variable <- v[3]
    if (variable %in% c(compartments, gdata_names, ldata_names,
                        v0_names)) {
        stop("Variable name already exists in 'u0', 'gdata', 'ldata' or 'v0'.",
             call. = FALSE)
    }

    tokens <- substr(x, attr(m[[1]], "match.length")[1] + 1, nchar(x))
    tokens <- remove_spaces(tokens)
    tokens <- tokenize(tokens)

    ## Find compartments in variable.
    i <- unique(match(tokens, compartments))
    variable_compartments <- compartments[i[!is.na(i)]]

    tokens <- rewrite_tokens(tokens, compartments, ldata_names,
                             gdata_names, v0_names, use_enum)

    list(variable = variable,
         tokens = tokens,
         type = type,
         compartments = variable_compartments)
}

parse_variables <- function(variables, compartments, ldata_names,
                            gdata_names, v0_names, use_enum) {
    if (length(variables) == 0)
        return(list())

    variables <- lapply(variables, function(x) {
        parse_variable(x, compartments, ldata_names, gdata_names,
                       v0_names, use_enum)
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
        x$code <- paste0(x$tokens, collapse = "")
        x$tokens <- NULL
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

    ## Character vector that will contain the sorted variables.
    i <- character(0)

    ## Find variables which have no dependencies.
    j <- colnames(x)[which(colSums(x) == 0)]
    if (length(j) == 0)
        stop("Invalid dependencies between variables.", call. = FALSE)

    m <- x
    while (length(j)) {
        i <- c(i, j[1])
        m <- m[, -which(colnames(m) == j[1]), drop = FALSE]
        m[j[1], ] <- 0L
        k <- colnames(m)[which(colSums(m) == 0)]
        j <- c(j, setdiff(k, j))
        j <- j[-1]
    }

    if (ncol(m))
        stop("Invalid dependencies between variables.", call. = FALSE)

    x[i, i, drop = FALSE]
}

parse_transitions <- function(transitions, compartments, ldata_names,
                              gdata_names, v0_names, use_enum) {
    ## Determine for each transition whether it is a variable or not.
    i <- vapply(transitions, is_variable, logical(1), USE.NAMES = FALSE)

    ## Extract the variables from the transitions.
    variables <- parse_variables(transitions[i], compartments,
                                 ldata_names, gdata_names, v0_names,
                                 use_enum)

    ## Extract the propensites from the transitions.
    propensities <- parse_propensities(transitions[!i], variables,
                                       compartments, ldata_names,
                                       gdata_names, v0_names,
                                       use_enum)

    list(propensities = propensities, variables = variables)
}

##' Extract variable names from data
##'
##' @param x data to extract the variable names from. Varible names
##'     can be empty, i.e., "", however, the empty variable names will
##'     be removed from the return value.
##' @param is_vector_ok TRUE if x can be a numeric vector, else FALSE.
##' @return character vector containting the variables with name and
##'     their enumeration value as attribute.
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

    ## Add enumeration value.
    value <- seq_along(lbl) - 1L
    n_values <- length(value)

    ## Keep only non-empty variable names.
    i <- which(nchar(lbl) > 0)
    lbl <- lbl[i]
    if (length(i)) {
        attr(lbl, "value") <- value[i]
        attr(lbl, "n_values") <- n_values
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
    rownames(S) <- as.character(compartments)
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

##' Model parser to define new models for \code{SimInf}
##'
##' Describe your model using a simple string syntax in
##' R. \code{mparse} parses this description, generates model-specific
##' C code, and returns a \code{\linkS4class{SimInf_model}} object
##' ready for simulation.
##'
##' @param transitions Character vector defining the state
##'     transitions.  Each element follows the format
##'     \code{"Source -> Propensity -> Target"}.
##'     \itemize{
##'       \item \strong{Syntax}: \code{"X -> rate_expr -> Y"} moves
##'         individuals from compartment \code{X} to \code{Y} with
##'         rate \code{rate_expr}.
##'       \item \strong{Empty Set}: Use \code{@} for the empty set
##'         (e.g., \code{"I -> mu*I -> @"} for death, or
##'         \code{"@ -> lambda -> S"} for birth).
##'       \item \strong{Variables}: Define intermediate variables
##'         using \code{<-}.  Example:
##'         \code{"N <- S + I + R"}. Variables can be used in
##'         propensity expressions. Order does not matter.
##'       \item \strong{Types}: By default, variables are
##'         \code{double}.  Use \code{(int)} to force integer type
##'         (e.g., \code{"(int)N <- S+I+R"}).
##'     }
##'     Example: \code{c("S -> beta*S*I/N -> I", "I -> gamma*I -> R",
##'     "N <- S+I+R")}.
##'
##' @param compartments Character vector of compartment names (e.g.,
##'     \code{c("S", "I", "R")}).
##'
##' @param ldata Optional local data (node-specific parameters). Can
##'     be:
##'     \itemize{
##'       \item A \code{data.frame} with one row per node (columns =
##'       parameters).
##'       \item A numeric matrix where columns are nodes and row names
##'       are parameters.
##'       \item A named vector (for single-node models).
##'     }
##'
##' @param gdata Optional global data (common to all nodes). Can be:
##'     \itemize{
##'       \item A named numeric vector (names identify parameters).
##'       \item A one-row \code{data.frame}.
##'     }
##'     Unnamed parameters in a vector must be referenced by index in
##'     the transitions.
##'
##' @param u0 Initial state vector. Can be a \code{data.frame},
##'     matrix, or named vector.  See
##'     \code{\linkS4class{SimInf_model}} for details.
##'
##' @param v0 optional data with the initial continuous state in each
##'     node. \code{v0} can be specified as a \code{data.frame} with
##'     one row per node, as a numeric matrix where column \code{v0[,
##'     j]} contains the initial state vector for the node \code{j},
##'     or as a named vector when the model only contains one node. If
##'     \code{v0} is specified as a \code{data.frame}, each column is
##'     one parameter. If \code{v0} is specified as a matrix, the row
##'     names identify the parameters. If \code{v0} is specified as a
##'     named vector, the names identify the parameters. The
##'     \sQuote{v} vector is passed as an argument to the transition
##'     rate functions and the post time step function. The continuous
##'     state can be updated in the post time step function.
##'
##' @template tspan-param
##'
##' @param events A \code{data.frame} with the scheduled
##'     events. Default is \code{NULL} i.e. no scheduled events in the
##'     model.
##'
##' @param E Optional select matrix for events. Can be a \code{matrix}
##'     or \code{data.frame}.  Defines which compartments are affected
##'     by events and their sampling weights.  See
##'     \code{\linkS4class{SimInf_events}} for detailed logic on
##'     \code{data.frame} columns (\code{compartment}, \code{select},
##'     \code{value}).  Default is \code{NULL} (no events).
##'
##' @param N Optional shift matrix for events. Can be a \code{matrix}
##'     or \code{data.frame}.  Defines how individuals are moved
##'     between compartments during events.  See
##'     \code{\linkS4class{SimInf_events}} for detailed logic on
##'     \code{data.frame} columns (\code{compartment}, \code{shift},
##'     \code{value}).  Default is \code{NULL} (no events).
##'
##' @param pts_fun optional character vector with C code for the post
##'     time step function. The C code should contain only the body of
##'     the function i.e. the code between the opening and closing
##'     curly brackets.
##'
##' @param use_enum Logical. If \code{TRUE}, generates C enumeration
##'     constants for parameters found in the \code{u}, \code{v},
##'     \code{gdata}, and \code{ldata} vectors to facilitate manual C
##'     code modification or use in \code{pts_fun}.  The name of each
##'     enumeration constant is the upper-case version of the
##'     corresponding parameter name (e.g., \code{beta} becomes
##'     \code{BETA}).  Additionally, the following constants are
##'     automatically generated:
##'     \itemize{
##'       \item \code{N_COMPARTMENTS_U}: The number of discrete
##'       compartments.
##'       \item \code{N_COMPARTMENTS_V}: The number of continuous
##'       state variables.
##'     }
##'     These constants facilitate indexing of the \code{u} and
##'     \code{v} vectors in C code.  Note that \code{N_COMPARTMENTS_U}
##'     and \code{N_COMPARTMENTS_V} cannot be used as compartment or
##'     variable names. Default is \code{FALSE} (parameters accessed
##'     by integer indices).
##'
##' @return a \code{\linkS4class{SimInf_model}} object
##' @export
##' @examples
##' ## For reproducibility, set the seed.
##' set.seed(22)
##'
##' ## Create an SIR model with a defined population size variable.
##' model <- mparse(
##'   transitions = c(
##'     "S -> beta * S * I / N -> I",
##'     "I -> gamma * I -> R",
##'     "N <- S + I + R"
##'   ),
##'   compartments = c("S", "I", "R"),
##'   gdata = c(beta = 0.16, gamma = 0.077),
##'   u0 = data.frame(S = 100, I = 1, R = 0),
##'   tspan = 1:100
##' )
##'
##' ## Run and plot the result.
##' result <- run(model)
##' plot(result)
mparse <- function(transitions = NULL, compartments = NULL, ldata = NULL,
                   gdata = NULL, u0 = NULL, v0 = NULL, tspan = NULL,
                   events = NULL, E = NULL, N = NULL, pts_fun = NULL,
                   use_enum = FALSE) {
    ## Check transitions
    if (!is.vector(transitions, mode = "character") ||
        any(nchar(transitions) == 0)) {
        stop("'transitions' must be specified in a character vector.",
             call. = FALSE)
    }

    ## Check u0 and compartments
    u0 <- check_u0(u0, compartments)

    ## Add enumeration value to compartments.
    attr(compartments, "value") <- seq_along(compartments) - 1L
    attr(compartments, "n_values") <- length(compartments)

    ## Extract variable names from data.
    ldata_names <- variable_names(ldata, nrow(u0) == 1L)
    gdata_names <- variable_names(gdata, TRUE)
    v0_names <- variable_names(v0, nrow(u0) == 1L)

    if (any(duplicated(c(compartments, gdata_names, ldata_names, v0_names)))) {
        stop("'u0', 'gdata', 'ldata' and 'v0' have names in common.",
             call. = FALSE)
    }

    if (any(c("N_COMPARTMENTS_U", "N_COMPARTMENTS_V") %in%
            c(compartments, gdata_names, ldata_names, v0_names))) {
        stop("Invalid compartment or variable name.",
             call. = FALSE)
    }

    ## Parse transitions
    transitions <- parse_transitions(transitions, compartments,
                                     ldata_names, gdata_names,
                                     v0_names, use_enum)

    S <- state_change_matrix(transitions$propensities, compartments)
    G <- dependency_graph(transitions$propensities, S)

    ## Generate C code.
    C_code <- C_code_mparse(transitions, pts_fun, compartments,
                            ldata_names, gdata_names, v0_names,
                            use_enum)

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
                 C_code = C_code)
}
