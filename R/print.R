## This file is part of SimInf, a framework for stochastic
## disease spread simulations.
##
## Copyright (C) 2015 Pavol Bauer
## Copyright (C) 2017 -- 2019 Robin Eriksson
## Copyright (C) 2015 -- 2019 Stefan Engblom
## Copyright (C) 2015 -- 2025 Stefan Widgren
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

print_title <- function(title) {
    cat("\n")
    cat(paste0(title, "\n"))
    cat(gsub(" ", "-", sprintf("%*s\n", nchar(title), "")))
}

##' Summarise trajectory.
##' @noRd
summary_trajectory <- function(object, compartments) {
    if (is_trajectory_empty(object)) {
        cat(" - Empty, please run the model first\n")
    } else {
        qq <- lapply(compartments, function(compartment) {
            x <- as.numeric(trajectory(object, compartment, format = "matrix"))
            qq <- stats::quantile(x)
            qq <- c(qq[1L:3L], mean(x), qq[4L:5L])
        })
        qq <- do.call("rbind", qq)
        colnames(qq) <- c("Min.", "1st Qu.", "Median",
                          "Mean", "3rd Qu.", "Max.")
        rownames(qq) <- paste0(" ", compartments)
        print.table(qq, digits = 3)
    }
}

summary_output_matrix <- function(object, title, names) {
    if (length(names) > 0) {
        print_title(title)
        summary_trajectory(object, names)
    }
}

summary_matrix <- function(x) {
    qq <- t(apply(x, 1, function(xx) {
        qq <- stats::quantile(xx)
        c(qq[1L:3L], mean(xx), qq[4L:5L])
    }))
    colnames(qq) <- c("Min.", "1st Qu.", "Median",
                      "Mean", "3rd Qu.", "Max.")
    rownames(qq) <- paste0(" ", rownames(x))

    ## Summarise named parameters.
    i <- which(nchar(rownames(x)) > 0)
    if (length(i) > 0) {
        print.table(qq[i, , drop = FALSE], digits = 3)
    }

    ## Summarise parameters without a name.
    fmt <- " Number of parameters without a name: %i\n"
    i <- which(nchar(rownames(x)) == 0)
    if (length(i) > 0)
        cat(sprintf(fmt, length(i)))
    if (is.null(rownames(x)))
        cat(sprintf(fmt, nrow(x)))
}

summary_vector <- function(x) {
    ## Summarise named parameters.
    i <- which(nchar(names(x)) > 0)
    if (length(i) > 0) {
        xx <- data.frame(Parameter = names(x)[i], Value = x[i])
        print.data.frame(xx, right = FALSE, row.names = FALSE)
    }

    ## Summarise parameters without a name.
    fmt <- " Number of parameters without a name: %i\n"
    i <- which(nchar(names(x)) == 0)
    if (length(i) > 0)
        cat(sprintf(fmt, length(i)))
    if (is.null(names(x)))
        cat(sprintf(fmt, length(x)))

    ## No parameters in the vector.
    if (!length(x))
        cat(" - None\n")
}

##' Summarise model parameters
##' @noRd
summary_data_matrix <- function(x, title) {
    if (!is.null(rownames(x))) {
        print_title(title)

        ## If all columns are identical, then print the data in a
        ## 'Parameter Value' format.
        d <- sum(duplicated(x, MARGIN = 2)) + 1
        if (ncol(x) == d) {
            summary_vector(x[, 1, drop = TRUE])
        } else {
            summary_matrix(x)
        }
    }
}

##' Summarise global model parameters
##' @noRd
summary_gdata <- function(object) {
    print_title("Global data")
    summary_vector(object@gdata)
}

summary_events <- function(object) {
    print_title("Scheduled events")

    if (length(object@events@event) > 0) {
        ## Summarise exit events
        i <- which(object@events@event == 0L)
        cat(sprintf(" Exit: %i\n", length(i)))

        ## Summarise enter events
        i <- which(object@events@event == 1L)
        cat(sprintf(" Enter: %i\n", length(i)))

        ## Summarise internal transfer events
        i <- which(object@events@event == 2L)
        cat(sprintf(" Internal transfer: %i\n", length(i)))

        ## Summarise external transfer events
        i <- which(object@events@event == 3L)
        cat(sprintf(" External transfer: %i\n", length(i)))

        if (length(i) > 0) {
            ## Summarise network
            cat("\nNetwork summary\n")
            cat("---------------\n")
            id <- indegree(object)
            od <- outdegree(object)
            qq_id <- stats::quantile(id)
            qq_id <- c(qq_id[1L:3L], mean(id), qq_id[4L:5L])
            qq_od <- stats::quantile(od)
            qq_od <- c(qq_od[1L:3L], mean(od), qq_od[4L:5L])
            qq <- rbind(qq_id, qq_od)
            colnames(qq) <- c("Min.", "1st Qu.", "Median",
                              "Mean", "3rd Qu.", "Max.")
            rownames(qq) <- c(" Indegree:", " Outdegree:")
            print.table(qq, digits = 3)
        }
    } else {
        cat(" - None\n")
    }
}

##' @noRd
summary_transitions <- function(object) {
    print_title("Transitions")
    cat(paste0(" ", rownames(object@G), collapse = "\n"), sep = "\n")
}

##' Brief summary of \code{SimInf_model}
##'
##' @param object The SimInf_model \code{object}
##' @return None (invisible 'NULL').
##' @include SimInf_model.R
##' @export
##' @examples
##' ## Create an 'SIR' model with 10 nodes and initialise
##' ## it to run over 100 days.
##' model <- SIR(u0 = data.frame(S = rep(99, 10),
##'                              I = rep(1, 10),
##'                              R = rep(0, 10)),
##'              tspan = 1:100,
##'              beta = 0.16,
##'              gamma = 0.077)
##'
##' ## Brief summary of the model
##' model
##'
##' ## Run the model and save the result
##' result <- run(model)
##'
##' ## Brief summary of the result. Note that 'U' and 'V' are
##' ## non-empty after running the model.
##' result
setMethod(
    "show",
    signature(object = "SimInf_model"),
    function(object) {
        ## The model name
        cat(sprintf("Model: %s\n", as.character(class(object))))
        cat(sprintf("Number of nodes: %i\n", n_nodes(object)))
        if (n_replicates(object) > 1L)
            cat(sprintf("Number of replicates: %i\n", n_replicates(object)))
        cat(sprintf("Number of transitions: %i\n", n_transitions(object)))
        show(object@events)

        if (length(object@gdata))
            summary_gdata(object)
        if (ncol(object@ldata))
            summary_data_matrix(object@ldata, "Local data")
        summary_output_matrix(object, "Continuous state variables",
                              rownames(object@v0))
        summary_output_matrix(object, "Compartments",
                              rownames(object@S))

        invisible(object)
    }
)

##' Detailed summary of a \code{SimInf_model} object
##'
##' @param object The \code{SimInf_model} object
##' @param ... Additional arguments affecting the summary produced.
##' @return None (invisible 'NULL').
##' @include SimInf_model.R
##' @export
setMethod(
    "summary",
    signature(object = "SimInf_model"),
    function(object, ...) {
        ## The model name
        cat(sprintf("Model: %s\n", as.character(class(object))))

        ## Nodes
        cat(sprintf("Number of nodes: %i\n", n_nodes(object)))

        ## Replicates
        if (n_replicates(object) > 1L)
            cat(sprintf("Number of replicates: %i\n", n_replicates(object)))

        summary_transitions(object)
        summary_gdata(object)
        summary_data_matrix(object@ldata, "Local data")
        summary_events(object)
        summary_output_matrix(object, "Continuous state variables",
                              rownames(object@v0))
        summary_output_matrix(object, "Compartments",
                              rownames(object@S))

        invisible(NULL)
    }
)
