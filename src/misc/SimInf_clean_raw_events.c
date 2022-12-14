/*
 * This file is part of SimInf, a framework for stochastic
 * disease spread simulations.
 *
 * Copyright (C) 2015 -- 2022 Stefan Widgren
 *
 * SimInf is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * SimInf is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <https://www.gnu.org/licenses/>.
 */

#include <R.h>
#include <Rdefines.h>
#include <R_ext/Visibility.h>
#include "SimInf.h"
#include "SimInf_openmp.h"

/**
 * Find the longest path through the events.
 *
 * @param event integer vector with the event type. Each entry must
 *        contain one of '0' (exit), '1' (enter) or '3' (external
 *        transfer event, i.e., movement).
 * @param time integer vector with the time for each event.
 * @param node integer vector with the node that the event operates
 *        on.
 * @param dest integer vector with the destination node for an
 *        external transfer event i.e.of proposals to generate Not
 *        used for the other event types.
 * @param keep integer vector for results with 1 for each event to
 *        keep, else 0.
 */
static void
SimInf_find_longest_path(
    int *event,
    int *time,
    int *node,
    int *dest,
    int *keep,
    int n)
{
    int path[n]; /* One-based indices to the events in the current
                  * search. */
    int longest_path = 0;
    int must_exit = 0;

    /* If one of the events is an exit event, then the last event in
     * the path must be an exit event. */
    for (int i = 0; i < n && must_exit == 0; i++) {
        if (event[i] == EXIT_EVENT)
            must_exit = 1;
    }

    /* Iterate over all events to identify an enter event that begins
     * each path. */
    for (int begin = 0; begin < n; begin++) {
        if (event[begin] == ENTER_EVENT) {
            int depth = 1;

            /* Clear the path. */
            memset(&path[0], 0, n * sizeof(int));

            /* Initialize the path with the first event that must be
             * an enter event. This is the root for the search. */
            path[depth - 1] = begin + 1;

            /* Check if this enter event might be the longest path,
             * for example, if there are no more events. */
            if (must_exit == 0 && longest_path == 0) {
                longest_path = 1;
                keep[path[0] - 1] = 1;
            }

            /* Perform a depth first search of the events to find the
             * longest path. */
            while (depth > 0 &&
                   depth < (n - begin) &&
                   longest_path < (n - begin))
            {
                int i = path[depth - 1] - 1;
                int from = event[i] == ENTER_EVENT ? node[i] : dest[i];

                /* Continue the search from a previous search at this
                 * depth? */
                if (path[depth] > 0) {
                    i = path[depth] - 1;
                    path[depth] = 0;
                }

                /* Find an event that is consistent with 'from' in the
                 * previous event. */
                for (int j = i + 1; j < n && path[depth] == 0; j++) {
                    if (time[j] > time[i] &&
                        from == node[j] &&
                        from != dest[j] &&
                        (event[j] == EXIT_EVENT || event[j] == EXTERNAL_TRANSFER_EVENT))
                    {
                        path[depth] = j + 1;
                        if (!(must_exit && event[j] == EXTERNAL_TRANSFER_EVENT) &&
                            (depth + 1) > longest_path)
                        {
                            longest_path = depth + 1;
                            memset(keep, 0, n * sizeof(int));
                            for (int k = 0; k < longest_path; k++)
                                keep[path[k] - 1] = 1;
                        }
                    }
                }

                if (path[depth] == 0) {
                    /* No new event found at this depth, move up in
                     * the search tree. */
                    depth -= 1;
                } else if (event[path[depth] - 1] == EXIT_EVENT) {
                    /* The last event is an exit event, move up in the
                     * search tree. */
                    path[depth] = 0;
                    depth -= 1;
                } else {
                    /* Go down in the search tree. */
                    depth += 1;
                }
            }
        }
    }
}

/**
 * Utility function to clean raw events.
 *
 * @param id integer vector with an unique identifier for each
 *        individual.
 * @param event integer vector with the event type. Each entry must
 *        contain one of '0' (exit), '1' (enter) or '3' (external
 *        transfer event, i.e., movement).
 * @param time integer vector with the time for each event.
 * @param node integer vector with the node that the event operates
 *        on.
 * @param dest integer vector with the destination node for an
 *        external transfer event i.e.of proposals to generate Not
 *        used for the other event types.
 * @return a logical vector with TRUE for each event to keep, else
 *         FALSE.
 */
SEXP attribute_hidden
SimInf_clean_raw_events(
    SEXP id,
    SEXP event,
    SEXP time,
    SEXP node,
    SEXP dest)
{
    int *ptr_id = INTEGER(id);
    int *ptr_event = INTEGER(event);
    int *ptr_time = INTEGER(time);
    int *ptr_node = INTEGER(node);
    int *ptr_dest = INTEGER(dest);
    R_xlen_t len = XLENGTH(id);
    SEXP keep;
    int *ptr_keep;

    /* Use all available threads in parallel regions. */
    SimInf_set_num_threads(-1);

    /* Check that the input vectors have an identical length > 0. */
    if (len < 1)
        Rf_error("'id' must be an integer vector with length > 0.");
    if (XLENGTH(event) != len)
        Rf_error("'event' must be an integer vector with length %i.", len);
    if (XLENGTH(time) != len)
        Rf_error("'time' must be an integer vector with length %i.", len);
    if (XLENGTH(node) != len)
        Rf_error("'node' must be an integer vector with length %i.", len);
    if (XLENGTH(dest) != len)
        Rf_error("'dest' must be an integer vector with length %i.", len);

    for (R_xlen_t i = 0; i < len; i++) {
        switch (ptr_event[i]) {
        case EXIT_EVENT:
        case ENTER_EVENT:
        case EXTERNAL_TRANSFER_EVENT:
            break;
        default:
            Rf_error("'event[%i]' is invalid.", i + 1);
        }
    }

    PROTECT(keep = Rf_allocVector(LGLSXP, len));
    ptr_keep = LOGICAL(keep);

    /* The default is to drop all events. */
    memset(ptr_keep, 0, len * sizeof(int));

    #ifdef _OPENMP
    #  pragma omp parallel num_threads(SimInf_num_threads())
    #  pragma omp single
    #endif
    for (R_xlen_t i = 0, j = 0; i < len; i++) {
        /* Check for last event or a new individual. */
        if (i == (len - 1) || ptr_id[i] != ptr_id[i + 1]) {
            #ifdef _OPENMP
            #  pragma omp task
            #endif
            SimInf_find_longest_path(
                &ptr_event[j],
                &ptr_time[j],
                &ptr_node[j],
                &ptr_dest[j],
                &ptr_keep[j],
                i - j + 1);

            j = i + 1;
        }
    }

    UNPROTECT(1);

    return keep;
}
