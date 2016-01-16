/*
 *  SimInf, a framework for stochastic disease spread simulations
 *  Copyright (C) 2015  Pavol Bauer
 *  Copyright (C) 2015 - 2016  Stefan Engblom
 *  Copyright (C) 2015 - 2016  Stefan Widgren
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include "siminf.h"

/**
 * Report error
 *
 * @param err The error code.
 */
void siminf_error(int err)
{
    switch (err) {
    case SIMINF_ERR_NEGATIVE_STATE:
        Rf_error("Negative state detected.");
        break;
    case SIMINF_ERR_ALLOC_MEMORY_BUFFER:
        Rf_error("Unable to allocate memory buffer");
        break;
    case SIMINF_INVALID_EDGE_PROBABILITY:
        Rf_error("Invalid 'p_edge': Must be in interval 0 < p_edge < 1");
        break;
    case SIMINF_INVALID_SEED_VALUE:
        Rf_error("Invalid 'seed' value");
        break;
    case SIMINF_INVALID_THREADS_VALUE:
        Rf_error("Invalid 'threads' value");
        break;
    case SIMINF_ERR_V_IS_NOT_FINITE:
        Rf_error("The continuous state 'v' is not finite.");
        break;
    case SIMINF_ERR_SAMPLE_SELECT:
        Rf_error("Unable to sample individuals for event.");
        break;
    default:
        Rf_error("Unknown error code.");
    }
}
