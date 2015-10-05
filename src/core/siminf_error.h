/*
 *  siminf, a framework for stochastic disease spread simulations
 *  Copyright (C) 2015  Pavol Bauer
 *  Copyright (C) 2015  Stefan Engblom
 *  Copyright (C) 2015  Stefan Widgren
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

#ifndef INCLUDE_SIMINF_ERROR_H
#define INCLUDE_SIMINF_ERROR_H

/* Error constants */
#define SIMINF_ERR_NEGATIVE_STATE          1
#define SIMINF_ERR_ALLOC_MEMORY_BUFFER     2
#define SIMINF_ERR_SPLIT_EVENTS            3
#define SIMINF_UNDEFINED_EVENT             4
#define SIMINF_INVALID_EDGE_PROBABILITY    5
#define SIMINF_INVALID_SEED_VALUE          6
#define SIMINF_INVALID_THREADS_VALUE       7

/* Definition of error function */
void siminf_error(int err);

#endif
