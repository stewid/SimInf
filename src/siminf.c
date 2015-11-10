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

#include <R_ext/Rdynload.h>

#include "siminf_events.h"
#include "SISe.h"
#include "SISe3.h"

/**
* List model run functions.
*/
static const R_CallMethodDef callMethods[] =
{
    {"siminf_scheduled_events", (DL_FUNC)&siminf_scheduled_events, 5},
    {"SISe_run", (DL_FUNC)&SISe_run, 3},
    {"SISe3_run", (DL_FUNC)&SISe3_run, 3},
    {NULL, NULL, 0}
};

/**
* Register routines to R.
*
* @param info Information about the DLL being loaded
*/
void
R_init_siminf(DllInfo *info)
{
    R_registerRoutines(info, NULL, callMethods, NULL, NULL);
}
