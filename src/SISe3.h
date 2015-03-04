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

#ifndef INCLUDE_SISe3_h
#define INCLUDE_SISe3_h

double SISe3_S_1_to_I_1(const int *x, double t, const double *data, int sd);
double SISe3_S_2_to_I_2(const int *x, double t, const double *data, int sd);
double SISe3_S_3_to_I_3(const int *x, double t, const double *data, int sd);
double SISe3_I_1_to_S_1(const int *x, double t, const double *data, int sd);
double SISe3_I_2_to_S_2(const int *x, double t, const double *data, int sd);
double SISe3_I_3_to_S_3(const int *x, double t, const double *data, int sd);

int SISe3_post_time_step(const int *x, int src, double t, double *data);

#endif
