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

#include <stdlib.h>

#include "siminf.h"
#include "siminf_vec.h"

/**
 * Free memory for integer vector
 *
 * @param vector The vector to free memory for
 * @return void
*/
void siminf_vec_free(struct siminf_vec *vector)
{
    if (vector->buf)
        free(vector->buf);
    vector->buf = NULL;
    vector->size = 0;
    vector->capacity = 0;
}

/**
 * Add value at the end
 *
 * Adds a new value at the end of the vector, after its current last
 * value. This increases the vector size by one. If and only if, the
 * new size is greater than the current capacity, the allocated memory
 * is reallocated.
 * @param vector The vector to reserve memory for
 * @param value The value to add to the vector
 * @return 0 if Ok, else error code.
 */
int siminf_vec_push_back(struct siminf_vec *vector, int value)
{
    if (vector->size == vector->capacity) {
        int err = siminf_vec_reserve(vector, 2 * vector->capacity);
        if (err)
            return err;
    }

    vector->buf[vector->size++] = value;

    return 0;
}

/**
 * Reserves memory for an integer vector
 *
 * @param vector The vector to reserve memory for
 * @param capacity The new capacity of the vector
 * @return 0 if Ok, else error code.
*/
int siminf_vec_reserve(struct siminf_vec *vector, size_t capacity)
{
    int *tmp = NULL;

    /* Allocate at least a vector of length 1 */
    if (capacity == 0 && vector->capacity == 0)
        capacity = 1;
    if (capacity <= vector->capacity)
        return 0;

    tmp = realloc(vector->buf, capacity * sizeof(int));
    if (!tmp)
        return SIMINF_ERR_ALLOC_MEMORY_BUFFER;

    vector->buf = tmp;
    vector->capacity = capacity;

    return 0;
}
