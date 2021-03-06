/*
    Copyright (C) 2016 Vincent Delecroix

    This file is part of e-antic

    e-antic is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include <e-antic/renf_elem.h>
#include <string.h>

void renf_elem_print_pretty(const renf_elem_t a, const char * var, const renf_t nf, slong n)
{
    char * res = renf_elem_get_str_pretty(a, var, nf, n);
    fprintf(stdout, res);
    flint_free(res);
}
