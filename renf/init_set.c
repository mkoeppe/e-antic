/*
    Copyright (C) 2016 Vincent Delecroix

    This file is part of e-antic

    e-antic is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/


#include <e-antic/renf.h>
#include <e-antic/poly_extra.h>

void renf_init_set(renf_t dest, renf_t src)
{
    nf_init(dest->nf, src->nf->pol);
    fmpz_poly_init(dest->der);
    fmpz_poly_set(dest->der, src->der);
    dest->prec = src->prec;
    arb_init(dest->emb);
    arb_set(dest->emb, src->emb);
}
