/*
    Copyright (C) 2016 Vincent Delecroix

    This file is part of e-antic

    e-antic is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/


#include "renf_elem.h"

void renf_elem_floor(fmpz_t a, renf_elem_t b, renf_t nf)
{
    arf_t cl, cr;
    slong prec;

#ifdef DEBUG
    printf("[renf_elem_floor]: nf with pol "); fmpq_poly_print_pretty(nf->nf->pol, "x"); printf("\n");
    printf("[renf_elem_floor]: embedding "); arb_printd(nf->emb, 10); printf("\n");
    printf("[renf_elem_floor]: b = "); renf_elem_print_pretty(b, nf, "a", 10); printf("\n");
#endif

    arf_init(cl);
    arf_init(cr);
    prec = nf->prec;
    do{
        arb_get_interval_arf(cl, cr, b->emb, prec);
#ifdef DEBUG
        printf("[floor] cl = "); arf_printd(cl, 30); printf("\n");
        printf("[floor] cr = "); arf_printd(cr, 30); printf("\n");
#endif
        arf_floor(cl, cl);
        arf_floor(cr, cr);
#ifdef DEBUG
        printf("[floor] floor(cl) = "); arf_printd(cl, 30); printf("\n");
        printf("[floor] floor(cr) = "); arf_printd(cr, 30); printf("\n");
#endif
        if(arf_equal(cl,cr))
        {
            arf_get_fmpz(a, cl, ARF_RND_NEAR);
            arf_clear(cl);
            arf_clear(cr);
            return;
        }
        prec *= 2;
        if(arf_bits(arb_midref(nf->emb)) < prec)
            renf_refine_embedding(nf, 2 * prec);
        if(2 * arf_bits(arb_midref(b->emb)) < prec)
            renf_elem_set_evaluation(b, nf, prec);
    }while (1);
}
