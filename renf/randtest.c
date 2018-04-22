/*
    Copyright (C) 2016 Vincent Delecroix

    This file is part of e-antic

    e-antic is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3.0 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include <flint/fmpz_vec.h>
#include <e-antic/renf.h>
#include <e-antic/poly_extra.h>

void arb_from_interval(arb_t x, const fmpz_t c, const slong k, const slong prec)
{
    /* we build the ball that gives exactly (c 2^k, (c+1) 2^k) */
    /*  center: (2c+1) 2^(k-1) */
    /*  radius: 2^(k-1)        */
    if (prec <= 0 || prec < fmpz_bits(c) + 2)
    {
        fprintf(stderr, "negative precision");
        flint_abort();
    }

    arb_set_fmpz(x, c);
    arb_mul_2exp_si(x, x, 1);
    arb_add_si(x, x, 1, prec);
    arb_mul_2exp_si(x, x, k-1);
    arb_add_error_2exp_si(x, k-1);
}

void renf_randtest(renf_t nf, flint_rand_t state, slong len, mp_bitcnt_t bits)
{
    fmpz_poly_t p;
    fmpq_poly_t p2;
    fmpz * c_array;
    slong * k_array;
    slong n_interval, n_exact;
    slong prec;
    ulong i;
    arb_t emb;

    /* compute a random irreducible polynomial */
    if (len <= 1)
    {
        fprintf(stderr, "ERROR (renf_randtest): got length < 2\n");
        flint_abort();
    }

    fmpz_poly_init(p);
    do{
        fmpz_poly_randtest_irreducible(p, state, len, bits);
    }while(!fmpz_poly_has_real_root(p));

    /* pick a random real root */
    c_array = _fmpz_vec_init(p->length);
    k_array = malloc((p->length) * sizeof(slong));
    n_interval = 0;

    fmpz_poly_isolate_real_roots(NULL, &n_exact,
            c_array, k_array, &n_interval, p);

    if (n_interval == 0)
    {
        fprintf(stderr, "Runtime error\n");
        flint_abort();
    }

    i = n_randint(state, n_interval);

    /* construct the associated number field */
    prec = fmpz_bits(c_array+i) + 2 + n_randint(state, 100);

    arb_init(emb);
    arb_from_interval(emb, c_array+i, k_array[i], prec);

    fmpq_poly_init(p2);
    fmpq_poly_set_fmpz_poly(p2, p);

    /* renf init might not be happy with our ball !! */
    renf_init(nf, p2, emb, prec);

    _fmpz_vec_clear(c_array, p->length);
    free(k_array);
    fmpz_poly_clear(p);
    fmpq_poly_clear(p2);
    arb_clear(emb);
}
