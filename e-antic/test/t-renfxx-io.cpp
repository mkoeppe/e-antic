/*
    Copyright (C) 2017 Vincent Delecroix

    This file is part of e-antic

    e-antic is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include <e-antic/renfxx.h>
#include <sstream>
#include <cstdlib>

using namespace std;

int main(void)
{
    renf_elem_class f;
    cout << f << "\n";

    {
        // By default, read QQ encapsulated in renf_elem_class
        istringstream is("42");
        is >> f;
        if (f != 42) {
            cerr << "FAIL: C++ read rational" << endl;
            abort();
        }
    }

    {
        /* x^2 - 2 */ /* copied from t-init.c */
        fmpq_poly_t p;
        arb_t emb;
        renf_t nf;

        fmpq_poly_init(p);
        arb_init(emb);

        fmpq_poly_set_coeff_si(p, 0, -2);
        fmpq_poly_set_coeff_si(p, 2, 1);
        arb_set_d(emb, 1.414213562373095);
        arb_add_error_2exp_si(emb, -20);
        renf_init(nf, p, emb, 64);

        istringstream is("(4*a + 2)"); // 4 sqrt(2) + 2
        is >> set_renf(nf); // >> f;  /// COMMENTED OUT BECAUSE NOT IMPLEMENTED.
        // FIXME: Check f.
        
        renf_class NF(nf);
        cout << NF;

        renf_clear(nf);
        fmpq_poly_clear(p);
        arb_clear(emb);
    }
    {
    istringstream is("min_poly 3  -2 0 1 embedding 1.4+/-0.1");
    renf_class NF;
    is >> NF;
    
    
    cout << NF;

    fmpq_poly_t p;
    fmpq_poly_init(p);
    fmpq_poly_set_str(p, "2  1 7");
    //renf_elem_class elem(p);   /// FIXME: Not implemented
    renf_elem_class elem(NF.get_renf());
    elem = p;
    cout << "A pretty element: " << elem << endl;
    }
    
    {
        istringstream is("min_poly 6  -3 0 0 0 0 1 embedding 1+/-1");
        renf_class NF;
        is >> NF;
        
        
        cout << NF;

        fmpq_poly_t p;
        fmpq_poly_init(p);
        fmpq_poly_set_str(p, "4  1 7 3 2");
        //renf_elem_class elem(p);   /// FIXME: Not implemented
        renf_elem_class elem(NF.get_renf());
        elem = p;
        cout << "A pretty element: " << elem << endl;
    }
    return 0;
}
