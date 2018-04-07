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
    istringstream is("min_poly (a2-2) embedding 1.4+/-0.1");
    renf_class NF;
    is >> NF;
    
    
    cout << NF;
    
    cout << "-----------" << endl;

    fmpq_poly_t p;
    fmpq_poly_init(p);
    fmpq_poly_set_str(p, "2  1 7");
    //renf_elem_class elem(p);   /// FIXME: Not implemented
    renf_elem_class elem(NF.get_renf());
    elem = p;
    cout << "A pretty element: " << elem << endl;
    }
    
    {
        istringstream is("min_poly (a^5-3) embedding 1+/-1");
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
        fmpq_poly_set_str(p, "1  77");
        elem = p;
        cout << "A rational pretty element: " << elem << endl;
        
        {
 
            
            renf_elem_class elem(NF.get_renf());            
            stringstream inout;   
            inout >> set_renf(NF.get_renf());

            inout << "( -a^2- 3*a2 + 5a)";
            inout >> elem;
            cout << "Wonderful " << elem << endl;
            inout << "(-1/10+1/10+a3-5*a  + 1/7 a4)";
            inout >> elem;
            cout << "Wonderful " << elem << endl;
            inout << "(5)";
            inout >> elem;
            cout << "Wonderful " << elem << endl;
            inout << "(-5/ 10  in blabla  )";
            inout >> elem;
            cout << "Wonderful " << elem << endl;
            inout << "(    a^0-1  i )";
            inout >> elem;
            cout << "Wonderful " << elem << endl;
            inout << "(a)";
            inout >> elem;
            cout << "Wonderful " << elem << endl;
            inout << "(-3  - 7 +1000)";
            inout >> elem;
            cout << "Wonderful " << elem << endl;
            inout << "(       -3-7 +1000)";
            inout >> elem;
            cout << "Wonderful " << elem << endl;
            
            inout << "(-10x/2*a^8)";
            try{
                inout >> elem;
                cout << "Done 1" << endl;
            }
            catch(const std::ios_base::failure& e){
                cout << "Caught 1" << endl;
            }
            inout << "(-12*a^-8)";
            try{
                inout >> elem;
            }
            catch(const std::ios_base::failure& e){
                cout << "Caught 2" << endl;
            }
            inout << "( 1a a)";
            try{
                inout >> elem;
                cout << "Done 3" << endl;
            }
            catch(const std::ios_base::failure& e){
                cout << "Caught 3" << endl;
            }
            inout << "( 1a^1/5 )";
            try{
                inout >> elem;
                cout << "Done 4" << endl;
            }
            catch(const std::ios_base::failure& e){
                cout << "Caught 4" << endl;
            }
            inout << "( 1*-a )";
            try{
                inout >> elem;
                cout << "Done 5" << endl;
            }
            catch(const std::ios_base::failure& e){
                cout << "Caught 5" << endl;
            }
            inout << "()";
            try{
                inout >> elem;
                cout << "Done 6" << endl;
            }
            catch(const std::ios_base::failure& e){
                cout << "Caught 6a" << endl;
            }
            inout << "(+)";
            try{
                inout >> elem;
                cout << "Done 6a" << endl;
            }
            catch(const std::ios_base::failure& e){
                cout << "Caught 6" << endl;
            }
            inout << "(3 3)";
            try{
                inout >> elem;
                cout << "Done 7" << endl;
            }
            catch(const std::ios_base::failure& e){
                cout << "Caught 7" << endl;
            }
            inout << "(3* 3)";
            try{
                inout >> elem;
                cout << "Done 7a" << endl;
            }
            catch(const std::ios_base::failure& e){
                cout << "Caught 7a" << endl;
            }
            inout << "(3*3)";
            try{
                inout >> elem;
                cout << "Done 7b" << endl;
            }
            catch(const std::ios_base::failure& e){
                cout << "Caught 7b" << endl;
            }
            inout << "(3^ 5)";
            try{
                inout >> elem;
                cout << "Done 8" << endl;
            }
            catch(const std::ios_base::failure& e){
                cout << "Caught 8" << endl;
            }
            inout << "(3^5)";
            try{
                inout >> elem;
                cout << "Done 8a" << endl;
            }
            catch(const std::ios_base::failure& e){
                cout << "Caught 8a" << endl;
            }
            inout << "(a^*)";
            try{
                inout >> elem;
                cout << "Done 9" << endl;
            }
            catch(const std::ios_base::failure& e){
                cout << "Caught 9" << endl;
            }
            inout << "(a^a)";
            try{
                inout >> elem;
                cout << "Done 9a" << endl;
            }
            catch(const std::ios_base::failure& e){
                cout << "Caught 9a" << endl;
            }

        }
        
        {
            istringstream is("min_poly (a^2-2) embedding 1+/-1");
            renf_class NF;
            is >> NF;
            
            renf_elem_class elem(NF.get_renf()); 
            
            stringstream inout;   
            inout >> set_renf(NF);
            
            inout << "(a+1)";
            inout >> elem;
            cout << "Wonderful " << elem << endl;
            
            try{
                inout << "(a^3)";
                cout << "Done 10" << endl;
            }
            catch(const std::ios_base::failure& e){
                cout << "Caught 10" << endl;
            }
            
            
        }
    }
    return 0;
}
