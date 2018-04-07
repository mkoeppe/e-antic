/*  This is a -*- C++ -*- header file.

    Copyright (C) 2016 Vincent Delecroix

    This file is part of e-antic

    e-antic is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#ifndef NF_EMB_ELEMXX_H
#define NF_EMB_ELEMXX_H

#include <iostream>
#include <string>
#include <gmpxx.h>
#include <vector>
#include <sstream>


#include "renf.h"
#include "renf_elem.h"

class renf_class
{
    renf_srcptr nf;
    renf_class(const renf_class &);
    renf_class& operator = (const renf_class&);
    renf_class& operator = (renf_srcptr a_nf) { 
        if (nf) {
            renf_clear(nf); 
            flint_free(nf);
        }
        nf = (renf_srcptr) flint_malloc(sizeof(renf_t));
        renf_init_set(nf, a_nf);
        return *this;
    }
public:
    renf_class() : nf(0) {}
    renf_class(renf_srcptr a_nf) { 
        nf = (renf_srcptr) flint_malloc(sizeof(renf_t));
        renf_init_set(nf, a_nf); 
    }
    ~renf_class() {
        if (nf) {
            renf_clear(nf); 
            flint_free(nf);
        }
    }
    renf_srcptr get_renf() const { return nf; }

    friend std::ostream& operator << (std::ostream &, const renf_class&);
    friend std::istream& operator >> (std::istream &, renf_class&);
};

class renf_elem_class
{
private:
    mutable renf_srcptr nf;
    mutable renf_elem_t a;
    mutable fmpq_t b;
public:
    // construction as zero in a given number field
    renf_elem_class(renf_t nf);

    // construction as integers or rationals
    renf_elem_class(const int=0);            // also default constructor
    renf_elem_class(const unsigned int);
    renf_elem_class(const long);
    renf_elem_class(const unsigned long);
    renf_elem_class(const mpz_class&);
    renf_elem_class(const mpq_class&);
    renf_elem_class(const fmpz_t&);
    renf_elem_class(const fmpq_t&);

    // copy constructor
    renf_elem_class(const renf_elem_class&);

    ~renf_elem_class();

    // access to attribute
    bool is_fmpq(void);
    fmpq * get_fmpq(void);
    renf_elem_srcptr get_renf_elem(void);

    // assignment
    renf_elem_class& operator = (const fmpz_t&);
    renf_elem_class& operator = (const fmpq_t&);
    renf_elem_class& operator = (const fmpq_poly_t&);
    renf_elem_class& operator = (const renf_elem_class&);

    // unary operations
    renf_elem_class operator-() const;
    renf_elem_class operator+() const;

    // binary operations
    renf_elem_class operator + (const renf_elem_class&) const;
    renf_elem_class operator - (const renf_elem_class&) const;
    renf_elem_class operator * (const renf_elem_class&) const;
    renf_elem_class operator / (const renf_elem_class&) const;
    renf_elem_class& operator += (const renf_elem_class&);
    renf_elem_class& operator -= (const renf_elem_class&);
    renf_elem_class& operator *= (const renf_elem_class&);
    renf_elem_class& operator /= (const renf_elem_class&);
    bool operator == (const renf_elem_class&) const;
    bool operator != (const renf_elem_class&) const;
    bool operator >= (const renf_elem_class&) const;
    bool operator <= (const renf_elem_class&) const;
    bool operator > (const renf_elem_class&) const;
    bool operator < (const renf_elem_class&) const;

    bool is_zero() { return renf_elem_is_zero(this->a, this->nf); };
    bool is_one() { return renf_elem_is_one(this->a, this->nf); };

    // input, output
    // I/O manipulator that stores a renf in an input stream
    // for use by operator >>.
    friend std::ios_base& set_renf(std::ios_base &, renf_t);
    friend std::ostream& operator << (std::ostream &, const renf_elem_class&);
    friend std::istream& operator >> (std::istream &, renf_elem_class&);
    void print();

    // macro for declaration of
    // - binary operations
    // - inplace binary operations
    // - comparisons
    // - assignment
    #define __renf_ops(TYP) \
    renf_elem_class operator + (const TYP) const; \
    renf_elem_class operator - (const TYP) const; \
    renf_elem_class operator * (const TYP) const; \
    renf_elem_class operator / (const TYP) const; \
    friend renf_elem_class operator + (const TYP, const renf_elem_class &); \
    friend renf_elem_class operator - (const TYP, const renf_elem_class &); \
    friend renf_elem_class operator * (const TYP, const renf_elem_class &); \
    friend renf_elem_class operator / (const TYP, const renf_elem_class &); \
    renf_elem_class& operator += (const TYP); \
    renf_elem_class& operator -= (const TYP); \
    renf_elem_class& operator *= (const TYP); \
    renf_elem_class& operator /= (const TYP); \
    bool operator == (const TYP) const; \
    bool operator != (const TYP) const; \
    bool operator >= (const TYP) const; \
    bool operator <= (const TYP) const; \
    bool operator > (const TYP) const; \
    bool operator < (const TYP) const; \
    renf_elem_class& operator = (const TYP);

    __renf_ops(mpq_class&);
    __renf_ops(mpz_class&);
    __renf_ops(unsigned long);
    __renf_ops(long);
    __renf_ops(unsigned int);
    __renf_ops(int);
    #undef __renf_ops
};

inline renf_elem_class::renf_elem_class(renf_t k)
{
    nf = k;
    renf_elem_init(a, nf);
    renf_elem_zero(a, nf);
}
inline renf_elem_class::renf_elem_class(int x)
{
    nf = NULL;
    fmpq_init(b);
    fmpq_set_si(b, x, 1);
}
inline renf_elem_class::renf_elem_class(unsigned int x)
{
    nf = NULL;
    fmpq_init(b);
    fmpz_set_ui(fmpq_numref(b), x);
    fmpz_one(fmpq_denref(b));
}
inline renf_elem_class::renf_elem_class(long x)
{
    nf = NULL;
    fmpq_init(b);
    fmpq_set_si(b, x, 1);
}
inline renf_elem_class::renf_elem_class(unsigned long x)
{
    nf = NULL;
    fmpq_init(b);
    fmpz_set_ui(fmpq_numref(b), x);
    fmpz_one(fmpq_denref(b));
}
inline renf_elem_class::renf_elem_class(const mpz_class &x)
{
    nf = NULL;
    fmpq_init(b);
    fmpz_set_mpz(fmpq_numref(b), x.__get_mp());
    fmpz_set_si(fmpq_denref(b), 1);
}
inline renf_elem_class::renf_elem_class(const mpq_class &x)
{
    nf = NULL;
    fmpq_init(b);
    fmpq_set_mpq(b, x.__get_mp());
}
inline renf_elem_class::renf_elem_class(const fmpz_t& x)
{
    nf = NULL;
    fmpq_init(b);
    fmpz_set(fmpq_numref(b), x);
    fmpz_one(fmpq_denref(b));
}
inline renf_elem_class::renf_elem_class(const fmpq_t& x)
{
    nf = NULL;
    fmpq_init(b);
    fmpq_set(b, x);
}
inline renf_elem_class::renf_elem_class(const renf_elem_class& x)
{
    nf = x.nf;
    if (nf == NULL)
    {
        fmpq_init(b);
        fmpq_set(b, x.b);
    }
    else
    {
        renf_elem_init(a, nf);
        renf_elem_set(a, x.a, nf);
    }
}


inline renf_elem_class::~renf_elem_class(void)
{
    if (nf == NULL) fmpq_clear(b);
    else renf_elem_clear(a, nf);
}


inline bool renf_elem_class::is_fmpq(void)
{
    return (nf == NULL);
}

inline fmpq * renf_elem_class::get_fmpq(void)
{
    if(not is_fmpq()) throw 42;
    else return b;
}
inline renf_elem_srcptr renf_elem_class::get_renf_elem(void)
{
    if(is_fmpq()) throw 42;
    else return a;
}


inline renf_elem_class& renf_elem_class::operator = (const int n)
{
    if (nf == NULL) fmpq_set_si(b, n, 1);
    else renf_elem_set_si(a, n, nf);
    return *this;
}
inline renf_elem_class& renf_elem_class::operator = (const unsigned int n)
{
    fmpz_t x;
    fmpz_init(x);
    fmpz_set_ui(x, n);
    *this = x;
    fmpz_clear(x);
    return *this;
}
inline renf_elem_class& renf_elem_class::operator = (const long n)
{
    if (nf == NULL) fmpq_set_si(b, n, 1);
    else renf_elem_set_si(a, n, nf);
    return *this;
}
inline renf_elem_class& renf_elem_class::operator = (const unsigned long n)
{
    fmpz_t x;
    fmpz_init(x);
    fmpz_set_ui(x, n);
    *this = x;
    fmpz_clear(x);
    return *this;
}
inline renf_elem_class& renf_elem_class::operator = (const fmpz_t& z)
{
    if (nf == NULL)
    {
        fmpz_set(fmpq_numref(b), z);
        fmpz_one(fmpq_denref(b));
    }
    else renf_elem_set_fmpz(a, z, nf);
    return *this;
}
inline renf_elem_class& renf_elem_class::operator = (const fmpq_t& q)
{
    if (nf == NULL) fmpq_set(b, q);
    else renf_elem_set_fmpq(a, q, nf);
    return *this;
}
inline renf_elem_class& renf_elem_class::operator = (const fmpq_poly_t& p)
{
    if (nf == NULL) throw 42;
    renf_elem_set_fmpq_poly(a, p, nf);
    return *this;
}
inline renf_elem_class& renf_elem_class::operator = (const mpz_class& z)
{
    fmpz_t x;
    fmpz_init(x);
    fmpz_set_mpz(x, z.get_mpz_t());
    *this = x;
    fmpz_clear(x);
    return *this;
}
inline renf_elem_class& renf_elem_class::operator = (const mpq_class& q)
{
    fmpq_t x;
    fmpq_init(x);
    fmpq_set_mpq(x, q.get_mpq_t());
    *this = x;
    fmpq_clear(x);
    return *this;
}
inline renf_elem_class& renf_elem_class::operator=(const renf_elem_class &x)
{
    if (x.nf == NULL)
    {
        if (nf != NULL)
        {
            renf_elem_clear(a, nf);
            fmpq_init(b);
            nf = NULL;
        }
        fmpq_set(b, x.b);
    }
    else if (nf == NULL)
    {
        nf = x.nf;
        fmpq_clear(b);
        renf_elem_init(a, nf);
        renf_elem_set(a, x.a, nf);
    }
    else
    {
        renf_elem_clear(a, nf);
        nf = x.nf;
        renf_elem_init(a, nf);
        renf_elem_set(a, x.a, nf);
    }

    return *this;
}

// I/O


inline std::string shorten_exact_renf_string(const std::string& long_form){

    std::string short_form;
    for(size_t i=0;i<long_form.size()-1;++i){
        if(long_form[i]=='i'){
            i+=2;
            short_form+="[";
            continue;
        }
        short_form+=long_form[i];
    }
    short_form+="]";
    return short_form;    
}

inline std::string shorten_renf_string(const std::string& long_form){

    std::string short_form;
    bool skip=false;
    bool bracket_read=false;
    for(size_t i=0;i<long_form.size();++i){
        if(long_form[i]=='i')
            skip=true;
        if(long_form[i]=='['){
            bracket_read=true;
            skip=false;
        }
        if(bracket_read && long_form[i]==' '){
            skip=true;            
        }
        if(!skip)
            short_form+=long_form[i];        
    }
    if(!bracket_read)
        return shorten_exact_renf_string(long_form);
    short_form+="?]";
    return short_form;    
}

inline std::ostream& operator<<(std::ostream & os, const renf_elem_class& a)
{
    char * res;
    if (a.nf == NULL) {
        if(fmpz_is_one(fmpq_denref(a.b))){
            res = fmpq_get_str(NULL, 10, a.b);
            os << res;
        }
        else{
            arb_t emb;
            char * res1;
            arb_init(emb);
            arb_set_fmpq(emb,a.b,23);
            res1=arb_get_str(emb,5,0);
            res = fmpq_get_str(NULL, 10, a.b);
            std::string total=res;
            total+=" in ";
            total+=+res1;
            std::string short_output=shorten_renf_string(total);
            // os << "(" << res << " in " << res1 << ")";
            os << "(" << short_output << ")";
            flint_free(res1);
        }
    }
    else {
        
        if(renf_elem_is_integer(a.a,a.nf)){
            res = nf_elem_get_str_pretty(a.a->elem, "a", a.nf->nf);
            os  << res;
        }
        else{
            if(true){ // renf_elem_output_short
                renf_elem_set_evaluation(a.a,a.nf,23);
                res = renf_elem_get_str_pretty(a.a, "a", a.nf, 5);
                std::string short_output=shorten_renf_string(res);
                os << "(" << short_output << ")";
            }
            else{
                renf_elem_set_evaluation(a.a,a.nf,23);
                res = renf_elem_get_str_pretty(a.a, "a", a.nf, 5);
                os << "(" << res << ")";
                }
        }
    }
    flint_free(res);
    return os;
}

inline std::ostream& operator<<(std::ostream & os, const renf_class& nf)
{
    char *res, *res1;
    res=fmpq_poly_get_str_pretty(nf.nf->nf->pol,"a");
    res1=arb_get_str(nf.nf->emb,64,0);
    os << "min_poly "<< "(" << res << ")" << " embedding " << res1 << std::endl;
    flint_free(res);
    flint_free(res1);
    return os;
}

inline std::vector<mpq_class> terms_to_vector(std::vector<std::string> term_strings){
    
    std::vector<mpq_class> result;
    
    for(size_t i=0;i<term_strings.size();++i){
        
        // std::cout << "Doing " << term_strings[i] << std::endl;       

        bool has_content=false;
        bool digit_then_spaces=false;
        
        std::string purified;
        
        for(size_t j=0;j<term_strings[i].size();++j){
            char test=term_strings[i][j];
            if(test=='a' || isdigit(test))
                has_content=true;
            if(isspace(test))
                continue;
            if(!isdigit(test))
                digit_then_spaces=false;
            else{
                if(digit_then_spaces)
                    throw std::ios_base::failure("Error in reading number field element: space separates digits");
                if(j<term_strings[i].size()-1 && isspace(term_strings[i][j+1]))
                    digit_then_spaces=true;                
            }            
            purified+=test;               
        }
        if(!has_content)
            throw std::ios_base::failure("Error in reading number field element: empty term or illegal character in it");

        bool a_read=false;
        bool a_just_read=false;
        bool caret_read=false;
        // bool caret_just_read=false;
        bool star_read=false;
        // bool star_just_read=false;
        bool last_read_digit=false;
        
        std::string coeff_string, exp_string;
        mpq_class coeff=1;
        long sign=1;
        int expo=0;
        
        for(size_t j=0;j<purified.size();++j){
            char test=purified[j];
           
            if(test=='a'){
                if(a_read){
                    // std::cout << "Double a" << std::endl;
                    throw std::ios_base::failure("Error in reading number field element: double a");
                }
                a_read=true;
                a_just_read=true;
                expo=1;
                continue;
            }
            
            if(test=='^'){
                if(!a_just_read || caret_read || j==purified.size()-1 || !isdigit(purified[j+1]) )
                    throw std::ios_base::failure("Error in reading number field element: double ^ or not between a and expo");
                caret_read=true;
                continue;
            }
            
            if(test=='*'){
                if(star_read || !last_read_digit || j==purified.size()-1 || purified[j+1]!='a')
                    throw std::ios_base::failure("Error in reading number field element: double * or * not between coeff and a");
                star_read=true;
                continue;
            }
                    
            if(!a_read){
                if(test=='+') // no leading + allowed for mpq_class
                    continue;
                if(test=='-'){
                    sign=-1;
                    continue;
                }                    
                coeff_string+=test;
            }
            else
                exp_string+=test;
            
            if(test!='a')
                a_just_read=false;
            
            if(isdigit(test))
                last_read_digit=true;
            else
                last_read_digit=false;
            
        }
                
        // std::cout << "Coeff_string " << coeff_string << " Exp_string " << exp_string << std::endl;
        if(coeff_string.size()>0){
            try{
                coeff=mpq_class(coeff_string);
            }
            catch(const std::invalid_argument& e) {
                throw std::ios_base::failure("Error in reading number field element: invalid coefficient "+coeff_string);                
            }            
            coeff=mpq_class(coeff_string);
        }
        coeff*=sign;
        coeff.canonicalize();
        mpz_class mpz_expo;
        if(exp_string.size()>0){
            try{
            mpz_expo=mpz_class(exp_string);
            }
            catch(const std::invalid_argument& e) {
                throw std::ios_base::failure("Error in reading number field element: invalid exponent  "+exp_string);                
            }
            expo=mpz_expo.get_si();
        }
        if(expo<0)
            throw std::ios_base::failure("Error in reading number field element: negative exponent");
        
        if(result.size()<=expo)
            result.resize(expo+1);
        result[expo]+=coeff;
        // std::cout << "Coeff " << coeff << " Exp " << expo << std::endl;
    }
    
    return result;
}

inline std::vector<mpq_class> poly_components(std::string poly_string){
    
    bool start=true;
    std::vector<std::string> term_strings;
    std::string current;
    std::vector<mpq_class> result;
    
    // std::cout << "Given sting " << poly_string << " " << std::endl;
    
    for(size_t i=0;i<poly_string.size();++i){
        if(poly_string[i]=='+'|| poly_string[i]=='-'){
            if(!start){
                term_strings.push_back(current);
                current.clear();
            }
        }
        if(!isspace(poly_string[i]))
            start=false;
        current+=poly_string[i];
    }
    term_strings.push_back(current);    
    
    /* for(size_t i=0; i<term_strings.size();++i)
        std::cout << "i " << i << " ---- " << term_strings[i] << std::endl;*/        
        
    return terms_to_vector(term_strings);
}

inline void vector2fmpq_poly(fmpq_poly_t flp, const std::vector<mpq_class>& poly_vector){
    
    slong n= (slong) poly_vector.size();

    fmpq_poly_fit_length(flp,n);
    for(size_t i=0;i<poly_vector.size();++i){
        fmpq_poly_set_coeff_mpq(flp,(slong) i, poly_vector[i].get_mpq_t());
    }

}

inline std::istream& operator>>(std::istream & is, renf_class& a)
{
    char c;
    std::string s;
    is >> s;   
    if(s!="min_poly")
        throw std::ios_base::failure("Error in reading number field: expected keyword min_poly");
    is >> std::ws;
    c=is.peek();
    if(c!='(')
        throw std::ios_base::failure("Error in reading number field: min_poly does not start with (");
    is >> c;
    
    std::string mp_string;
    while(true){
        c=is.peek();
        if(c==')'){
            is.get(c);
            break;
        }
        is.get(c);
        if(is.fail())
            throw std::ios_base::failure("Error in reading number field: min_poly not terminated by )");
        mp_string+=c;               
    }
    std::vector<mpq_class> mp_vector=poly_components(mp_string);
    fmpq_poly_t inpoly;
    fmpq_poly_init(inpoly);
    vector2fmpq_poly(inpoly,mp_vector);
    /* int error = fmpq_poly_set_str(inpoly,t.c_str());
    if (error)
         throw std::ios_base::failure("Error in reading number field: invalid polynomial " + t);*/
    
    is >> s;
    if(s!="embedding")
        throw std::ios_base::failure("Error in reading number field: expected keyword embedding");
    is >> std::ws;
    std::string u;
    c=is.peek();
    if(c=='['){
        while(true){
            is.get(c);
            u+=c;
            if(c==']')
                break;
        }
    }
    else{
        is >> u;
    }
    arb_t emb;
    arb_init(emb);
    int error=arb_set_str(emb,u.c_str(),10);
    if(error)
        throw std::ios_base::failure("Error in reading number field: bad formatting of embedding " + u );
    renf_t nf;
    renf_init(nf,inpoly,emb,64);
    a=nf; 
    return is;
}

struct set_renf {
    renf_srcptr _nf;  // Does not belong to us.
    set_renf(renf_srcptr nf) { _nf = nf; }
    set_renf(const renf_class &NF) { _nf = NF.get_renf(); }
    static int xalloc();
};

inline int set_renf::xalloc()
{
    static int xa = std::ios_base::xalloc();
    return xa;
}

inline std::istream& operator>>(std::istream & is, const set_renf &sr)
{
    is.iword(set_renf::xalloc()) = (long) sr._nf;
    return is;
}



inline std::istream& operator>>(std::istream & is, renf_elem_class& a)
{
    renf *nf = (renf *) is.iword(set_renf::xalloc());
    
    if (!nf) {
        // If no number field has been set, use rational input.
        mpq_class x;
        is >> x;
        a = x;
    }
    else {
        long degree=fmpq_poly_degree(nf->nf->pol);
        
        is >> std::ws;
        char c=is.peek();
        if(c!='('){
            mpq_class x;
            is >> x;
            a = x;
        }
        else{
            bool error_par=false;
            bool skip=false;
            is.get(c);
            std::string poly_string;
            while(true){
                is.get(c);
                if(c==')')
                    break;
                if(!is.good())
                    throw std::ios_base::failure("Error in reading number field element: unexpected end of input");
                if(c=='(')
                    error_par=true;
                if(c=='i' || c=='[')
                    skip=true;
                if(!skip)
                    poly_string+=c;
            }
            
            if(error_par)
                throw std::ios_base::failure("Error in reading number field element: double )");              
                
            std::vector<mpq_class> poly_vector=poly_components(poly_string);
            if(poly_vector.size()>= degree+1)
                throw std::ios_base::failure("Error in reading number field element: nonreduced element read");  
            
            fmpq_poly_t flp;
            fmpq_poly_init(flp);
            vector2fmpq_poly(flp,poly_vector);
            
            // Set up element in correct number field
            renf_elem_class a1(nf);           
            a1=flp;

            // Fill result
            a=a1;
        }
    }
    return is;
}


inline renf_elem_class renf_elem_class::operator-() const
{
    renf_elem_class ans(*this);
    if (nf == NULL) fmpq_neg(ans.b, ans.b);
    else renf_elem_neg(ans.a, ans.a, ans.nf);
    return ans;
}
inline renf_elem_class renf_elem_class::operator+() const
{
    return *this;
}


static inline void renf_elem_fmpq_add(renf_elem_t a, fmpq_t b, renf_elem_t c, renf_t d) {renf_elem_add_fmpq(a, c, b, d);}
static inline void renf_elem_fmpq_mul(renf_elem_t a, fmpq_t b, renf_elem_t c, renf_t d) {renf_elem_mul_fmpq(a, c, b, d);}

#define __renf_elem_op(OP, INOP, FUN1, FUN2, FUN3, FUN4) \
inline renf_elem_class& renf_elem_class::operator INOP (const renf_elem_class & other) \
{                                         \
    if (nf != NULL)                       \
    {                                     \
        if (nf == other.nf)               \
            FUN1(a, a, other.a, nf);      \
        else if (other.nf == NULL)        \
            FUN2(a, a, other.b, nf);      \
        else                              \
            throw 42;                     \
    }                                     \
    else if (other.nf != NULL)            \
    {                                     \
        /* promote to nf elt */           \
        nf = other.nf;                    \
        renf_elem_init(a, nf);            \
        renf_elem_set_fmpq(a, b, nf);     \
        fmpq_clear(b);                    \
        FUN1(a, a, other.a, nf);          \
    }                                     \
    else                                  \
    {                                     \
        /* all rationals */               \
        FUN4(b, b, other.b);              \
    }                                     \
    return *this;                         \
}                                         \
inline renf_elem_class renf_elem_class::operator OP (const renf_elem_class & other) const \
{                                         \
    renf_elem_class ans(*this);           \
    ans INOP other;                       \
    return ans;                           \
}
__renf_elem_op(+, +=, renf_elem_add, renf_elem_add_fmpq, renf_elem_fmpq_add, fmpq_add);
__renf_elem_op(*, *=, renf_elem_mul, renf_elem_mul_fmpq, renf_elem_fmpq_mul, fmpq_mul);
__renf_elem_op(-, -=, renf_elem_sub, renf_elem_sub_fmpq, renf_elem_fmpq_sub, fmpq_sub);
__renf_elem_op(/, /=, renf_elem_div, renf_elem_div_fmpq, renf_elem_fmpq_div, fmpq_div);
#undef __renf_elem_op


#define __renf_elem_op(TYP, OP, INOP) \
inline renf_elem_class& renf_elem_class::operator INOP (const TYP other) \
{                                         \
    renf_elem_class x(other);             \
    *this INOP x;                         \
    return *this;                         \
}                                         \
inline renf_elem_class renf_elem_class::operator OP (const TYP other) const \
{                                         \
    renf_elem_class x(*this);             \
    x INOP other;                         \
    return x;                             \
}                                         \
inline renf_elem_class operator OP (const TYP a, const renf_elem_class& b) \
{                                         \
    renf_elem_class x(a);                 \
    x INOP b;                             \
    return x;                             \
}

__renf_elem_op(int, +, +=);
__renf_elem_op(unsigned int, +, +=);
__renf_elem_op(long, +, +=);
__renf_elem_op(unsigned long, +, +=);
__renf_elem_op(mpz_class&, +, +=);
__renf_elem_op(mpq_class&, +, +=);


__renf_elem_op(int, -, -=);
__renf_elem_op(unsigned int, -, -=);
__renf_elem_op(long, -, -=);
__renf_elem_op(unsigned long, -, -=);
__renf_elem_op(mpz_class&, -, -=);
__renf_elem_op(mpq_class&, -, -=);

__renf_elem_op(int, *, *=);
__renf_elem_op(unsigned int, *, *=);
__renf_elem_op(long, *, *=);
__renf_elem_op(unsigned long, *, *=);
__renf_elem_op(mpz_class&, *, *=);
__renf_elem_op(mpq_class&, *, *=);

__renf_elem_op(int, /, /=);
__renf_elem_op(unsigned int, /, /=);
__renf_elem_op(long, /, /=);
__renf_elem_op(unsigned long, /, /=);
__renf_elem_op(mpz_class&, /, /=);
__renf_elem_op(mpq_class&, /, /=);

#undef __renf_elem_op




inline bool renf_elem_class::operator==(const renf_elem_class& other) const
{
    if (this->nf != NULL)
    {
        if (this->nf == other.nf)
            return renf_elem_equal(this->a, other.a, this->nf);
        else if (other.nf == NULL)
            return renf_elem_cmp_fmpq(this->a, other.b, this->nf) == 0;
        else
            throw 42;
    }
    else if (other.nf == NULL)
    {
        return fmpq_equal(this->b, other.b);
    }
    else
    {
        return renf_elem_cmp_fmpq(other.a, this->b, other.nf) == 0;
    }
}

inline bool renf_elem_class::operator>(const renf_elem_class & other) const
{
    if (this->nf != NULL)
    {
        if (this->nf == other.nf)
            return renf_elem_cmp(this->a, other.a, this->nf) > 0;
        else if (other.nf == NULL)
            return renf_elem_cmp_fmpq(this->a, other.b, this->nf) > 0;
        else
            throw 42;
    }
    else if (other.nf == NULL)
    {
        return fmpq_cmp(this->b, other.b) > 0;
    }
    else
    {
        return renf_elem_cmp_fmpq(other.a, this->b, other.nf) < 0;
    }
}

#define __other_ops(TYP) \
inline bool renf_elem_class::operator == (const TYP other) const \
{                               \
    renf_elem_class x(other);   \
    return (*this) == x;        \
}                               \
inline bool renf_elem_class::operator > (const TYP other) const \
{                               \
    renf_elem_class x(other);   \
    return (*this) > x;         \
}
__other_ops(int);
__other_ops(unsigned int);
__other_ops(long);
__other_ops(unsigned long);
__other_ops(mpz_class&);
__other_ops(mpq_class&);
#undef __other_ops

#define __other_ops(TYP) \
inline bool renf_elem_class::operator!=(const TYP other) const {return not (*this == other);} \
inline bool renf_elem_class::operator>=(const TYP other) const {return *this == other || *this > other;} \
inline bool renf_elem_class::operator<=(const TYP other) const {return not (*this > other);} \
inline bool renf_elem_class::operator< (const TYP other) const {return not (*this >= other);}

__other_ops(int);
__other_ops(unsigned int);
__other_ops(long);
__other_ops(unsigned long);
__other_ops(renf_elem_class&);
__other_ops(mpz_class&);
__other_ops(mpq_class&);
#undef __other_ops

#endif
