// This file is part of TriWild, a software for generating linear/curved triangle meshes.
//
// Copyright (C) 2019 Yixin Hu <yixin.hu@nyu.edu>
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
//

#ifndef TRIWILD_RATIONAL_H
#define TRIWILD_RATIONAL_H

#ifdef FLOAT_TETWILD_USE_GMP

// ---------------------------------------------------------------------------
// GMP-backed exact rational arithmetic
// ---------------------------------------------------------------------------
#include <gmp.h>
#include <iostream>

namespace triwild {

class Rational {
public:
    mpq_t value;

    void canonicalize() { mpq_canonicalize(value); }
    int  get_sign()     { return mpq_sgn(value); }

    Rational()               { mpq_init(value); mpq_set_d(value, 0); }
    Rational(double d)       { mpq_init(value); mpq_set_d(value, d); }
    Rational(const mpq_t& v) { mpq_init(value); mpq_set(value, v); }
    Rational(const Rational& o) { mpq_init(value); mpq_set(value, o.value); }

    ~Rational() { mpq_clear(value); }

    Rational& operator=(const Rational& x) {
        if (this != &x) mpq_set(value, x.value);
        return *this;
    }
    Rational& operator=(double x) { mpq_set_d(value, x); return *this; }

    friend Rational operator-(const Rational& x)
        { Rational r; mpq_neg(r.value, x.value); return r; }
    friend Rational operator+(const Rational& x, const Rational& y)
        { Rational r; mpq_add(r.value, x.value, y.value); return r; }
    friend Rational operator-(const Rational& x, const Rational& y)
        { Rational r; mpq_sub(r.value, x.value, y.value); return r; }
    friend Rational operator*(const Rational& x, const Rational& y)
        { Rational r; mpq_mul(r.value, x.value, y.value); return r; }
    friend Rational operator/(const Rational& x, const Rational& y)
        { Rational r; mpq_div(r.value, x.value, y.value); return r; }

    friend Rational pow(const Rational& x, int p) {
        Rational r = x;
        for (int i = 1; i < std::abs(p); i++) r = r * x;
        if (p < 0) return Rational(1.0) / r;
        return r;
    }

    friend bool operator< (const Rational& r, const Rational& s) { return mpq_cmp(r.value, s.value) <  0; }
    friend bool operator> (const Rational& r, const Rational& s) { return mpq_cmp(r.value, s.value) >  0; }
    friend bool operator<=(const Rational& r, const Rational& s) { return mpq_cmp(r.value, s.value) <= 0; }
    friend bool operator>=(const Rational& r, const Rational& s) { return mpq_cmp(r.value, s.value) >= 0; }
    friend bool operator==(const Rational& r, const Rational& s) { return  mpq_equal(r.value, s.value); }
    friend bool operator!=(const Rational& r, const Rational& s) { return !mpq_equal(r.value, s.value); }

    double to_double() { return mpq_get_d(value); }

    friend std::ostream& operator<<(std::ostream& os, const Rational& r) {
        os << mpq_get_d(r.value); return os;
    }
};

} // namespace triwild

#else // FLOAT_TETWILD_USE_GMP

// ---------------------------------------------------------------------------
// Fallback: long double arithmetic
//
// Uses 80-bit extended precision (x86/x86_64 with GCC/Clang) which is
// sufficient to stabilise the AMIPS energy evaluation for elements with
// energy up to ~10^8 without requiring GMP.  The same public API as the
// GMP implementation is preserved so the rest of the codebase is unchanged.
// ---------------------------------------------------------------------------
#include <cmath>
#include <iostream>

namespace triwild {

class Rational {
public:
    long double value = 0.0L;

    void canonicalize() {}
    int  get_sign() const {
        if (value > 0.0L) return  1;
        if (value < 0.0L) return -1;
        return 0;
    }

    Rational() = default;
    Rational(double d)  : value(static_cast<long double>(d)) {}
    Rational(int i)     : value(static_cast<long double>(i)) {}
    Rational(long long i) : value(static_cast<long double>(i)) {}
    Rational(const Rational&) = default;
    ~Rational() = default;

    Rational& operator=(const Rational&) = default;
    Rational& operator=(double x) { value = static_cast<long double>(x); return *this; }

    friend Rational operator-(const Rational& x)
        { Rational r; r.value = -x.value; return r; }
    friend Rational operator+(const Rational& x, const Rational& y)
        { Rational r; r.value = x.value + y.value; return r; }
    friend Rational operator-(const Rational& x, const Rational& y)
        { Rational r; r.value = x.value - y.value; return r; }
    friend Rational operator*(const Rational& x, const Rational& y)
        { Rational r; r.value = x.value * y.value; return r; }
    friend Rational operator/(const Rational& x, const Rational& y)
        { Rational r; r.value = x.value / y.value; return r; }

    friend Rational pow(const Rational& x, int p) {
        if (p == 0) return Rational(1.0);
        Rational r = x;
        for (int i = 1; i < std::abs(p); i++) r = r * x;
        if (p < 0) { Rational one(1.0); return one / r; }
        return r;
    }

    friend bool operator< (const Rational& r, const Rational& s) { return r.value <  s.value; }
    friend bool operator> (const Rational& r, const Rational& s) { return r.value >  s.value; }
    friend bool operator<=(const Rational& r, const Rational& s) { return r.value <= s.value; }
    friend bool operator>=(const Rational& r, const Rational& s) { return r.value >= s.value; }
    friend bool operator==(const Rational& r, const Rational& s) { return r.value == s.value; }
    friend bool operator!=(const Rational& r, const Rational& s) { return r.value != s.value; }

    double to_double() const { return static_cast<double>(value); }

    friend std::ostream& operator<<(std::ostream& os, const Rational& r) {
        os << r.to_double(); return os;
    }
};

} // namespace triwild

#endif // FLOAT_TETWILD_USE_GMP

#endif // TRIWILD_RATIONAL_H
