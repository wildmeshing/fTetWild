#pragma once

#include <gmp.h>
#include <mpfr.h>
#include <iostream>

namespace floatTetWild
{

class Multiprecision
{
public:
	mpfr_t value;

	static const int prec = 32;
	//https://www.mpfr.org/mpfr-current/mpfr.html#Rounding-Modes
	static const mpfr_rnd_t rnd = MPFR_RNDN;

	int get_sign()
	{
		return mpfr_sgn(value);
	}

	int get_prec_bits()
	{
		return mpfr_get_prec(value);
	}

	Multiprecision()
	{
		mpfr_init2(value, prec);
		mpfr_set_d(value, 0, rnd);
	}

	Multiprecision(double d)
	{
		mpfr_init2(value, prec);
		mpfr_set_d(value, d, rnd);
	}

	Multiprecision(const mpfr_t &v_)
	{
		mpfr_init2(value, mpfr_get_prec(v_));
		mpfr_set(value, v_, rnd);
	}

	Multiprecision(const Multiprecision &other)
	{
		mpfr_init2(value, other.prec);
		mpfr_set(value, other.value, rnd);
	}

	~Multiprecision()
	{
		mpfr_clear(value);
	}

	friend Multiprecision operator+(const Multiprecision &x, const Multiprecision &y)
	{
		Multiprecision r_out;

		mpfr_add(r_out.value, x.value, y.value, rnd);

		return r_out;
	}

	friend Multiprecision operator-(const Multiprecision &x, const Multiprecision &y)
	{
		Multiprecision r_out;

		mpfr_sub(r_out.value, x.value, y.value, rnd);

		return r_out;
	}

	friend Multiprecision operator*(const Multiprecision &x, const Multiprecision &y)
	{
		Multiprecision r_out;

		mpfr_mul(r_out.value, x.value, y.value, rnd);

		return r_out;
	}

	friend Multiprecision operator/(const Multiprecision &x, const Multiprecision &y)
	{
		Multiprecision r_out;

		mpfr_div(r_out.value, x.value, y.value, rnd);

		return r_out;
	}

	Multiprecision &operator=(const Multiprecision &x)
	{
		if (this == &x)
			return *this;

		//mpfr_init2(value, prec);

		mpfr_set(value, x.value, rnd);

		return *this;
	}

	Multiprecision &operator=(const double x)
	{
		//mpfr_init2(value, prec);

		mpfr_set_d(value, x, rnd);

		return *this;
	}

	//> < ==

	friend bool operator<(const Multiprecision &r, const Multiprecision &r1)
	{
		return mpfr_cmp(r.value, r1.value) < 0;
	}

	friend bool operator>(const Multiprecision &r, const Multiprecision &r1)
	{
		return mpfr_cmp(r.value, r1.value) > 0;
	}

	friend bool operator<=(const Multiprecision &r, const Multiprecision &r1)
	{
		return mpfr_cmp(r.value, r1.value) <= 0;
	}

	friend bool operator>=(const Multiprecision &r, const Multiprecision &r1)
	{
		return mpfr_cmp(r.value, r1.value) >= 0;
	}

	friend bool operator==(const Multiprecision &r, const Multiprecision &r1)
	{
		return mpfr_cmp(r.value, r1.value) == 0;
	}

	friend bool operator!=(const Multiprecision &r, const Multiprecision &r1)
	{
		return mpfr_cmp(r.value, r1.value) != 0;
	}

	//to double

	double to_double() const
	{
		return mpfr_get_d(value, rnd);
	}

	//<<

	friend std::ostream &operator<<(std::ostream &os, const Multiprecision &r)
	{
		os << r.to_double();

		return os;
	}

	friend Multiprecision sqrt(const Multiprecision &mp)
	{
		Multiprecision res;

		mpfr_sqrt(res.value, mp.value, rnd);

		return res;
	}

	friend Multiprecision cbrt(const Multiprecision &mp)
	{
		Multiprecision res;

		mpfr_cbrt(res.value, mp.value, rnd);

		return res;
	}
};

} // namespace fastEnvelope
