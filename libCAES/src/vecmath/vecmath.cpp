#include "vecmath.h"

#include <cmath>

/*
 * Functions borrowed from the Cephes math library
 * http://www.netlib.org/cephes/
 */
double polevl(const double x, const double coef[], const int32_t N)
{
	double ans;
	int32_t i;
	const double *p;

	p = coef;
	ans = *p++;
	i = N;

	do
		ans = ans * x  +  *p++;
	while (--i);

	return ans;
}

double p1evl(const double x, const double coef[], const int32_t N)
{
	double ans;
	int i;
	const double *p;

	p = coef;
	ans = x + *p++;
	i = N - 1;

	do
		ans = ans * x  + *p++;
	while (--i);

	return ans;
}

namespace ECHMET {
namespace CAES {

/*
 * Constants borrowed from the Cephes math library
 * http://www.netlib.org/cephes/
 */
const unsigned short VecMathCommon::PExp10[] = {
	0x2dd4,0xf306,0xfd75,0x3fa4,
	0x5934,0x74c5,0x7d94,0x4027,
	0x49e4,0x0503,0x6b7a,0x4079,
	0x4a01,0x8e13,0xb479,0x40a2,
};

const unsigned short VecMathCommon::QExp10[] = {
	0xca08,0xce51,0x45fd,0x4055,
	0x7782,0xefd6,0xe05e,0x4093,
	0xf6e2,0x650d,0x3f37,0x40a0,
};

const unsigned short VecMathCommon::PLog10[] = {
	0x974f,0x6a5f,0x09a7,0x3f08,
	0x5a1a,0xd979,0xe7ee,0x3fdf,
	0x74c9,0xc66c,0x40a2,0x401a,
	0x411d,0x7e3d,0xc9a9,0x403d,
	0xdcdc,0x64eb,0x4e6d,0x404e,
	0xd312,0x2519,0x5e12,0x404c,
	0x3130,0x89b1,0xe3a5,0x4033
};

const unsigned short VecMathCommon::QLog10[] = {
	0xd0a2,0x0dfb,0x1016,0x402e,
	0x79c7,0x47ae,0xaf6d,0x4054,
	0x455a,0xa44b,0x9542,0x406b,
	0x10d7,0x2983,0x3411,0x4073,
	0x3423,0x2a8d,0xde94,0x406a,
	0xc9c8,0x4e89,0xd578,0x404d
};

void alignedFree(void *ptr)
{
	_mm_free(ptr);
}

double VecMathCommon::cephes_frexp(const double x, int32_t *pw2) noexcept
{
union
	{
	double y;
	unsigned short sh[4];
	} u;
int32_t i;
int32_t k;
short *q;

u.y = x;
q = (short *)&u.sh[3];

/* find the exponent (power of 2) */
i  = ( *q >> 4) & 0x7ff;
if( i != 0 )
	goto ieeedon;

if( u.y == 0.0 )
	{
	*pw2 = 0;
	return( 0.0 );
	}


/* Handle denormal number. */
do
	{
	u.y *= 2.0;
	i -= 1;
	k  = ( *q >> 4) & 0x7ff;
	}
while( k == 0 );
i = i + k;

ieeedon:

i -= 0x3fe;
*pw2 = i;
*q &= 0x800f;
*q |= 0x3fe0;
return( u.y );
}

double VecMathCommon::cephes_ldexp(const double x, int32_t pw2) noexcept
{
union
	{
	double y;
	unsigned short sh[4];
	} u;
short *q;
int32_t e;

u.y = x;

q = (short *)&u.sh[3];
while( (e = (*q & 0x7ff0) >> 4) == 0 )
	{
	if( u.y == 0.0 )
		{
		return( 0.0 );
		}
/* Input is denormal. */
	if( pw2 > 0 )
		{
		u.y *= 2.0;
		pw2 -= 1;
		}
	if( pw2 < 0 )
		{
		if( pw2 < -53 )
			return(0.0);
		u.y /= 2.0;
		pw2 += 1;
		}
	if( pw2 == 0 )
		return(u.y);
	}

e += pw2;

/* Handle denormalized results */
if( e < 1 )
	{
	if( e < -53 )
		return(0.0);
	*q &= 0x800f;
	*q |= 0x10;
	/* For denormals, significant bits may be lost even
	   when dividing by 2.  Construct 2^-(1-e) so the result
	   is obtained with only one multiplication.  */
	u.y *= cephes_ldexp(1.0, e-1);
	return(u.y);
	}
else
	{
	*q &= 0x800f;
	*q |= (e & 0x7ff) << 4;
	return(u.y);
	}
}

double VecMathCommon::exp10m_single(double x) noexcept
{
	double px, xx;
	short n;

	if (std::isnan(x))
		return x;
	if (x > MAXL10)
		return DBL_INF;

	if( x < -MAXL10 )
		return 0.0;

	x = -x;

	/* Express 10**x = 10**g 2**n
	 *   = 10**g 10**( n log10(2) )
	 *   = 10**( g + n log10(2) )
	 */
	px = std::floor(LOG210 * x + 0.5 );
	n = px;
	x -= px * LG102A;
	x -= px * LG102B;

	/* rational approximation for exponential
	 * of the fractional part:
	 * 10**x = 1 + 2x P(x**2)/( Q(x**2) - P(x**2) )
	 */
	xx = x * x;
	px = x * polevl( xx, (double *)PExp10, 3 );
	x =  px/( p1evl( xx, (double *)QExp10, 3 ) - px );
	x = 1.0 + 2.0 * x;

	/* multiply by power of 2 */
	x = cephes_ldexp( x, n );

	return x;
}

double VecMathCommon::mlog10_single(double x) noexcept
{
	double z;
	double y;
	int32_t e;

	if ( std::isnan(x))
		return x;

	if (x == DBL_INF )
		return x;

	x = cephes_frexp(x, &e);

	if (x < SQRTH_LOG10) {
		e -= 1;
		x = 2 * x - 1.0; /*  2x - 1  */
	} else
		x = x - 1.0;

	/* rational form */
	z = x*x;
	y = x * ( z * polevl( x, (double *)PLog10, 6 ) / p1evl( x, (double *)QLog10, 6 ) );
	y = y - 0.5 * z;   /*  y - 0.5 * x**2  */

	/* multiply log of fraction by log10(e)
	 * and base 2 exponent by log10(2)
	 */
	z = (x + y) * L10EB_LOG10;  /* accumulate terms in order of size */
	z += y * L10EA_LOG10;
	z += x * L10EA_LOG10;
	z += e * L102B_LOG10;
	z += e * L102A_LOG10;

	return -z;
}

} // namespace CAES
} // namespace ECHMET
