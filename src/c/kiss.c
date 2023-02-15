/*
** Minimal skeleton for using the Marsaglia KISS generator
** from Fortran 77/90/95 code.
**
** Usage (from Fortran):
**	integer kiss32, mscw, seed, state(4), ui
**	double precision uf, urcw
**
**	ui = kiss32()
**
**	uf = urcw1()
**
**	seed = mscw()
**
**  state(1) = 129281
**	state(2) = 362436069
**	state(3) = 123456789
**	state(4) = seed
**	set_kiss32(state)
**
**	get_kiss32(state)
**
** Assumptions:
**	Fortran integer == C int == 32 or more bits
**	C long int == 64 or more bits
*/

#include <limits.h>			/* for CHAR_BIT */
#include <stdint.h>			/* for uint_least32_t etc types */
#include <time.h>			/* for clock() */
#include <unistd.h>			/* for getpid() */

#define HAVE_UINT_LEAST64_T

typedef unsigned int UINT_T;
typedef uint_least32_t UINT_LEAST32_T;
typedef uint_least64_t UINT_LEAST64_T;

#define TWO_TO_32	4294967296.0

/*
** DBL_MAXRANLL_INV   = 1 / (2**64 - 1)
**		      = 2**(-64) * (1 / (1 - 2**(-64)))
**		      = 2**(-64) * (1 + 2**(-64) + (2**(-64))**2 + (2**(-64))**3 + ...)
**	             ~= 2**(-64)         -- when long double precision < 64 bits
**
** LDBL_MAXRANLL_INV  = 1 / (2**128 - 1)
**		      = 2**(-128) * (1 / (1 - 2**(-128)))
**		      = 2**(-128) * (1 + 2**(-128) + (2**(-128))**2 + (2**(-128))**3 + ...)
**	             ~= 2**(-128)         -- when long double precision < 128 bits
**
** Apple compilers refuse to compile the constant expression, so we
** rewrite it as a hexadecimal floating-point constant, which the compiler
** fortunately supports.
*/

static const double DBL_MAXRANLL_INV		= 1.0 / 18446744073709551615.0; /* 1 / (2**64 - 1) */

#define DEFAULT_C	((UINT_LEAST32_T)129281L)
#define DEFAULT_JSR	((UINT_LEAST32_T)362436069L)
#define DEFAULT_SEED	((UINT_LEAST32_T)871119182L)
#define DEFAULT_X	((UINT_LEAST32_T)123456789L)
#define DEFAULT_Y	DEFAULT_SEED

#define KISS32			kiss32_
#define RAND32()		KISS32()

static UINT_LEAST32_T c		= DEFAULT_C;
static UINT_LEAST32_T jsr	= DEFAULT_JSR;
static UINT_LEAST32_T x		= DEFAULT_X;
static UINT_LEAST32_T y		= DEFAULT_Y;

void get_kiss32_(UINT_LEAST32_T state[])
{
    state[0] = c;
    state[1] = jsr;
    state[2] = x;
    state[3] = y;
}

void set_kiss32_(UINT_LEAST32_T state[])
{
    c   = state[0];
    jsr = state[1];
    x   = state[2];
    y   = state[3];
}

UINT_LEAST32_T
KISS32(void)
{   /* result in [0, 2^{32} - 1] */
    /* Original 32-bit KISS generator, requiring 64-bit internal arithmetic */
    /* According to http://www.aptech.com/papers/rndKMi.pdf, the range of */
    /* the related KISS + Monster generator is [0, 2**(32) - 1], and period */
    /* is about 10**(8859).  Marsaglia says that the period of the KISS */
    /* generator is about 2**(123), or about 10**(37).  However, its range */
    /* is not specified in any publication that I can find.  Instrumentation */
    /* of the generator in exp/tstkis.c demonstrates that the range is  */
    /* [0,0xffffffff], which is maximal for a 32-bit generator. */

#if defined(HAVE_UINT_LEAST64_T)
    const UINT_LEAST64_T a = (UINT_LEAST64_T)333333314L;
    UINT_LEAST64_T t;
#else
    const UINT_LEAST32_T a = (UINT_LEAST32_T)333333314L;
    UINT_LEAST32_T t[2];
#endif

    UINT_LEAST32_T result;

    y = (UINT_LEAST32_T)69069L * y + (UINT_LEAST32_T)12345L;

#if MAX_UINT_LEAST32 != 0xffffffffL
    y &= (UINT_LEAST32_T)0xffffffffL;
#endif

    jsr ^= (jsr << 13);
    jsr ^= (jsr >> 17);
    jsr ^= (jsr << 5);

#if defined(HAVE_UINT_LEAST64_T)
    t = a * x + c;
    c = (t >> 32);
    x = t + c;
#else
    umul64(t, a, x);
    uadd32(t, c);
    c = t[0];
    uadd32(t, c);
    x = t[1];
#endif

    if (x < c)
    {
	x++;
	c++;
    }
    x = ~x + 1;

#if MAX_UINT_LEAST32 != 0xffffffffL
    x &= (UINT_LEAST32_T)0xffffffffL;
#endif

    result = x + y + jsr;

#if MAX_UINT_LEAST32 != 0xffffffffL
    result &= (UINT_LEAST32_T)0xffffffffL;
#endif

    return (result);
}

UINT_T
mscw_(void)
{   /* make unique generator seed */
    UINT_T n, p, seed, t;
    static UINT_T k = 0;
    static const UINT_T c = 0xfeedface;
    static const unsigned int half_wordsize =
	(unsigned int)(CHAR_BIT * sizeof(UINT_T) / 2);

    n = (UINT_T)(1 | (--k & 0x0f));

    p = (UINT_T)getpid();
    p ^= (p << half_wordsize);
    p ^= p << n;

    t = (UINT_T)time((time_t *)NULL) ^ (UINT_T)clock();
    t ^= t << half_wordsize;
    t ^= t << n;

    seed = c ^ k ^ p ^ t;

    return (seed);
}

static UINT_LEAST32_T
(lrancw)(void)
{					/* result in [0, 2^{32} - 1] */
    return (RAND32());
}

double
urcw1_(void)
{					/* result in [0,1] */
    double result;

    result = TWO_TO_32 * (double)lrancw();
    result += (double)lrancw();
    result *= DBL_MAXRANLL_INV;

    if (result > 1.0)			/* rare, but necessary */
	result = 1.0;

    return (result);
}

void
c_vurcw1_(double rvect[], int *n)
{
    for (int i = 0; i < *n; i++)
    {
        rvect[i] = urcw1_();
    }
}

void
c_vurcw2_(double rvect[], int n)
{
    for (int i = 0; i < n; i++)
    {
        rvect[i] = urcw1_();
    }
}

