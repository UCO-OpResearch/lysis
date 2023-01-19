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
**      state(1) = 129281
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

#ifndef KISS_H
#define KISS_H

#include <limits.h>			/* for CHAR_BIT */
#include <stdint.h>			/* for uint_least32_t etc types */
#include <time.h>			/* for clock() */
//#include <unistd.h>			/* for getpid() */

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

// static const double DBL_MAXRANLL_INV		= 1.0 / 18446744073709551615.0; /* 1 / (2**64 - 1) */

#define DEFAULT_C	((UINT_LEAST32_T)129281L)
#define DEFAULT_JSR	((UINT_LEAST32_T)362436069L)
#define DEFAULT_SEED	((UINT_LEAST32_T)871119182L)
#define DEFAULT_X	((UINT_LEAST32_T)123456789L)
#define DEFAULT_Y	DEFAULT_SEED

#define KISS32			kiss32_
#define RAND32()		KISS32()


//static UINT_LEAST32_T c		= DEFAULT_C;
//static UINT_LEAST32_T jsr	= DEFAULT_JSR;
//static UINT_LEAST32_T x		= DEFAULT_X;
//static UINT_LEAST32_T y		= DEFAULT_Y;

void get_kiss32_(UINT_LEAST32_T state[]);

void set_kiss32_(UINT_LEAST32_T state[]);

UINT_LEAST32_T KISS32(void);

UINT_T
mscw_(void);

static UINT_LEAST32_T
(lrancw)(void);

double
urcw1_(void);

void
vurcw1_double rvect[], int n);

#endif