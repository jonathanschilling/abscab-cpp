#ifndef CEL_H
#define CEL_H

#include <float.h>
#include <limits.h>
#include <math.h>

namespace abscab_c {
namespace cel {

// C99 does not define M_PI
#ifndef M_PI
#define M_PI (3.14159265358979323846)
#endif

/** half of pi */
const double PI_2 = M_PI / 2.0;

/** sqrt of machine epsilon */
const double SQRT_EPS = sqrt(DBL_EPSILON);

/**
 * Compute the complete elliptic integral introduced in
 * "Numerical Calculation of Elliptic Integrals and Elliptic Functions. III"
 * by R. Bulirsch in "Numerische Mathematik" 13, 305-315 (1969):
 * cel(k_c, p, a, b) =
 * \int_0^{\pi/2} \frac{a \cos^2{\varphi} + b \sin^2{\varphi}}
 *                     {  \cos^2{\varphi} + p \sin^2{\varphi}}
 *                \frac{\mathrm{d}\varphi}
 *                     {\sqrt{\cos^2{\varphi} + k_c^2 \sin^2{\varphi}}}
 * @param k_c parameter k_c of cel(); absolute value must not be 0
 * @param p   parameter p of cel()
 * @param a   parameter a of cel()
 * @param b   parameter b of cel()
 * @return the value of cel(k_c, p, a, b)
 */
double cel(double k_c, double p, double a, double b);

} // namespace cel
} // namespace abscab_c

#endif // CEL_H
