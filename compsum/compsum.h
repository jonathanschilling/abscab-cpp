#ifndef COMPSUM_H
#define COMPSUM_H

#include <math.h>

/**
 * Add a single contribution to the sum.
 * The compensated sum is obtained by summing the final values of s, cs and ccs
 * after this method has been called for all contributions.
 *
 * @param contribution contribution to add to the sum
 * @param compSum[3]: {s, cs, ccs}: target for output
 */
void compAdd(double contribution, double *compSum);

#endif // COMPSUM_H
