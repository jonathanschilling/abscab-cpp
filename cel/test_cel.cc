
#include <gtest/gtest.h>

#include "abscab-c/util/util.h"

#include "abscab-c/cel/cel.h"

/** test case for cel() implementation as described in section 4.2 of the 1969 Bulirsch article */
TEST(TestCel, CheckCel) {
	double tolerance = 1.0e-15;

	double k_c = 0.1;
	double p1 =  4.1;
	double p2 = -4.1;
	double a = 1.2;
	double b = 1.1;

	double cel1 =  1.5464442694017956;
	double cel2 = -6.7687378198360556e-1;

	double c1 = cel(k_c, p1, a, b);
	double c2 = cel(k_c, p2, a, b);

//	double ra1 = fabs(cel1 - c1)/(1.0 + fabs(cel1));
//	double ra2 = fabs(cel2 - c2)/(1.0 + fabs(cel2));
//	printf("case 1: rel/abs deviation = %g\n", ra1);
//	printf("case 2: rel/abs deviation = %g\n", ra2);

	EXPECT_EQ(assertRelAbsEquals(cel1, c1, tolerance), 0);
	EXPECT_EQ(assertRelAbsEquals(cel2, c2, tolerance), 0);
}
