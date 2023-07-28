
#include <gtest/gtest.h>

#include "abscab-cpp/util/util.hh"
#include "abscab-cpp/abscab/abscab.hh"

namespace abscab {

/** test method for Straight Wire Segment methods */
TEST(TestAbscab, CheckStraightWireSegment) {

	int rowsRp = 0;
	int colsRp = 0;
	double **test_points_rp = loadColumnsFromFile("abscab-cpp/abscab/resources/testPointsRpStraightWireSegment.dat", &rowsRp, &colsRp);
	if (rowsRp < 1) {
		printf("error: need at least one row of test point coordinates for rho'\n");
		FAIL();
	}
	if (colsRp < 1) {
		printf("error: need at least one column of test point coordinates for rho'\n");
		FAIL();
	}

	int rowsZp = 0;
	int colsZp = 0;
	double **test_points_zp = loadColumnsFromFile("abscab-cpp/abscab/resources/testPointsZpStraightWireSegment.dat", &rowsZp, &colsZp);
	if (rowsZp < 1) {
		printf("error: need at least one row of test point coordinates for z'\n");
		FAIL();
	}
	if (colsZp < 1) {
		printf("error: need at least one column of test point coordinates for z'\n");
		FAIL();
	}

	if (rowsRp != rowsZp) {
		printf("error: number of rows of test point coordinates has to agree between rho' (%d) and z' (%d)\n", rowsRp, rowsZp);
		FAIL();
	}
	if (colsRp != colsZp) {
		printf("error: number of columns of test point coordinates has to agree between rho' (%d) and z' (%d)\n", colsRp, colsZp);
		FAIL();
	}

	int numCases = rowsRp;

	int rowsAZRef = 0;
	int colsAZRef = 0;
	double **aZRef = loadColumnsFromFile("abscab-cpp/abscab/resources/StraightWireSegment_A_z_ref.dat", &rowsAZRef, &colsAZRef);
	if (rowsAZRef < 1) {
		printf("error: need at least one row of reference values for A_z\n");
		FAIL();
	}
	if (colsAZRef < 1) {
		printf("error: need at least one column of reference values for A_z\n");
		FAIL();
	}

	if (numCases != rowsAZRef) {
		printf("error: number of reference values for A_z (%d) has to match number of test cases (%d)\n", rowsAZRef, rowsRp);
		FAIL();
	}

	int rowsBPhiRef = 0;
	int colsBPhiRef = 0;
	double **bPhiRef = loadColumnsFromFile("abscab-cpp/abscab/resources/StraightWireSegment_B_phi_ref.dat", &rowsBPhiRef, &colsBPhiRef);
	if (rowsBPhiRef < 1) {
		printf("error: need at least one row of reference values for B_phi\n");
		FAIL();
	}
	if (colsBPhiRef < 1) {
		printf("error: need at least one column of reference values for B_phi\n");
		FAIL();
	}

	if (numCases != rowsBPhiRef) {
		printf("error: number of reference values for B_phi (%d) has to match number of test cases (%d)\n", rowsBPhiRef, rowsRp);
		FAIL();
	}

	double toleranceAZ   = 1.0e-15;
	double toleranceBPhi = 1.0e-15;

	for (int i = 0; i < numCases; ++i) {
		int status = 0;

		double rp = test_points_rp[0][i];
		double zp = test_points_zp[0][i];

		// compute values using C implementation to test
		double aZ   = straightWireSegment_A_z(rp, zp);
		double bPhi = straightWireSegment_B_phi(rp, zp);

		int aZStatus = assertRelAbsEquals(aZRef[0][i], aZ, toleranceAZ);
		EXPECT_EQ(aZStatus, 0);
		if (aZStatus) {
			printf("error: mismatch at Straight Wire Segment A_z test case %d\n", i);
			printf("  rho' = %.17e\n", rp);
			printf("    z' = %.17e\n", zp);
			printf("  ref A_z = %+.17e\n", aZRef[0][i]);
			printf("  act A_z = %+.17e\n", aZ);
		}
		status |= aZStatus;

		int bPhiStatus = assertRelAbsEquals(bPhiRef[0][i], bPhi, toleranceBPhi);
		EXPECT_EQ(aZStatus, 0);
		if (bPhiStatus) {
			printf("error: mismatch at Straight Wire Segment B_phi test case %d\n", i);
			printf("  rho' = %.17e\n", rp);
			printf("    z' = %.17e\n", zp);
			printf("  ref B_phi = %+.17e\n", bPhiRef[0][i]);
			printf("  act B_phi = %+.17e\n", bPhi);
		}
		status |= bPhiStatus;

		if (status) {
			break;
		}
	}

	// free data loaded from text files
	free(test_points_rp[0]); free(test_points_rp);
	free(test_points_zp[0]); free(test_points_zp);
	free(aZRef[0]);          free(aZRef);
	free(bPhiRef[0]);        free(bPhiRef);
}

/** test method for Circular Wire Loop methods */
TEST(TestAbscab, CheckCircularWireLoop) {

	int rowsRp = 0;
	int colsRp = 0;
	double **test_points_rp = loadColumnsFromFile("abscab-cpp/abscab/resources/testPointsRpCircularWireLoop.dat", &rowsRp, &colsRp);
	if (rowsRp < 1) {
		printf("error: need at least one row of test point coordinates for rho'\n");
		FAIL();
	}
	if (colsRp < 1) {
		printf("error: need at least one column of test point coordinates for rho'\n");
		FAIL();
	}

	int rowsZp = 0;
	int colsZp = 0;
	double **test_points_zp = loadColumnsFromFile("abscab-cpp/abscab/resources/testPointsZpCircularWireLoop.dat", &rowsZp, &colsZp);
	if (rowsZp < 1) {
		printf("error: need at least one row of test point coordinates for z'\n");
		FAIL();
	}
	if (colsZp < 1) {
		printf("error: need at least one column of test point coordinates for z'\n");
		FAIL();
	}

	if (rowsRp != rowsZp) {
		printf("error: number of rows of test point coordinates has to agree between rho' (%d) and z' (%d)\n", rowsRp, rowsZp);
		FAIL();
	}
	if (colsRp != colsZp) {
		printf("error: number of columns of test point coordinates has to agree between rho' (%d) and z' (%d)\n", colsRp, colsZp);
		FAIL();
	}

	int numCases = rowsRp;

	int rowsAPhiRef = 0;
	int colsAPhiRef = 0;
	double **aPhiRef = loadColumnsFromFile("abscab-cpp/abscab/resources/CircularWireLoop_A_phi_ref.dat", &rowsAPhiRef, &colsAPhiRef);
	if (rowsAPhiRef < 1) {
		printf("error: need at least one row of reference values for A_phi\n");
		FAIL();
	}
	if (colsAPhiRef < 1) {
		printf("error: need at least one column of reference values for A_phi\n");
		FAIL();
	}

	if (numCases != rowsAPhiRef) {
		printf("error: number of reference values for A_phi (%d) has to match number of test cases (%d)\n", rowsAPhiRef, rowsRp);
		FAIL();
	}

	int rowsBRhoRef = 0;
	int colsBRhoRef = 0;
	double **bRhoRef = loadColumnsFromFile("abscab-cpp/abscab/resources/CircularWireLoop_B_rho_ref.dat", &rowsBRhoRef, &colsBRhoRef);
	if (rowsBRhoRef < 1) {
		printf("error: need at least one row of reference values for B_rho\n");
		FAIL();
	}
	if (colsBRhoRef < 1) {
		printf("error: need at least one column of reference values for B_rho\n");
		FAIL();
	}

	if (numCases != rowsBRhoRef) {
		printf("error: number of reference values for B_rho (%d) has to match number of test cases (%d)\n", rowsBRhoRef, rowsRp);
		FAIL();
	}

	int rowsBZRef = 0;
	int colsBZRef = 0;
	double **bZRef = loadColumnsFromFile("abscab-cpp/abscab/resources/CircularWireLoop_B_z_ref.dat", &rowsBZRef, &colsBZRef);
	if (rowsBZRef < 1) {
		printf("error: need at least one row of reference values for B_z\n");
		FAIL();
	}
	if (colsBZRef < 1) {
		printf("error: need at least one column of reference values for B_z\n");
		FAIL();
	}

	if (numCases != rowsBZRef) {
		printf("error: number of reference values for B_z (%d) has to match number of test cases (%d)\n", rowsBZRef, rowsRp);
		FAIL();
	}

	double toleranceAPhi = 1.0e-15;
	double toleranceBRho = 1.0e-13;
	double toleranceBZ   = 1.0e-14;

	for (int i = 0; i < numCases; ++i) {
		int status = 0;

		double rp = test_points_rp[0][i];
		double zp = test_points_zp[0][i];

		// compute values using C implementation to test
		double aPhi = circularWireLoop_A_phi(rp, zp);
		double bRho = circularWireLoop_B_rho(rp, zp);
		double bZ   = circularWireLoop_B_z(rp, zp);

		int aPhiStatus = assertRelAbsEquals(aPhiRef[0][i], aPhi, toleranceAPhi);
		EXPECT_EQ(aPhiStatus, 0);
		if (aPhiStatus) {
			printf("error: mismatch at Circular Wire Loop A_phi test case %d\n", i);
			printf("  rho' = %.17e\n", rp);
			printf("    z' = %.17e\n", zp);
			printf("  ref A_phi = %+.17e\n", aPhiRef[0][i]);
			printf("  act A_phi = %+.17e\n", aPhi);
		}
		status |= aPhiStatus;

		int bRhoStatus = assertRelAbsEquals(bRhoRef[0][i], bRho, toleranceBRho);
		EXPECT_EQ(bRhoStatus, 0);
		if (bRhoStatus) {
			printf("error: mismatch at Circular Wire Loop B_rho test case %d\n", i);
			printf("  rho' = %.17e\n", rp);
			printf("    z' = %.17e\n", zp);
			printf("  ref B_rho = %+.17e\n", bRhoRef[0][i]);
			printf("  act B_rho = %+.17e\n", bRho);
		}
		status |= bRhoStatus;

		int bZStatus = assertRelAbsEquals(bZRef[0][i], bZ, toleranceBZ);
		EXPECT_EQ(bZStatus, 0);
		if (bZStatus) {
			printf("error: mismatch at Circular Wire Loop B_z test case %d\n", i);
			printf("  rho' = %.17e\n", rp);
			printf("    z' = %.17e\n", zp);
			printf("  ref B_z = %+.17e\n", bZRef[0][i]);
			printf("  act B_z = %+.17e\n", bZ);
		}
		status |= bZStatus;
	}

	// free data loaded from text files
	free(test_points_rp[0]); free(test_points_rp);
	free(test_points_zp[0]); free(test_points_zp);
	free(aPhiRef[0]);        free(aPhiRef);
	free(bRhoRef[0]);        free(bRhoRef);
	free(bZRef[0]);          free(bZRef);
}

TEST(TestAbscab, CheckMagneticFieldInfiniteLineFilament) {
	double tolerance = 1.0e-15;

	// Demtroeder 2, Sec. 3.2.2 ("Magnetic field of a straight wire")
	// B(r) = mu_0 * I / (2 pi r)
	// Test this here with:
	// I = 123.0 A
	// r = 0.132 m
	// => B = 0.186 mT
	double current = 123.0;
	double r = 0.132;
	double bPhiRef = MU_0 * current / (2.0 * M_PI * r);
//	printf("ref bPhi = %.5e\n", bPhiRef);

	double vertices[] = {
			0.0, 0.0, -1.0e6,
			0.0, 0.0,  1.0e6
	};

	double evalPos[] = {
			r, 0.0, 0.0
	};

	// y component is B_phi
	double magneticField[3];
	magneticFieldPolygonFilament(2, vertices, current, 1, evalPos, magneticField);
	double bPhi = magneticField[1];
//	printf("act bPhi = %.5e\n", bPhi);

//	double relAbsErr = fabs(bPhi - bPhiRef) / (1.0 + fabs(bPhiRef));
//	printf("raErr = %.5e\n", relAbsErr);

	EXPECT_EQ(assertRelAbsEquals(bPhiRef, bPhi, tolerance), 0);
}

TEST(TestAbscab, CheckBPhiInfiniteLineFilament) {
	double tolerance = 1.0e-15;

	// Demtroeder 2, Sec. 3.2.2 ("Magnetic field of a straight wire")
	// B(r) = mu_0 * I / (2 pi r)
	// Test this here with:
	// I = 123.0 A
	// r = 0.132 m
	// => B = 0.186 mT
	double current = 123.0;
	double r = 0.132;
	double bPhiRef = MU_0 * current / (2.0 * M_PI * r);
//	printf("ref bPhi = %.5e\n", bPhiRef);

	// half the length of the wire segment
	double halfL = 1e6;
	double L = 2*halfL;
	double rhoP = r / L;
	double zP = halfL / L;
	double bPhi = MU_0 * current / (4.0 * M_PI * L) * straightWireSegment_B_phi(rhoP, zP);
//	printf("act bPhi = %.5e\n", bPhi);

//	double relAbsErr = fabs(bPhi - bPhiRef) / (1.0 + fabs(bPhiRef));
//	printf("raErr = %.5e\n", relAbsErr);

	EXPECT_EQ(assertRelAbsEquals(bPhiRef, bPhi, tolerance), 0);
}

TEST(TestAbscab, CheckMagneticFieldInsideLongCoil) {
	double tolerance = 1.0e-4;

	// Demtroeder 2, Sec. 3.2.3 ("Magnetic field of a long coil")
	// B_z = mu_0 * n * I
	// where n is the winding density: n = N / L
	// of a coil of N windings over a length L
	// Example (which is tested here):
	// n = 1e3 m^{-1}
	// I = 10 A
	// => B = 0.0126T
	double bZRef = 0.0126;

	int N = 50000; // windings
	double L = 50.0; // total length of coil in m
	double n = N/L;

	double current = 10.0; // A
	double radius = 1.0; // m

	double bZ = 0.0;
	for (int i = 0; i < N; ++i) {

		// axial position of coil
		double z0 = -L/2.0 + (i + 0.5) / n;

		// compute magnetic field
		double prefac = MU_0 * current / (M_PI * radius);
		double bZContrib = prefac * circularWireLoop_B_z(0.0, z0);

//		printf("coil %d at z0 = % .3e => contrib = %.3e\n", i, z0, bZContrib);

		bZ += bZContrib;
	}
//	printf("B_z = %.5e\n", bZ);

//	double relAbsErr = fabs(bZ - bZRef) / (1.0 + fabs(bZRef));
//	printf("raErr = %.5e\n", relAbsErr);

	EXPECT_EQ(assertRelAbsEquals(bZRef, bZ, tolerance), 0);
}

// check that previous field contributions are preserved
// for circular wire loop methods
TEST(TestAbscab, CheckPreseveranceCircularWireLoop) {

	double tolerance = 1.0e-15;

	double center[3] = { 1.0, 2.0, 3.0 };
	double normal[3] = { 4.0, 5.0, 6.0 };
	double radius = 7.89;

	double current_1 = 12.3;
	double current_2 = 45.6;

	double evalPos[3] = { 7.0, 8.0, 9.0 };

	// reference: add two contributions
	double vectorPotential_ref[3] = { 0.0, 0.0, 0.0 };
	vectorPotentialCircularFilament(center, normal, radius, current_1 + current_2, 1, evalPos, vectorPotential_ref);

	// to be tested: add 
	double vectorPotential[3] = { 0.0, 0.0, 0.0 };
	vectorPotentialCircularFilament(center, normal, radius, current_1, 1, evalPos, vectorPotential);
	vectorPotentialCircularFilament(center, normal, radius, current_2, 1, evalPos, vectorPotential);

	EXPECT_EQ(assertRelAbsEquals(vectorPotential_ref[0], vectorPotential[0], tolerance), 0);
	EXPECT_EQ(assertRelAbsEquals(vectorPotential_ref[1], vectorPotential[1], tolerance), 0);
	EXPECT_EQ(assertRelAbsEquals(vectorPotential_ref[2], vectorPotential[2], tolerance), 0);

	// and again for the magnetic field ...
	double magneticField_ref[3] = { 0.0, 0.0, 0.0 };
	magneticFieldCircularFilament(center, normal, radius, current_1 + current_2, 1, evalPos, magneticField_ref);

	double magneticField[3] = { 0.0, 0.0, 0.0 };
	magneticFieldCircularFilament(center, normal, radius, current_1, 1, evalPos, magneticField);
	magneticFieldCircularFilament(center, normal, radius, current_2, 1, evalPos, magneticField);

	EXPECT_EQ(assertRelAbsEquals(magneticField_ref[0], magneticField[0], tolerance), 0);
	EXPECT_EQ(assertRelAbsEquals(magneticField_ref[1], magneticField[1], tolerance), 0);
	EXPECT_EQ(assertRelAbsEquals(magneticField_ref[2], magneticField[2], tolerance), 0);
}

TEST(TestAbscab, CheckPreseverancePolygonFilamentArrays) {

	double tolerance = 1.0e-15;

	int numVertices = 4;
	double vertices[numVertices * 3] = {
		0.0, 0.0, 0.0,
		1.0, 2.0, 3.0,
		4.0, 5.0, 6.0,
		7.0, 8.0, 9.0
	};

	double current_1 = 1.23;
	double current_2 = 4.56;

	for (int numProcessors = 1; numProcessors < 3; ++numProcessors) {
		for (int compensatedSummationCase = 0; compensatedSummationCase < 2; ++compensatedSummationCase) {
			bool useCompensatedSummation = (compensatedSummationCase == 0);
			for (int numEvalPos = 1; numEvalPos < 6; numEvalPos += 4) {

				int numBytes = numEvalPos * 3 * sizeof(double);

				double *evalPos = (double *) malloc(numBytes);
				if (evalPos == NULL) {
					FAIL() << "could not allocate evalPos";
				}
				for (int idxEvalPos = 0; idxEvalPos < numEvalPos; ++idxEvalPos) {
					evalPos[idxEvalPos * 3 + 0] = -(idxEvalPos + 1);
					evalPos[idxEvalPos * 3 + 1] = -(idxEvalPos + 1);
					evalPos[idxEvalPos * 3 + 2] = -(idxEvalPos + 1);
				}

				double *vectorPotential_ref = (double *) malloc(numBytes);
				if (vectorPotential_ref == NULL) {
					FAIL() << "could not allocate vectorPotential_ref";
				}
				memset(vectorPotential_ref, 0, numBytes);

				vectorPotentialPolygonFilament(numVertices, vertices, current_1 + current_2, numEvalPos, evalPos, vectorPotential_ref, numProcessors, useCompensatedSummation);

				double *vectorPotential = (double *) malloc(numBytes);
				if (vectorPotential == NULL) {
					FAIL() << "could not allocate vectorPotential";
				}
				memset(vectorPotential, 0, numBytes);

				vectorPotentialPolygonFilament(numVertices, vertices, current_1, numEvalPos, evalPos, vectorPotential, numProcessors, useCompensatedSummation);
				vectorPotentialPolygonFilament(numVertices, vertices, current_2, numEvalPos, evalPos, vectorPotential, numProcessors, useCompensatedSummation);

				for (int idxEvalPos = 0; idxEvalPos < numEvalPos; ++idxEvalPos) {
					EXPECT_EQ(assertRelAbsEquals(vectorPotential_ref[idxEvalPos * 3 + 0], vectorPotential[idxEvalPos * 3 + 0], tolerance), 0);
					EXPECT_EQ(assertRelAbsEquals(vectorPotential_ref[idxEvalPos * 3 + 1], vectorPotential[idxEvalPos * 3 + 1], tolerance), 0);
					EXPECT_EQ(assertRelAbsEquals(vectorPotential_ref[idxEvalPos * 3 + 2], vectorPotential[idxEvalPos * 3 + 2], tolerance), 0);
				}

				free(vectorPotential);
				free(vectorPotential_ref);
				free(evalPos);
			}
		}
	}
}

TEST(TestAbscab, CheckPreseverancePolygonFilamentVertexSupplier) {

	double tolerance = 1.0e-15;

	constexpr int numVertices = 4;
	void (*vertexSupplier)(int, double *) = [](int i, double *point) {
		switch (i) {
		case 0:
			point[0] = 0.0; point[1] = 0.0; point[2] = 0.0;
			break;
		case 1:
			point[0] = 1.0; point[1] = 2.0; point[2] = 3.0;
			break;
		case 2:
			point[0] = 4.0; point[1] = 5.0; point[2] = 6.0;
			break;
		case 3:
			point[0] = 7.0; point[1] = 8.0; point[2] = 9.0;
			break;
		default:
			std::stringstream msg;
			msg << "vertexSupplier only provides " << numVertices << " vertices";
			FAIL() << msg.str();
			break;
		}
	};

	double current_1 = 1.23;
	double current_2 = 4.56;

	for (int numProcessors = 1; numProcessors < 3; ++numProcessors) {
		for (int compensatedSummationCase = 0; compensatedSummationCase < 2; ++compensatedSummationCase) {
			bool useCompensatedSummation = (compensatedSummationCase == 0);
			for (int numEvalPos = 1; numEvalPos < 6; numEvalPos += 4) {

				int numBytes = numEvalPos * 3 * sizeof(double);

				double *evalPos = (double *) malloc(numBytes);
				if (evalPos == NULL) {
					FAIL() << "could not allocate evalPos";
				}
				for (int idxEvalPos = 0; idxEvalPos < numEvalPos; ++idxEvalPos) {
					evalPos[idxEvalPos * 3 + 0] = -(idxEvalPos + 1);
					evalPos[idxEvalPos * 3 + 1] = -(idxEvalPos + 1);
					evalPos[idxEvalPos * 3 + 2] = -(idxEvalPos + 1);
				}

				double *vectorPotential_ref = (double *) malloc(numBytes);
				if (vectorPotential_ref == NULL) {
					FAIL() << "could not allocate vectorPotential_ref";
				}
				memset(vectorPotential_ref, 0, numBytes);

				vectorPotentialPolygonFilament(numVertices, vertexSupplier, current_1 + current_2, numEvalPos, evalPos, vectorPotential_ref, numProcessors, useCompensatedSummation);

				double *vectorPotential = (double *) malloc(numBytes);
				if (vectorPotential == NULL) {
					FAIL() << "could not allocate vectorPotential";
				}
				memset(vectorPotential, 0, numBytes);

				vectorPotentialPolygonFilament(numVertices, vertexSupplier, current_1, numEvalPos, evalPos, vectorPotential, numProcessors, useCompensatedSummation);
				vectorPotentialPolygonFilament(numVertices, vertexSupplier, current_2, numEvalPos, evalPos, vectorPotential, numProcessors, useCompensatedSummation);

				for (int idxEvalPos = 0; idxEvalPos < numEvalPos; ++idxEvalPos) {
					EXPECT_EQ(assertRelAbsEquals(vectorPotential_ref[idxEvalPos * 3 + 0], vectorPotential[idxEvalPos * 3 + 0], tolerance), 0);
					EXPECT_EQ(assertRelAbsEquals(vectorPotential_ref[idxEvalPos * 3 + 1], vectorPotential[idxEvalPos * 3 + 1], tolerance), 0);
					EXPECT_EQ(assertRelAbsEquals(vectorPotential_ref[idxEvalPos * 3 + 2], vectorPotential[idxEvalPos * 3 + 2], tolerance), 0);
				}

				free(vectorPotential);
				free(vectorPotential_ref);
				free(evalPos);
			}
		}
	}
}



} // namespace abscab

