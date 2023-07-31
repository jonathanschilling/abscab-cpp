#ifndef UTIL_H
#define UTIL_H

#include <ctype.h>
#include <stdio.h>
#include <stdbool.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

namespace abscab {

/**
 * Check if two values are approximately equal within a prescribed tolerance.
 * For values much smaller than 1, this is similar to a comparison of the
 * absolute values. For values much greater than 1, this is similar to a
 * comparison of the relative values.
 *
 * This method is described in Gill, Murray & Wright, "Practical Optimization" (1984).
 *
 * @param expected  expected result
 * @param actual    actual result
 * @param tolerance relative or absolute tolerance on the mismatch between the
 *                  expected and the actual values
 * @return 0 if the values match within the prescribed tolerance; 1 otherwise
 */
int assertRelAbsEquals(double expected, double actual, double tolerance);

/**
 * Error metric for computing the "number of matching digits" between two numbers.
 */
double errorMetric(double ref, double act);

// https://stackoverflow.com/a/122721
// Note: This function returns a pointer to a substring of the original string.
// If the given string was allocated dynamically, the caller must not overwrite
// that pointer with the returned value, since the original pointer must be
// deallocated using the same allocator with which it was allocated.  The return
// value must NOT be deallocated using free() etc.
char* trim_whitespace(char *str);

/**
 * Load columns of data from a text file.
 * The number of rows found will be written into numRows.
 * The number of columns found will be written into numColumns.
 * The return value is an array of arrays (pointer to an array of pointers),
 * where one entry corresponds to one column of the data in the file.
 * Rows in the file starting with "#" are not counted and not parsed.
 */
double** loadColumnsFromFile(const char *filename, int *numRows, int *numColumns);

void dumpToFile(int numCols, int numRows, double *data, char *filename);

} // namespace abscab

#endif // UTIL_H
