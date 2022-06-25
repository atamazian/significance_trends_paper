/** \file main.c
 * NAME:        mhist
 *
 * DESCRIPTION: reads data from file, calculates histogram, probability
 *              distribution or exceedance probability using n bins from a to
 *              b and prints result to stdout.
 *
 * INPUTS:      infile     - text file with one-column double data x
 *              double a   - minimum value of x
 *              double b   - maximum value of x
 *              size_t n   - number of bins
 *
 * OUTPUTS:     stdout - two-column double data of the following form:
 *              [no flag]: [bin upper]  [frequency]
 *              -p mode: [bin center] [probability]
 *              -c mode: [bin center] [cumulative probability]
 *              -e mode: [bin center] [exceedance probability]
 *
 * USAGE:       mhist [infile] [options]
 *               -p estimate probability density function
 *               -c estimate cumulative distribution function
 *               -e estimate exceedance probability
 *               -a start point
 *               -b end point
 *               -n number of bins
 *               -w bin width
 *
 * AUTHOR:      Araik Tamazian
 *
 ******************************************************************************/

#include <getopt.h>
#include <gsl/gsl_histogram.h>
#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/** Count number of lines in a file
*
* @param fname is a source file name
*
* @return number of lines in source file
*/
int get_file_lcount(const char *fname) {
  int l = 0;

  /* Open the file */
  FILE *fin = fopen(fname, "r");

  /* Check if this file can be opened */
  if (fin == NULL) {
    printf("Error: unable to open file %s\n", fname);
    exit(1);
  }

  /* Count number of lines in the file */
  for (char c = 0; (c = fgetc(fin)) && c != EOF; l += c == '\n')
    ;

  /* Close the file */
  fclose(fin);

  return l;
}

/**
* Read data from file
*
* @param fname is a source file name
* @param l is a length of source file
*
* @return data [double array] from source file
*
*/
double *read_file(const char *fname, size_t l) {

  FILE *fin = fopen(fname, "r");

  /* Check if this file can be opened */
  if (fin == NULL) {
    printf("error: unable to open file %s\n", fname);
    exit(1);
  }

  /* Read data from the file to array */
  double *x = malloc(l * sizeof(double));

  for (size_t i = 0; i < l; i++)
    fscanf(fin, "%lf", &x[i]);

  fclose(fin);

  return x;
}

/**
* Calculates histogram
*
* @param x is a source data
* @param l is a length of source data
* @param a is a left point of histogram range
* @param b is a riht point of histogram range
* @param n is a number of bins
*
* @return Nothing
*
*/
void calculate_hist(const double *x, int l, double a, double b, size_t n) {

  /* Create histogram and set its bins */
  gsl_histogram *h = gsl_histogram_alloc(n);
  gsl_histogram_set_ranges_uniform(h, a, b);

  for (size_t i = 0; i < l; i++) {
    gsl_histogram_increment(h, x[i]);
  }

  double lw[1], up[1];

  for (size_t i = 0; i < n; i++) {
    gsl_histogram_get_range(h, i, lw, up);
    printf("%lf %lf\n", up[0], gsl_histogram_get(h, i));
  }

  /* Remove histogram data from memory */
  gsl_histogram_free(h);
}

/**
* Calculates probability density estimate
*
* @param x is a source data
* @param l is a length of source data
* @param a is a left point of histogram range
* @param b is a riht point of histogram range
* @param n is a number of bins
*
* @return Nothing
*
*/
void calculate_pdf(const double *x, int l, double a, double b, size_t n) {

  /* Create histogram and set its bins */
  gsl_histogram *h = gsl_histogram_alloc(n);
  gsl_histogram_set_ranges_uniform(h, a, b);

  for (size_t i = 0; i < l; i++) {
    gsl_histogram_increment(h, x[i]);
  }

  double lw[1], up[1];

  /* Calculate bin width */
  double d = (b - a) / ((double)n);

  /* Normalize histogram */
  gsl_histogram_scale(h, 1.0 / ((double)l * d));

  for (size_t i = 0; i < n; i++) {
    gsl_histogram_get_range(h, i, lw, up);
    printf("%lf %lf\n", (lw[0] + up[0]) / 2.0, gsl_histogram_get(h, i));
  }

  gsl_histogram_free(h);
}

/**
* Calculates cumulative distribution function
*
* @param x is a source data
* @param l is a length of source data
* @param a is a left point of histogram range
* @param b is a riht point of histogram range
* @param n is a number of bins
*
* @return Nothing
*
*/
void calculate_cdf(const double *x, int l, double a, double b, size_t n) {

  /* Create histogram and set its bins */
  gsl_histogram *h = gsl_histogram_alloc(n);
  gsl_histogram_set_ranges_uniform(h, a, b);

  for (size_t i = 0; i < l; i++) {
    gsl_histogram_increment(h, x[i]);
  }

  gsl_histogram_scale(h, 1.0 / (double)l);

  double lw[1] = {0.0};
  double up[1] = {0.0};
  double cp = 0;
  for (size_t i = 0; i < n; i++) {
    printf("%lf %lf\n", (lw[0] + up[0]) / 2.0, cp);
    gsl_histogram_get_range(h, i, lw, up);
    cp += gsl_histogram_get(h, i);
  }

  gsl_histogram_free(h);
}

/**
* Calculates exceedance probability
*
* @param x is a source data
* @param l is a length of source data
* @param a is a left point of histogram range
* @param b is a riht point of histogram range
* @param n is a number of bins
*
* @return Nothing
*
*/
void calculate_edf(const double *x, int l, double a, double b, size_t n) {

  /* Create histogram and set its bins */
  gsl_histogram *h = gsl_histogram_alloc(n);
  gsl_histogram_set_ranges_uniform(h, a, b);

  for (size_t i = 0; i < l; i++) {
    gsl_histogram_increment(h, x[i]);
  }

  gsl_histogram_scale(h, 1.0 / (double)l);

  double lw[1] = {0.0};
  double up[1] = {0.0};
  double cp = 0;
  for (size_t i = 0; i < n; i++) {
    printf("%lf %lf\n", (lw[0] + up[0]) / 2.0, 1.0 - cp);
    gsl_histogram_get_range(h, i, lw, up);
    cp += gsl_histogram_get(h, i);
  }

  gsl_histogram_free(h);
}

/**
* Calculates min value
*
* @param x is a source data
* @param l is a length of source data
*
* @return Min value of x
*
*/
double get_min(const double *x, size_t l) {
  /* Find minimum x value */
  double x_min = x[0];
  for (size_t i = 0; i < l; i++) {
    if (x[i] < x_min) {
      x_min = x[i];
    }
  }
  return x_min;
}

/**
* Calculates max value
*
* @param x is a source data
* @param l is a length of source data
*
* @return Max value of x
*
*/
double get_max(const double *x, size_t l) {
  /* Find maximum x value */
  double x_max = x[0];
  for (size_t i = 0; i < l; i++) {
    if (x[i] > x_max) {
      x_max = x[i];
    }
  }
  return x_max;
}

/**
*
* Displays usage hint
*
* @return Nothing
*/
void print_usage() {
  printf("Usage: mhist [infile] [options]\n"
         " -p estimate probability density function\n"
         " -c estimate cumulative distribution function\n"
         " -e estimate exceedance probability\n"
         " -a start point\n"
         " -b end point\n"
         " -n number of bins. Must be greater than zero.\n"
         " -w bin width. Must be positive.\n");
}


int main(int argc, char *argv[]) {

  /* Disable error messages */
  opterr = 0;

  /* Check if file name is specified */
  if (argv[1] == NULL) {
    print_usage();
    exit(1);
  }

  /* Get input file name */
  char *infile;
  asprintf(&infile, "%s", argv[1]);

  /* Read program meters */
  int option = 0;
  char mode = 0;
  double a = NAN, b = NAN, d = NAN;
  size_t n = 0;
  while ((option = getopt(argc, argv, "spcea:b:n:w:g")) != -1) {
    switch (option) {
    case 'p':
    case 'c':
    case 'e':
      mode = option;
      break;
    case 'a':
      a = atof(optarg);
      break;
    case 'b':
      b = atof(optarg);
      break;
    case 'n':
      n = atoi(optarg);
      break;
    case 'w':
      d = atof(optarg);
      break;
    default:
      print_usage();
      exit(EXIT_FAILURE);
    }
  }

  /* Terminate program if not all necessary meters were set*/
  if (((n == 0) && isnan(d)) || ((n != 0) && !isnan(d))) {
    print_usage();
    exit(EXIT_FAILURE);
  }

  /* Check if bin width has a positive, non-zero value */
  if ((d <= 0) && (n == 0)) {
    print_usage();
    exit(EXIT_FAILURE);
  }

  /* Count lines in input file */
  size_t l = get_file_lcount(infile);

  /* Allocate memory for data */
  double *x = read_file(infile, l);

  /* Calculate a and/or b if not provided */
  if (isnan(a)) {
    a = get_min(x, l);
  }

  if (isnan(b)) {
    b = get_max(x, l);
  }

  /* Calculate number of bins in case of bin width provided */
  if (n == 0)
    n = (int)trunc((b - a) / d);
  
  switch(mode) {
    case 'p':
	  /* MODE -p. Calculate probability distribution and print it */  
      calculate_pdf(x, l, a, b, n);	
	  break;
    case 'c':
	  /* MODE -c. Calculate cumulative probability and print it */  
      calculate_cdf(x, l, a, b, n); 
	  break;
    case 'e':
	  /* MODE -e. Calculate exceedance probability and print it */ 
      calculate_edf(x, l, a, b, n);
	  break;
	default:
      /* DEFAULT MODE. Print histogram */
      calculate_hist(x, l, a, b, n);
  }

  free(infile);
  free(x);
}
