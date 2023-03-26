/*  SPECTRUM.C

%   Copyright 1997 E. I. du Pont de Nemours and Company              	 %
%                                                                      	 %
%   Permission to use, copy, modify, distribute, and sell this software and %
%   its documentation for any purpose is hereby granted without fee,        %
%   provided that the above Copyright notice appear in all copies and that  %
%   both that Copyright notice and this permission notice appear in         %
%   supporting documentation, and that the name of E. I. du Pont de Nemours %
%   and Company not be used in advertising or publicity pertaining to       %
%   distribution of the software without specific, written prior            %
%   permission.  E. I. du Pont de Nemours and Company makes no representations %
%   about the suitability of this software for any purpose.  It is provided %
%   "as is" without express or implied warranty.                            %
%  									    %
%   E. I. du Pont de Nemours and Company disclaims all warranties with regard %
%   to this software, including all implied warranties of merchantability   %
%   and fitness, in no event shall E. I. du Pont de Nemours and Company be  %
%   liable for any special, indirect or consequential damages or any        %
%   damages whatsoever resulting from loss of use, data or profits, whether %
%   in an action of contract, negligence or other tortious action, arising  %
%   out of or in connection with the use or performance of this software.   % 

    spectrum :

            This program takes a set of transtion frequencies and
        intensities and produces a data set suitable for plotting.  By
        choosing a reasonable valueof HWHM you can create a 
        realistic-looking specturm.

            The input file should have one transition frequency and 
        intensity per line.

         -i   <string>  Input file name of transition frequencies and intensities
         -o   <string>  Output file of a plottable spectrum

        [-l]            Lineshape
                             0 (Lorentzian),
                             1 (Gaussian),
                          or 2 (Delta function - default)
        [-h]  <number>  HWHM (Half-Width at Half-Maximum).  Ignored if the
                          lineshape is the delta function
        [-sc] <number>  Frequency scale factor.  Multiply all input frequencies
                          by this value.  (does not scale HWHM)
        [-sf] <number>  Starting frequency of the output file
        [-ef] <number>  Ending frequency
        [-fs] <number>  Frequency step size.  Ignored if the lineshape
	                  is the delta function
        [-d]            Sort output by decreasing frequency
                          (default is by increasing frequency)
        [-n]            Normalize intensities to 100
        [-nf] <number>  Normalization factor.  Divide all output intensities
                          by this value.
        [-wnf]          Write the normalization factor to a file.  (This allows
                          a factor to be shared among spectra for consistent
                          normalization.)  The name of the file is that of the
                          output file with \".nf\" appended.
        [-w]            After generating the spectrum from frequencies in wavenumbers,
                             write the ouput file with wavelengths in nanometers.

    Quantum mechanics programs such as Gaussian 92, DGAUSS, MNDO91, and MOPAC
    can calculate vibrational and other transition frequencies.  However, these
    are given as (frequency, intensity) pairs.  They have no width and cannot
    even be plotted directly as a delta function.  This program takes those data 
    and creates a more realistic spectrum which is suitable for plotting.
    
    The input file should have one frequency and one intensity per line, 
    separated by whitespace (spaces or tabs).  The output file has the same 
    format.

    This program uses only functions included in the ANSI C specification.

    Paul D. Soper
 
    DuPont
    Central Research & Development
    E328/123
    DUCOM 695-1757
    ESVAX::SOPERPD

Headers and function prototypes
*/


#include <ctype.h>
#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>
#include <time.h>

double    column_dmax(double** matrix, int n_rows, int col_index);
double    gauss(double rel_offset);   
double    lorentz(double rel_offset);   
double    next_higher(double d, double base);
double    next_lower(double d, double base);
double    order_of_magnitude(double d);

double**  calc_spectrum(double** in_spec, int n_lines, int lorentzian, int decreasing, double hwhm,double freq_step, double start_freq, double end_freq, int* n_freqs_ptr, int wavelengths, int exclude_transitions);
double**  delta_spectrum(double** in_spec, int n_lines, int decreasing, double start_freq, double end_freq, int* n_freqs_ptr);
double**  dmatrix_malloc(int nrows, int ncols);

int       count_frequencies(double** in_spec, int n_lines, double min_freq, double max_freq, double freq_step, int exclude_transitions);
int       dpcomp(const void* ptr1, const void* ptr2);
int       dpcomp_comp(const void* ptr1, const void* ptr2);
int       linecount(FILE* file);

void      calc_intensities(double** in_spec, int n_lines, double** out_spec, int n_freqs, int lorentzian, double hwhm);
void      calc_parameters(double** in_spec, int n_lines, int wavelengths, double hwhm, double* min_freq_ptr, double* max_freq_ptr, double* freq_step_ptr);
void      dmatrix_free(double** matrix, int nrows, int ncols);
void      dscale_column(double** matrix, int n_rows, int col_index, double scale);
void      limits(double low, double high, double* low_limit, double* high_limit);
void      normalize_spectrum(double** out_spec, int n_freqs, double n_value, int normalize, int write_norm_factor, char* outfile_name);
void      set_frequencies(double** in_spec, int n_lines, double** out_spec, int n_freqs, double min_freq, double max_freq, double freq_step, int exclude_transitions);

#define DZERO 1.0e-10


void print_usage(FILE* fp);



int main(int argc, char** argv) {

  char infile_name[101];
  char outfile_name[101];

  FILE* infile;
  FILE* outfile;

  double** in_spec;
  double** out_spec;

  int normalize, decreasing, lorentzian, n_lines, n_freqs, i;
  int write_norm_factor, wavelengths, lineshape, exclude_transitions;
  double start_freq, end_freq, freq_step, hwhm, n_value, freq_scale_factor;

  /* Set the default values of the command line parameters */

  infile_name[0] = '\0';
  outfile_name[0] = '\0';
  lineshape = 2;
  hwhm = 0.0;
  freq_scale_factor = 1.0;
  start_freq = -1.0;
  end_freq = -1.0;
  freq_step = 0.0;
  normalize = 0;
  n_value = 0.0;
  write_norm_factor = 0;
  decreasing = 0;
  wavelengths = 0;
  exclude_transitions = 0;

  /* Read and check the command line parameters.  Recall that
     atof() returns 0.0 if the argument contains an invalid
     character.  */

  if (argc == 1) {
    print_usage(stderr);
    exit(1);
  }

  for (i=1; i < argc; i++) {
    if (!strcmp(argv[i], "-i") && argc > i+1 && argv[i+1][0] != '-') {
      if (strlen(argv[i+1]) > 100) {
        fprintf(stderr, "The input file name is too long.\n");
        fprintf(stderr, "File names are limited to 100 characters\n");
        exit(1);
      }
      else strcpy(infile_name, argv[++i]);
    }
    else if (!strcmp(argv[i], "-o") && argc > i+1 && argv[i+1][0] != '-') {
      if (strlen(argv[i+1]) > 100) {
        fprintf(stderr, "The output file name is too long.\n");
        fprintf(stderr, "File names are limited to 100 characters\n");
        exit(1);
      }
      else strcpy(outfile_name, argv[++i]);
    }
    else if (!strcmp(argv[i], "-l") && argc > i+1) {
      lineshape = atoi(argv[++i]);
    }
    else if (!strcmp(argv[i], "-h") && argc > i+1) {
      hwhm = atof(argv[++i]);
    }
    else if (!strcmp(argv[i], "-sc") && argc > i+1) {
      freq_scale_factor = atof(argv[++i]);
    }
    else if (!strcmp(argv[i], "-sf") && argc > i+1) {
      start_freq = atof(argv[++i]);
    }
    else if (!strcmp(argv[i], "-ef") && argc > i+1) {
      end_freq = atof(argv[++i]);
    }
    else if (!strcmp(argv[i], "-fs") && argc > i+1) {
      freq_step = atof(argv[++i]);
      if (freq_step > DZERO) {
	exclude_transitions = 1;
      }
    }
    else if (!strcmp(argv[i], "-nf") && argc > i+1) {
      normalize = 1;
      n_value = atof(argv[++i]);
    }
    else if (!strcmp(argv[i], "-n")) {
      normalize = 1;
      n_value = 0.0;
    }
    else if (!strcmp(argv[i], "-wnf")) {
      write_norm_factor = 1;
    }
    else if (!strcmp(argv[i], "-d")) {
      decreasing = 1;
    }
    else if (!strcmp(argv[i], "-w")) {
      wavelengths = 1;
    }
    else {
      printf("  Unrecognized command line flag or parameter: %s\n", argv[i]);
      print_usage(stdout);
      exit(1);
    }
  }

  switch (lineshape) {
  case 0:
    lorentzian = 1;
    break;
  case 1:
    lorentzian = 0;
    break;
  case 2:
    lorentzian = 0;
    hwhm = 0.0;
    break;
  default:
    lorentzian = 0;
    hwhm = 0.0;
  }

  if (strlen(infile_name) == 0 || strlen(outfile_name) == 0) {
    fprintf(stderr, "\nPlease specify the input and output files.\n");  
    fprintf(stderr, "Both are required.\n");
    print_usage(stderr);
    exit(1);
  }
 
  /* Verify that both the input and output files can be opened. */

  infile = fopen(infile_name, "r");
  if (infile == NULL) {
    fprintf(stderr, "Error reading input file %s\n", infile_name);
    exit(1);
  }    

  outfile = fopen(outfile_name, "w");
  if (outfile == NULL) {
    fprintf(stderr, "Error opening output file %s\n", outfile_name);
    fclose(infile);
    exit(1);
  }

  /* Allocate memory for the input spectrum, read the input file,
     and then close it. */

  n_lines = linecount(infile);
  in_spec = dmatrix_malloc(n_lines, 2);
  if (in_spec == NULL) {
    fprintf(stderr, "Memory allocation error of in_spec in main()\n");
    exit(1);
  }

  rewind(infile);
  for (i=0; i < n_lines; i++) {
    fscanf(infile, "%lf %lf", &in_spec[i][0], &in_spec[i][1]);
    in_spec[i][0] *= freq_scale_factor;
  }
  fclose(infile);

  /* Call the function which calculates the output spectrum */

  out_spec = calc_spectrum(in_spec, n_lines, lorentzian, decreasing, hwhm, 
                           freq_step, start_freq, end_freq, &n_freqs,
			   wavelengths, exclude_transitions);

  if (out_spec == NULL) {
    fprintf(stderr, "Memory allocation error\n");
    fclose(outfile);
    dmatrix_free(in_spec, n_lines, 2);
    exit(1);
  }

  normalize_spectrum(out_spec, n_freqs, n_value, normalize, write_norm_factor, outfile_name);

  /* Write the output spectrum and close the output file. */

  rewind(outfile);
  if (wavelengths) {
    for (i=0; i < n_freqs; i++) {
      if (out_spec[i][0] > DZERO) {
	fprintf(outfile, "%20f\t%20f\n", 1.0e+07/out_spec[i][0], out_spec[i][1]);
      }
    }
  }
  else {
    for (i=0; i < n_freqs; i++) {
      fprintf(outfile, "%20f\t%20f\n", out_spec[i][0], out_spec[i][1]);
    }
  }
  fclose(outfile);

  /* Free allocated memory and exit */

  dmatrix_free(in_spec, n_lines, 2);
  dmatrix_free(out_spec, n_freqs, 2);

  exit(0);
}


/*
%==============================================================================%
%                                                                              %
%  void print_usage(FILE* fp);                                                 %
%                                                                              %
%  Print a description of the spectrum program and its input parameters        %
%  to (*fp).  In this program (*fp) is either stderr or stdout.                %
%                                                                              %
%  Side Effects:                                                               %
%                                                                              %
%    None                                                                      %
%                                                                              %
%  Requires:  <stdio.h>                                                        %
%                                                                              %
%==============================================================================%
*/

void print_usage(FILE* fp) {
    fprintf(fp, "spectrum\n");
    fprintf(fp, "      This program takes a set of transition frequencies and intensities\n");
    fprintf(fp, "  and produces a data set suitable for plotting.  By choosing a reasonable\n");
    fprintf(fp, "  value of HWHM you can create a realistic-looking spectrum.\n");
    fprintf(fp, "      The input file should have one transition frequency and intensity\n");
    fprintf(fp, "  per line.\n\n");
    fprintf(fp, "   -i   <string>  Input file of transition frequencies and intensities\n");
    fprintf(fp, "   -o   <string>  Output file of a plottable spectrum\n\n");
    fprintf(fp, "  [-l]            Lineshape\n");
    fprintf(fp, "                       0 (Lorentzian),\n");
    fprintf(fp, "                       1 (Gaussian),\n");
    fprintf(fp, "                    or 2 (Delta function - default)\n");
    fprintf(fp, "  [-h]  <number>  HWHM (Half-Width at Half-Maximum).  Ignored if the\n");
    fprintf(fp, "                    lineshape is the delta function\n");
    fprintf(fp, "  [-sc] <number>  Frequency scale factor.  Multiply all input frequencies\n");
    fprintf(fp, "                    by this value.  (does not scale HWHM)\n");
    fprintf(fp, "  [-sf] <number>  Starting frequency of the output file\n");
    fprintf(fp, "  [-ef] <number>  Ending frequency\n");
    fprintf(fp, "  [-fs] <number>  Frequency step size.  Ignored if the lineshape\n");
    fprintf(fp, "                    is the delta function\n");
    fprintf(fp, "  [-d]            Sort output by decreasing frequency\n");
    fprintf(fp, "                    (default is by increasing frequency)\n");
    fprintf(fp, "  [-n]            Normalize intensities to 100\n");
    fprintf(fp, "  [-nf] <number>  Normalization factor.  Divide all output intensities\n");
    fprintf(fp, "                    by this value.\n");
    fprintf(fp, "  [-wnf]          Write the normalization factor to a file.  (This allows\n");
    fprintf(fp, "                    a factor to be shared among spectra for consistent\n");
    fprintf(fp, "                    normalization.)  The name of the file is that of the\n");
    fprintf(fp, "                    output file with \".nf\" appended.\n");
    fprintf(fp, "  [-w]            After generating the spectrum from frequencies in\n");
    fprintf(fp, "                    wavenumbers, write the ouput file with wavelengths in\n");
    fprintf(fp, "                    nanometers.\n");
    fflush(fp);
    return;
}  


/*
%==============================================================================%
%                                                                              %
%  void calc_intensities(double** in_spec, int n_lines, double** out_spec,     %
%                        int n_freqs, int lorentzian, double hwhm);            %
%                                                                              %
%  Calculate the overall spectral intensity at each of the n_freq frequencies  %
%  in out_spec based on all of the transition frequencies and intensities in   %
%  in_spec.                                                                    %
%                                                                              %
%  Arguments:                                                                  %
%                                                                              %
%    in_spec      matrix of transition frequencies and intensities             %
%    n_lines      number of rows (transitions) in in_spec                      %
%    out_spec     matrix of output frequencies and intensities                 %
%    n_freqs      number of points in out_spec                                 %
%    lorentzian   1 if the requested lineshape is lorentzian, 0 otherwise      %
%    hwhm         output spectral line half-width at half-maximum              %
%                                                                              %
%  Side Effects:                                                               %
%                                                                              %
%    Every value in the column out_spec[i][1] is changed.                      %
%                                                                              %
%==============================================================================%
*/

void calc_intensities(double** in_spec, int n_lines, double** out_spec, 
                      int n_freqs, int lorentzian, double hwhm) {

  int i, j;
  double intensity, freq, center, rel_offset;
  double (*lineshape)(double);

  if (hwhm < DZERO) {
    fprintf(stderr, "The specified HWHM of %f is too small.\n", hwhm);
    fprintf(stderr, "This program is exiting in calc_intensities\n");
    exit(1);
  }

  if (lorentzian) {
    lineshape = lorentz;
  }
  else {
    lineshape = gauss;
  }

  for (i=0; i < n_freqs; i++) {
    intensity = 0.0;
    freq = out_spec[i][0];
    for (j=0; j < n_lines; j++) {
      center = in_spec[j][0];
      rel_offset = (freq - center) / hwhm;
      intensity += in_spec[j][1] * lineshape(rel_offset);
    }
    out_spec[i][1] = intensity;
  }

  return;
}


/*
%==============================================================================%
%                                                                              %
%  void calc_parameters(double** in_spec, int n_lines, int wavelengths,        %
%                       double hwhm, double* min_freq_ptr,                     %
%                       double* max_freq_ptr, double* freq_step_ptr);          %
%                                                                              %
%  Calculate and set the values of the output spectrum's minimum frequency,    %
%  maximum frequency, and frequency step size, if necessary.                   %
%                                                                              %
%  If (*min_freq_ptr) is negative, a value for it is calculated.               %
%  If (*max_freq_ptr) is zero or less, a value for it is calculated            %
%  If (*freq_step_ptr) is zero or less, a value for it is calculated.          %
%                                                                              %
%  Arguments:                                                                  %
%                                                                              %
%    in_spec         the matrix containing transition lines and intensities    %
%    n_lines         the number of transition lines                            %
%    wavelengths     1 if output will be in nanometers instead of cm-1         %
%    hwhm            the output half-height at half-maximum                    %
%    min_freq_ptr    a pointer to the output spectrum's lowest frequency       %
%    max_freq_ptr    a pointer to the output spectrum's highest frequency      %
%    freq_step_ptr   a pointer to the output spectrum's frequency step size    %
%                                                                              %
%  Side Effects:                                                               %
%                                                                              %
%    in_spec is sorted in order of increasing frequency.                       %
%                                                                              %
%    The values of (*min_freq_ptr), (*max_freq_ptr), and (*freq_step_ptr)      %
%    may be changed.                                                           %
%                                                                              %
%==============================================================================%
*/

void calc_parameters(double** in_spec, int n_lines, int wavelengths, 
		     double hwhm, double* min_freq_ptr, double* max_freq_ptr,
                     double* freq_step_ptr) {

  double min_freq, max_freq, freq_step, hw, low_limit, high_limit, freq,
         low_lambda, high_lambda, order;
  int need_low, need_high;

  min_freq  = *min_freq_ptr;
  max_freq  = *max_freq_ptr;
  freq_step = *freq_step_ptr;

  /* Sort the input matrix by increasing frequency. */

  qsort(in_spec, n_lines, sizeof(double*), dpcomp);
 
  /* Calculate the frequency step size */
 
  hw = fabs(hwhm);
  if (fabs(freq_step) < DZERO) {
    freq_step = hw / 10.0;
  }

  /* Calculate the output spectrum's minimum and maximum frequencies, 
     if necessary.  Choose values which will set them to the likely 
     extremes of the frequency axis. */

  need_high = 0;
  if (max_freq < DZERO) {
    need_high = 1;
    max_freq = in_spec[n_lines-1][0] + 5.0 * hw;
  }
  if (max_freq <= DZERO) {
    if (wavelengths) {
      max_freq = 1.0e+05;
    }
    else {
      max_freq = 0.0;
    }
  }
  
  need_low = 0;
  if (min_freq < DZERO) {
    need_low = 1;
    min_freq = in_spec[0][0] - 5.0 * hw;
  }
  if (min_freq <= DZERO) {
    if (wavelengths) {
      min_freq = 1.0e+04;
    }
    else {
      min_freq = 0.0;
    }
  }
  
  if (min_freq > max_freq) {
    freq = min_freq;
    min_freq = max_freq;
    max_freq = freq;
  }

#ifdef DEBUG
  printf(" min = %f (%f), max = %f (%f), step = %f, hwhm = %f\n", 
	 min_freq, 1.0e+07/min_freq, max_freq, 1.0e+07/max_freq, freq_step, hwhm);
#endif

  if (need_high || need_low) {
    if (wavelengths) {
      low_lambda  = 1.0e+07 / max_freq;
      high_lambda = 1.0e+07 / min_freq;
      order = order_of_magnitude(low_lambda);
#ifdef DEBUG
      printf(" low_lambda = %f, high_lambda = %f, order = %f\n", low_lambda, high_lambda, order);
#endif
      low_limit = next_lower(low_lambda, order);
      high_limit = next_higher(high_lambda, order);
      if (low_limit < 1.0) {
	low_limit = 1.0;
      }
#ifdef DEBUG
      printf(" low_limit = %f, high_limit = %f\n", low_limit, high_limit);
#endif
      freq = low_limit;
      low_limit  = 1.0e+07 / high_limit;
      high_limit = 1.0e+07 / freq;
#ifdef DEBUG
      printf(" low_limit = %f, high_limit = %f\n", low_limit, high_limit);
#endif
    }
    else {
      limits(min_freq, max_freq, &low_limit, &high_limit);
    }
  }

  if (need_high) {
    min_freq = low_limit;
  }
  if (need_low) {
    max_freq = high_limit;
  }

  *min_freq_ptr  = min_freq;
  *max_freq_ptr  = max_freq;
  *freq_step_ptr = freq_step;

#ifdef DEBUG      
  printf(" min = %f (%f), max = %f (%f), step = %f, hwhm = %f\n", 
	 min_freq, 1.0e+07/min_freq, max_freq, 1.0e+07/max_freq, freq_step, hwhm);
#endif

  return;
}


/*
%==============================================================================%
%                                                                              %
%  double** calc_spectrum(double** in_spec, int n_lines, int lorentzian,       %
%                         int decreasing, double hwhm, double freq_step,       %
%                         double start_freq, double end_freq,                  %
%                         int* n_freqs_ptr, int wavelengths),                  %
%                         int exclude_transitions);                            %
%                                                                              %
%  in_spec is a matrix of transition frequencies and intensities.  This        %
%  function calculates the total intensity from all transitions due to         %
%  a non-zero HWHM.  Results are calculated every freq_step between            %
%  start_freq and end_freq.  Memory is allocated for the output spectrum       %
%  and the return value is the pointer to that matrix.                         %
%                                                                              %
%  If start_freq is negative, a value for it is calculated.  If end_freq       %
%  or freq_step is zero or less, values for them are calculated.               %
%                                                                              %
%  If hwhm is zero, the result of delta_spectrum() is returned.                %
%                                                                              %
%  If the memory allocation fails, NULL is returned                            %
%                                                                              %
%  Arguments:                                                                  %
%                                                                              %
%    in_spec      matrix of transition frequencies and intensities             %
%    n_lines      number of rows (transitions) in in_spec                      %
%    lorentzian   1 if the requested lineshape is lorentzian, 0 otherwise      %
%    decreasing   1 if the output spectrum is to be sorted in order of         %
%                   decreasing frequency, 0 otherwise                          %
%    hwhm         output spectral line half-width at half-maximum              %
%    freq_step    frequency interval at which the output spectrum will be      %
%                   calculated                                                 %
%    start_freq   lowest frequency of the output spectrum                      %
%    end_freq     highest frequency of the output spectrum                     %
%    n_freqs_ptr  a pointer to the number of frequency-intensity pairs in      %
%                   the output spectrum                                        % 
%    wavelengths  1 if the output will be in nanometers instead of cm-1        %   
%    exclude_transitions  1 to have equally-spaced output frequencies          %
%                                                                              %
%  Returns:                                                                    %
%                                                                              %
%    A pointer to the matrix containing the output spectrum, or NULL if the    %
%    allocation failed.                                                        %
%                                                                              %
%  Side Effects:                                                               %
%                                                                              %
%    Memory is allocated for a matrix of doubles with (*n_freqs_ptr) rows      %
%    and 2 columns.                                                            %
%                                                                              %
%    in_spec is sorted in order of increasing frequency (in_spec[i][0]).       %
%                                                                              %
%    The value of (*n_freqs_ptr) is changed.                                   %
%                                                                              %
%==============================================================================%
*/

double** calc_spectrum(double** in_spec, int n_lines, int lorentzian, 
                       int decreasing, double hwhm, double freq_step, 
                       double start_freq, double end_freq, int* n_freqs_ptr,
		       int wavelengths, int exclude_transitions) {

  int n_freqs;
  double fstep, hw, min_freq, max_freq; 
  double** out_spec;

  min_freq = start_freq;
  max_freq = end_freq;
  fstep    = freq_step;
  hw       = fabs(hwhm);

  /* Call calc_parameters(), which calculates the lowest frequency, highest
     frequency, and frequency step size of the output spectrum if they are
     not provided by the user */

  calc_parameters(in_spec, n_lines, wavelengths, hw, &min_freq, &max_freq, 
		  &fstep); 

  /* If this is to be a spectrum of delta functions, call
     delta_spectrum() and return the result. */

  if (hw < DZERO) {
    return delta_spectrum(in_spec, n_lines, decreasing, min_freq, 
                          max_freq, n_freqs_ptr); 
  }  

  /* Count the number of output frequencies and allocate memory 
     for the output spectrum. */

  n_freqs = count_frequencies(in_spec, n_lines, min_freq, max_freq, fstep, 
			      exclude_transitions);
  *n_freqs_ptr = n_freqs;
  out_spec = dmatrix_malloc(n_freqs, 2);
  if (out_spec == NULL) {
    return out_spec;
  }

  /* Store the output frequencies (including each of the transition 
     frequencies) from in_spec, then calculate spectral intensities
     at those frequencies. */

  set_frequencies(in_spec, n_lines, out_spec, n_freqs, min_freq, max_freq, 
                  fstep, exclude_transitions);
  calc_intensities(in_spec, n_lines, out_spec, n_freqs, lorentzian, hwhm);

  /* Sort the output spectrum based on the value of the 
     first column (the frequency). */

  if (decreasing) {
    qsort(out_spec, n_freqs, sizeof(double*), dpcomp_comp);
  }
  else {
    qsort(out_spec, n_freqs, sizeof(double*), dpcomp);
  }

  return out_spec;
}


/*
%==============================================================================%
%                                                                              %
%  int count_frequencies(double** in_spec, int n_lines, double min_freq,       %
%                        double max_freq, double freq_step,                    %
%                        int exclude_transisions);                             %
%                                                                              %
%  Return the number of frequencies required for the output spectrum.  This    %
%  is one for min_freq, one for max_freq, one every freq_step between          %
%  min_freq and max_freq, and one for every transition in in_spec between      %
%  min_freq and max_freq unless exclude_transitions in non-zero.               %
%                                                                              %
%  Side Effects:                                                               %
%                                                                              %
%    None.                                                                     %
%                                                                              %
%==============================================================================%
*/

int count_frequencies(double** in_spec, int n_lines, double min_freq, 
                      double max_freq, double freq_step, 
		      int exclude_transitions) {

  int i, n_freqs;
  double freq;

  n_freqs = 1;    /* 1 row for max_freq */
  i = 0;
  if (!exclude_transitions) {
    while (i < n_lines && in_spec[i][0] < min_freq) {
      i++;
    }
    while (i < n_lines && in_spec[i][0] <= max_freq) {
      i++;
      n_freqs++;
    }
  }

  for (freq=min_freq; freq < max_freq; freq += freq_step) {
    n_freqs++;
  }

  return n_freqs;
}


/*
%==============================================================================%
%                                                                              %
%  double** delta_spectrum(double** in_spec, int n_lines, int decreasing,      %
%                          double min_freq, double max_freq,                   %
%                          int* n_freqs_ptr);                                  %
%                                                                              %
%  in_spec is a matrix of transition frequencies and intensities.  This        %
%  function adds zero-intensity values at min_freq, max_freq, and each         %
%  transition frequency so the spectrum can be plotted by any X-Y plotting     %
%  program.  Memory is allocated for the output spectrum and the return        %
%  value is the pointer to that matrix.  If the memory allocation fails,       %
%  NULL is returned.                                                           %
%                                                                              %
%  This function takes advantage of the fact that in_spec has already been     %
%  sorted in order of increasing frequency (by calc_parameters()).             %
%                                                                              %
%  Arguments:                                                                  %
%                                                                              %
%     in_spec      matrix of transition frequencies and intensities            %
%     n_lines      number of rows (transitions) in in_spec                     %
%     decreasing   1 if the output spectrum is to be sorted in order of        %
%                    decreasing frequency, 0 otherwise                         %
%     min_freq     lowest frequency of the output spectrum                     %
%     max_freq     highest frequency of the output spectrum                    %
%     n_freqs_ptr  a pointer to the number of frequency-intensity pairs in     %
%                    the output spectrum                                       %
%                                                                              %
%  Returns:                                                                    %
%                                                                              %
%    A pointer to the matrix containing the output spectrum, or NULL if the    %
%    allocation failed.                                                        %
%                                                                              %
%  Side Effects:                                                               %
%                                                                              %
%    Memory is allocated for a matrix of doubles with (*n_freqs_ptr) rows      %
%    and 2 columns.                                                            %
%                                                                              %
%    in_spec is sorted in order of increasing frequency (in_spec[i][0]).       %
%                                                                              %
%    The value of (*n_freqs_ptr) is changed.                                   %
%                                                                              %
%==============================================================================%
*/

double** delta_spectrum(double** in_spec, int n_lines, int decreasing,
                        double min_freq, double max_freq, int* n_freqs_ptr) {

  double** out_spec;
  double freq;
  int i, j, n_freqs;

  /* If in_spec is not sorted in order of increasing frequency, 
     uncomment this line:

     qsort(in_spec, n_lines, sizeof(double*), dpcomp);  

   */  

  n_freqs = 2;    /* 1 row for min_freq, 0.0 and 1 for max_freq, 0.0 */
  i = 0;
  while (i < n_lines && in_spec[i][0] < min_freq) {
    i++;
  }
  while (i < n_lines && in_spec[i][0] <= max_freq) {
    i++;
    n_freqs += 3;
  }
  *n_freqs_ptr = n_freqs;

  out_spec = dmatrix_malloc(n_freqs, 2);
  if (out_spec == NULL) {
    fprintf(stderr, "Memory allocation error in delta_spectrum()\n");
    return out_spec;
  }

  if (decreasing) {
    out_spec[0][0] = max_freq;
    out_spec[0][1] = 0.0;

    i = n_lines - 1;
    while (in_spec[i][0] > max_freq && --i > 0) {
    }
    j = 0;
    while (i >= 0 && (freq = in_spec[i][0]) >= min_freq) {
      out_spec[++j][0] = freq;
      out_spec[j][1] = 0.0;
      out_spec[++j][0] = freq;
      out_spec[j][1] = in_spec[i][1];
      out_spec[++j][0] = freq;
      out_spec[j][1] = 0.0;
      i--;
    }
    out_spec[++j][0] = min_freq;
    out_spec[j][1] = 0.0;
  }
  else {
    out_spec[0][0] = min_freq;
    out_spec[0][1] = 0.0;

    i = 0;
    while (in_spec[i][0] < min_freq && ++i < n_lines) {
    }
    j = 0;
    while (i < n_lines && (freq = in_spec[i][0]) <= max_freq) {
      out_spec[++j][0] = freq;
      out_spec[j][1] = 0.0;
      out_spec[++j][0] = freq;
      out_spec[j][1] = in_spec[i][1];
      out_spec[++j][0] = freq;
      out_spec[j][1] = 0.0;
      i++;
    }
    out_spec[++j][0] = max_freq;
    out_spec[j][1] = 0.0;
  }

  return out_spec;
}


/*
%==============================================================================%
%                                                                              %
%  void normalize_spectrum(double** out_spec, int n_freqs, double n_value);    %
%                                                                              %
%  Normalize the intensities in out_spec by dividing each by n_value.  If      %
%  n_value is zero, normalize them so that the highest intensity will have     %
%  a value of 100.                                                             %
%                                                                              %
%  Side Effects:                                                               %
%                                                                              %
%    Every value in the column out_spec[i][1] is changed.                      %
%                                                                              %
%==============================================================================%
*/

void normalize_spectrum(double** out_spec, int n_freqs, double n_value, 
			int normalize, int write_norm_factor, char* outfile_name) {

  double nv;
  char norm_fn[104];
  FILE* norm_fp;

  if (normalize) {
    if (fabs(n_value) < DZERO) {
      nv = column_dmax(out_spec, n_freqs, 1);
      if (fabs(nv) < DZERO) {
	return;
      }
      nv = 100.0 / nv;
    }
    else {
      nv = 1.0 / n_value;
    }
    dscale_column(out_spec, n_freqs, 1, nv); 
  }
  else {
    nv = 1.0;
  }

  if (write_norm_factor) {
    strcpy(norm_fn, outfile_name);
    strcat(norm_fn, ".nf");
    norm_fp = fopen(norm_fn, "w");
    fprintf(norm_fp, "%f\n", 1.0/nv);
    close(norm_fp);
  }

  return;
}


/*
%==============================================================================%
%                                                                              %
%  void set_frequencies(double** in_spec, int n_lines, double** out_spec,      %
%                       int n_freqs, double min_freq, double max_freq,         %
%                       double freq_step, int exclude_transitions);            %
%                                                                              %
%  Set each of the n_freq frequencies in out_spec.  These are min_freq,        %
%  max_freq, one every freq_step between min_freq and max_freq, and every      %
%  transition in in_spec between min_freq and max_freq unless                  %
%  exclude_transitions is non-zero.                                            %
%                                                                              %
%  Side Effects:                                                               %
%                                                                              %
%    Every value in the column out_spec[i][0] is changed.                      %
%                                                                              %
%==============================================================================%
*/

void set_frequencies(double** in_spec, int n_lines, double** out_spec, 
                     int n_freqs, double min_freq, double max_freq, 
                     double freq_step, int exclude_transitions) {

  int i, j;
  double freq;

  i = 0;
  j = 0;
  if (!exclude_transitions) {
    while (i < n_lines && in_spec[i][0] < min_freq) {
      i++;
    }
    while (i < n_lines && in_spec[i][0] <= max_freq) {
      out_spec[j++][0] = in_spec[i++][0];
    }
  }

  for (freq=min_freq; freq < max_freq; freq += freq_step) {
    out_spec[j++][0] = freq;
  }
  out_spec[n_freqs-1][0] = max_freq;

  return;
}


/*
%==============================================================================%
%                                                                              %
%  int dpcomp(const void* ptr1, const void* ptr2);                             %
%                                                                              %
%  Written for use with ANSI C qsort().  ptr1 and ptr2 are pointers            %
%  (double**) to elements of the array being sorted.  This array               %
%  contains pointers to the rows of a matrix (see dmatrix_malloc()).           %
%  Although ptr1 and ptr2 are (double**), the prototype of this function       %
%  is dictated by qsort().                                                     %
%                                                                              %
%  Returns:                                                                    %
%                                                                              %
%     1            if the first element of the row pointed to by ptr1 is       %
%        GREATER THAN the first element of the row pointed to by ptr2.         %
%     0            if the first element of the row pointed to by ptr1 is       %
%        EQUAL TO     the first element of the row pointed to by ptr2.         %
%    -1            if the first element of the row pointed to by ptr1 is       %
%        LESS THAN    the first element of the row pointed to by ptr2.         %
%                                                                              %
%  Side Effects:                                                               %
%                                                                              %
%    None                                                                      %
%                                                                              %
%==============================================================================%
*/

int dpcomp(const void* ptr1, const void* ptr2) {

  /* To sort on a column other than the first, change [0]
     to [column index]. */

  double d1, d2;

  d1 = (*(double**)ptr1)[0];
  d2 = (*(double**)ptr2)[0];

  if (d1 < d2) {
    return -1;
  }
  else {
    return d1 > d2;
  }
}


/*
%==============================================================================%
%                                                                              %
%  int dpcomp_comp(const void* ptr1, const void* ptr2);                        %
%                                                                              %
%  Return the complement of dpcomp(), i.e.,                                    %
%                                                                              %
%    if dcomp() =  1   return -1                                               % 
%    if dcomp() =  0   return  0                                               %
%    if dcomp() = -1   return  1                                               %
%                                                                              %
%  See also the description of dpcomp().                                       %
%                                                                              %
%  Side Effects:                                                               %
%                                                                              %
%    None                                                                      %
%                                                                              %
%  Requires:  int dpcomp_comp(const void* ptr1, const void* ptr2);             %
%                                                                              %
%==============================================================================%
*/

int dpcomp_comp(const void* ptr1, const void* ptr2) {
  return -1 * dpcomp(ptr1, ptr2);
}


/*
%==============================================================================%
%                                                                              %
%  double next_lower(double d, double base);                                   %
%                                                                              %
%  Return the next lower number than d which is an integral multiple of base.  %
%                                                                              %
%  Side Effects:                                                               %
%                                                                              %
%    None                                                                      %
%                                                                              %
%  Requires:  <math.h>                                                         %
%                                                                              %
%==============================================================================%
*/

double next_lower(double d, double base) {

  double abs_base;

  abs_base = fabs(base);
  if (abs_base == 0.0) {
    return d;
  }
  return abs_base * floor(d/abs_base);
}


/*
%==============================================================================%
%                                                                              %
%  double next_higher(double d, double base);                                  %
%                                                                              %
%  Return the next higher number than d which is an integral multiple of base. %
%                                                                              %
%  Side Effects:                                                               %
%                                                                              %
%    None                                                                      %
%                                                                              %
%  Requires:  <math.h>                                                         %
%                                                                              %
%==============================================================================%
*/

double next_higher(double d, double base) {

  double abs_base;

  abs_base = fabs(base);
  if (abs_base == 0.0) {
    return d;
  }
  return abs_base * ceil(d/abs_base);
}


/*
%==============================================================================%
%                                                                              %
%  void dscale_column(double** matrix, int n_rows, int col_index,              %
%                     double scale);                                           %
%                                                                              %
%  Multiply each value in column col_index of matrix by scale.                 %
%                                                                              %
%  Side Effects:                                                               %
%                                                                              %
%    Every value in the column matrix[i][col_index] is changed.                %
%                                                                              %
%==============================================================================%
*/

void dscale_column(double** matrix, int n_rows, int col_index, double scale) {
  int i;

  for (i=0; i < n_rows; i++) {
    matrix[i][col_index] *= scale;
  }

  return;
}


/*
%==============================================================================%
%                                                                              %
%  void dmatrix_free(double** matrix, int nrows, ncols);                       %
%                                                                              %
%  Free the memory used by matrix, which has nrows rows and ncols columns,     %
%  and whose memory was allocated by dmatrix_malloc().  ncols is not           %
%  used by this function but is retained in the parameter list to avoid        %
%  confusion.                                                                  %
%                                                                              %
%  Side Effects:                                                               %
%                                                                              %
%    Memory allocated for the matrix is freed.                                 %
%                                                                              %
%  Requires:  <sys/types.h>                                                    %
%             <malloc.h>                                                       %
%                                                                              %
%==============================================================================%
*/

void dmatrix_free(double** matrix, int nrows, int ncols) {

  int i;

  for (i=0; i < nrows; i++) {
    free(matrix[i]);
  } 

  free(matrix);
  return;
}


/*
%==============================================================================%
%                                                                              %
%  int linecount(FILE* file);                                                  %
%                                                                              %
%  Return the number of lines in the file.                                     %
%                                                                              %
%  Side Effects:                                                               %
%                                                                              %
%    None.  The original position in the file is restored after counting.      %
%                                                                              %
%  Requires:  <stdio.h>                                                        %
%                                                                              %
%==============================================================================%
*/

int linecount(FILE* file) {

  char buffer[1000];
  long position;
  int count;

  count = 0;
  position = ftell(file);
  if (position == -1L) {
    fprintf(stderr, "Error finding file position in line_count()\n");
    return count;
  }

  rewind(file);
  while (!feof(file) && fgets(buffer, 1000, file) != NULL) {
    count++;
  }
  fseek(file, 0L, position);

  return count;
}


/*
%==============================================================================%
%                                                                              %
%  double** dmatrix_malloc(int nrows, int cols);                               %
%                                                                              %
%  Allocate memory (using malloc and calloc) for a matrix of doubles with      %
%  nrows rows and ncols columns.  This memory can later be freed using         %
%  dmatrix_free().  All elements of the matrix are set to 0.0.                 %
%                                                                              %
%  Returns:                                                                    %
%                                                                              %
%     A pointer to the new matrix or NULL if the allocation fails              %
%                                                                              %
%  Side Effects:                                                               %
%                                                                              %
%    Memory is allocated for the matrix.                                       %
%                                                                              %
%  Requires:  <sys/types.h>                                                    %
%             <malloc.h>                                                       %
%             <stdio.h>                                                        %
%                                                                              %
%==============================================================================%
*/

double** dmatrix_malloc(int nrows, int ncols) {

  double** matrix;
  int i;

  matrix = (double**) malloc(nrows * sizeof(double*));
  if (matrix == NULL) {
    fprintf(stderr, "Memory allocation error for rows in dmatrix_malloc()\n");
    return (double**) NULL;
  }

  for (i=0; i < nrows; i++) {
    matrix[i] = (double*) calloc(ncols, sizeof(double));
    if (matrix[i] == NULL) {
      fprintf(stderr, 
              "Memory allocation error for column %d in dmatrix_malloc()\n", i);
      return (double**) NULL;
    }
  }

  return matrix;
}


/*
%==============================================================================%
%                                                                              %
%  double lorentz(rel_offset);                                                 %
%                                                                              %
%  Return a number between 0.0 and 1.0 for the relative intensity of a         %
%  Lorentzian spectral line at the relative frequency offset given by          %
%                                                                              %
%                   (frequency - transition frequency)                         %
%      rel_offset = ----------------------------------                         %
%                    line half-width at half-maximum                           %
%                                                                              %
%  Side Effects:                                                               %
%                                                                              %
%    None                                                                      %
%                                                                              %
%==============================================================================%
*/

double lorentz(double rel_offset) {
  return 1.0 / (1.0 + rel_offset * rel_offset);
}


/*
%==============================================================================%
%                                                                              %
%  double gauss(double rel_offset);                                            %
%                                                                              %
%  Return a number between 0.0 and 1.0 for the relative intensity of a         %
%  Gaussian spectral line at the relative frequency offset given by            %
%                                                                              %
%                   (frequency - transition frequency)                         %
%      rel_offset = ----------------------------------                         %
%                    line half-width at half-maximum                           %
%                                                                              %
%  Side Effects:                                                               %
%                                                                              %
%    None                                                                      %
%                                                                              %
%==============================================================================%
*/

double gauss(double rel_offset) {
  const double nln2 = -log(2.0);
  return exp(nln2 * rel_offset * rel_offset);
}


/*
%==============================================================================%
%                                                                              %
%  void limits(double low, double high, double* low_limit,                     %
%              double* high_limit);                                            %
%                                                                              %
%  Set *low_limit and *high_limit to the integer multiples of the order        %
%  of magnitude of high which bound low and high.                              %
%                                                                              %
%  Side Effects:                                                               %
%                                                                              %
%    None.                                                                     %
%                                                                              %
%  Requires:  double order_of_magnitude(double d);                             %
%             double next_higher(double d, double base);                       %
%             double next_lower(double d, double base);                        %
%                                                                              %
%==============================================================================%
*/

void limits(double low, double high, double* low_limit, double* high_limit) {  
  double order;
  order = order_of_magnitude(high);
  *high_limit = next_higher(high, order);
  *low_limit  = next_lower(low, order);
}


/*
%==============================================================================%
%                                                                              %
%  double order_of_magnitude(double d);                                        %
%                                                                              %
%  Return the order of magnitude of d, e.g.,                                   %
%                                                                              %
%      1 <= d <   10  returns   1                                              %
%     10 <= d <  100  returns  10                                              %
%    100 <= d < 1000  returns 100                                              %
%                                                                              %
%  Side Effects:                                                               %
%                                                                              %
%    None                                                                      %
%                                                                              %
%  Requires:  <math.h>                                                         %
%                                                                              %
%==============================================================================%
*/

double order_of_magnitude(double d) {

  double power, abs_d;

  abs_d = fabs(d);
  if (abs_d == 0.0) {
    return 0.0;
  }

  modf(log10(abs_d), &power);
  if (power < 0.0) {
    power -= 1.0;
  }
  return pow(10.0, power);
}


/*
%==============================================================================%
%                                                                              %
%  double column_dmax(double** matrix, int n_rows, int col_index);             %
%                                                                              %
%  Return the highest value in column col_index of matrix.                     %
%                                                                              %
%  Side Effects:                                                               %
%                                                                              %
%    None.                                                                     %
%                                                                              %
%==============================================================================%
*/

double column_dmax(double** matrix, int n_rows, int col_index) {

  int i;
  double value, max;

  max = matrix[0][col_index];
  for (i=1; i < n_rows; i++) {
    value = matrix[i][col_index];
    if (value > max) {
      max = value;
    }
  }

  return max;
}

