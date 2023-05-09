package com.example.myaudiofilter;

/**************************************************************************
 * Parks-McClellan algorithm for FIR filter design (C version)
 *-------------------------------------------------
 *  Copyright (c) 1995,1998  Jake Janovetz (janovetz@uiuc.edu)
 *
 *  This library is free software; you can redistribute it and/or
 *  modify it under the terms of the GNU Library General Public
 *  License as published by the Free Software Foundation; either
 *  version 2 of the License, or (at your option) any later version.
 *
 *  This library is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 *  Library General Public License for more details.
 *
 *  You should have received a copy of the GNU Library General Public
 *  License along with this library; if not, write to the Free
 *  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 *
 *************************************************************************/
//package com.ledin.firdesign;
/// <summary>
/// A static class containing a function (remez) that uses the Remez exchange
/// algorithm to design an FIR filter.
/// </summary>
public class Remez {
    // / <summary>
    // / Filter symmetry types
    // / </summary>
    enum SymmetryType {
        Negative, Positive
    };
    // / <summary>
    // / Density factor for creating the grid of frequencies from the specified
    // bands
    // / </summary>
    static final int GridDensity = 16;
    // / <summary>
    // / Maximum number of Remez exchange algorithm iterations
    // / </summary>
    final static int MaxIterations = 40;
    // IEEE-754 double precision epsilon
    public final static double Epsilon = Math.pow(2, -53);
    /*******************
     * CreateDenseGrid ================= Creates the dense grid of frequencies
     * from the specified bands. Also creates the Desired Frequency Response
     * function (D[]) and the Weight function (W[]) on that dense grid
     *
     *
     * INPUT: ------ int r - 1/2 the number of filter coefficients int numtaps -
     * Number of taps in the resulting filter int numband - Number of bands in
     * user specification double bands[] - User-specified band edges [2*numband]
     * double des[] - Desired response per band [numband] double weight[] -
     * Weight per band [numband] int symmetry - Symmetry of filter - used for
     * grid check
     *
     * OUTPUT: ------- int gridsize - Number of elements in the dense frequency
     * grid double Grid[] - Frequencies (0 to 0.5) on the dense grid [gridsize]
     * double D[] - Desired response on the dense grid [gridsize] double W[] -
     * Weight function on the dense grid [gridsize]
     *******************/
    static void CreateDenseGrid(int r, int numtaps, int numband,
                                double[] bands, double[] des, double[] weight, final int gridsize,
                                final double[] Grid, final double[] D, final double[] W,
                                SymmetryType symmetry) {
        int i, j, k, band;
        double delf, lowf, highf;
        delf = 0.5 / (GridDensity * r);
        /*
         * For differentiator, hilbert, symmetry is odd and Grid[0] = max(delf,
         * band[0])
         */
        if ((symmetry == SymmetryType.Negative) && (delf > bands[0]))
            bands[0] = delf;
        j = 0;
        for (band = 0; band < numband; band++) {
            System.out.println("j"+j);
            System.out.println("bands"+2*band);
            Grid[j] = bands[2 * band];
            lowf = bands[2 * band];
            highf = bands[2 * band + 1];
            k = (int) ((highf - lowf) / delf + 0.5); /* .5 for rounding */
            for (i = 0; i < k; i++) {
                D[j] = des[band];
                W[j] = weight[band];
                Grid[j] = lowf;
                lowf += delf;
                j++;
            }
            Grid[j - 1] = highf;
        }
        /*
         * Similar to above, if odd symmetry, last grid point can't be .5 - but,
         * if there are even taps, leave the last grid point at .5
         */
        if ((symmetry == SymmetryType.Negative)
                && (Grid[gridsize - 1] > (0.5 - delf)) && (numtaps % 2) != 0) {
            Grid[gridsize - 1] = 0.5 - delf;
        }
    }
    /********************
     * InitialGuess ============== Places Extremal Frequencies evenly throughout
     * the dense grid.
     *
     *
     * INPUT: ------ int r - 1/2 the number of filter coefficients int gridsize
     * - Number of elements in the dense frequency grid
     *
     * OUTPUT: ------- int Ext[] - Extremal indexes to dense frequency grid
     * [r+1]
     ********************/
    static void InitialGuess(int r, final int[] Ext, int gridsize) {
        int i;
        for (i = 0; i <= r; i++)
            Ext[i] = i * (gridsize - 1) / r;
    }
    /***********************
     * CalcParms ===========
     *
     *
     * INPUT: ------ int r - 1/2 the number of filter coefficients int Ext[] -
     * Extremal indexes to dense frequency grid [r+1] double Grid[] -
     * Frequencies (0 to 0.5) on the dense grid [gridsize] double D[] - Desired
     * response on the dense grid [gridsize] double W[] - Weight function on the
     * dense grid [gridsize]
     *
     * OUTPUT: ------- double ad[] - 'b' in Oppenheim & Schafer [r+1] double x[]
     * - [r+1] double y[] - 'C' in Oppenheim & Schafer [r+1]
     ***********************/
    static void CalcParms(int r, int[] Ext, double[] Grid, double[] D,
                          double[] W, final double[] ad, final double[] x, final double[] y) {
        int i, j, k, ld;
        double sign, xi, delta, denom, numer;
        /*
         * Find x[]
         */
        for (i = 0; i <= r; i++)
            x[i] = Math.cos(2 * Math.PI * Grid[Ext[i]]);
        /*
         * Calculate ad[] - Oppenheim & Schafer eq 7.132
         */
        ld = (r - 1) / 15 + 1; /* Skips around to avoid round errors */
        for (i = 0; i <= r; i++) {
            denom = 1.0;
            xi = x[i];
            for (j = 0; j < ld; j++) {
                for (k = j; k <= r; k += ld)
                    if (k != i)
                        denom *= 2.0 * (xi - x[k]);
            }
            if (Math.abs(denom) < 1e-9)
                denom = 1e-9;
            ad[i] = 1.0 / denom;
        }
        /*
         * Calculate delta - Oppenheim & Schafer eq 7.131
         */
        numer = denom = 0;
        sign = 1;
        for (i = 0; i <= r; i++) {
            numer += ad[i] * D[Ext[i]];
            denom += sign * ad[i] / W[Ext[i]];
            sign = -sign;
        }
        delta = numer / denom;
        sign = 1;
        /*
         * Calculate y[] - Oppenheim & Schafer eq 7.133b
         */
        for (i = 0; i <= r; i++) {
            y[i] = D[Ext[i]] - sign * delta / W[Ext[i]];
            sign = -sign;
        }
    }
    /*********************
     * ComputeA ========== Using values calculated in CalcParms, ComputeA
     * calculates the actual filter response at a given frequency (freq). Uses
     * eq 7.133a from Oppenheim & Schafer.
     *
     *
     * INPUT: ------ double freq - Frequency (0 to 0.5) at which to calculate A
     * int r - 1/2 the number of filter coefficients double ad[] - 'b' in
     * Oppenheim & Schafer [r+1] double x[] - [r+1] double y[] - 'C' in
     * Oppenheim & Schafer [r+1]
     *
     * OUTPUT: ------- Returns double value of A[freq]
     *********************/
    static double ComputeA(double freq, int r, double[] ad, double[] x,
                           double[] y) {
        int i;
        double xc, c, denom, numer;
        denom = numer = 0;
        xc = Math.cos(2 * Math.PI * freq);
        for (i = 0; i <= r; i++) {
            c = xc - x[i];
            if (c < 1.0e-7 && c > -1.0e-7) { // JAL 1/2/11 Changed from Math.abs
                // to speed it up
                numer = y[i];
                denom = 1;
                break;
            }
            c = ad[i] / c;
            denom += c;
            numer += c * y[i];
        }
        return numer / denom;
    }
    /************************
     * CalcError =========== Calculates the Error function from the desired
     * frequency response on the dense grid (D[]), the weight function on the
     * dense grid (W[]), and the present response calculation (A[])
     *
     *
     * INPUT: ------ int r - 1/2 the number of filter coefficients double ad[] -
     * [r+1] double x[] - [r+1] double y[] - [r+1] int gridsize - Number of
     * elements in the dense frequency grid double Grid[] - Frequencies on the
     * dense grid [gridsize] double D[] - Desired response on the dense grid
     * [gridsize] double W[] - Weight function on the desnse grid [gridsize]
     *
     * OUTPUT: ------- double E[] - Error function on dense grid [gridsize]
     ************************/
    static void CalcError(int r, double[] ad, double[] x, double[] y,
                          int gridsize, double[] Grid, double[] D, double[] W,
                          final double[] E) {
        int i;
        double A;
        for (i = 0; i < gridsize; i++) {
            A = ComputeA(Grid[i], r, ad, x, y);
            E[i] = W[i] * (D[i] - A);
        }
    }
    /************************
     * Search ======== Searches for the maxima/minima of the error curve. If
     * more than r+1 extrema are found, it uses the following heuristic (thanks
     * Chris Hanson): 1) Adjacent non-alternating extrema deleted first. 2) If
     * there are more than one excess extrema, delete the one with the smallest
     * error. This will create a non-alternation condition that is fixed by 1).
     * 3) If there is exactly one excess extremum, delete the smaller of the
     * first/last extremum
     *
     *
     * INPUT: ------ int r - 1/2 the number of filter coefficients int Ext[] -
     * Indexes to Grid[] of extremal frequencies [r+1] int gridsize - Number of
     * elements in the dense frequency grid double E[] - Array of error values.
     * [gridsize] OUTPUT: ------- int Ext[] - New indexes to extremal
     * frequencies [r+1]
     ************************/
    static void Search(int r, final int[] Ext, int gridsize, double[] E) {
        int i, j, k, l, extra; /* Counters */
        int up, alt;
        int[] foundExt = new int[10 * r]; /* Array of found extremals */
        k = 0;
        /*
         * Check for extremum at 0.
         */
        if (((E[0] > 0.0) && (E[0] > E[1])) || ((E[0] < 0.0) && (E[0] < E[1])))
            foundExt[k++] = 0;
        /*
         * Check for extrema inside dense grid
         */
        for (i = 1; i < gridsize - 1; i++) {
            if (((E[i] >= E[i - 1]) && (E[i] > E[i + 1]) && (E[i] > 0.0))
                    || ((E[i] <= E[i - 1]) && (E[i] < E[i + 1]) && (E[i] < 0.0)))
                foundExt[k++] = i;
        }
        /*
         * Check for extremum at 0.5
         */
        j = gridsize - 1;
        if (((E[j] > 0.0) && (E[j] > E[j - 1]))
                || ((E[j] < 0.0) && (E[j] < E[j - 1])))
            foundExt[k++] = j;
        /*
         * Remove extra extremals
         */
        extra = k - (r + 1);
        while (extra > 0) {
            if (E[foundExt[0]] > 0.0)
                up = 1; /* first one is a maxima */
            else
                up = 0; /* first one is a minima */
            l = 0;
            alt = 1;
            for (j = 1; j < k; j++) {
                if (Math.abs(E[foundExt[j]]) < Math.abs(E[foundExt[l]]))
                    l = j; /* new smallest error. */
                if ((up != 0) && (E[foundExt[j]] < 0.0))
                    up = 0; /* switch to a minima */
                else if ((up == 0) && (E[foundExt[j]] > 0.0))
                    up = 1; /* switch to a maxima */
                else {
                    alt = 0;
                    break; /* Ooops, found two non-alternating */
                } /* extrema. Delete smallest of them */
            } /* if the loop finishes, all extrema are alternating */
            /*
             * If there's only one extremal and all are alternating, delete the
             * smallest of the first/last extremals.
             */
            if ((alt != 0) && (extra == 1)) {
                if (Math.abs(E[foundExt[k - 1]]) < Math.abs(E[foundExt[0]]))
                    l = foundExt[k - 1]; /* Delete last extremal */
                else
                    l = foundExt[0]; /* Delete first extremal */
            }
            for (j = l; j < k; j++) /* Loop that does the deletion */
            {
                foundExt[j] = foundExt[j + 1];
            }
            k--;
            extra--;
        }
        for (i = 0; i <= r; i++) {
            Ext[i] = foundExt[i]; /* Copy found extremals to Ext[] */
        }
    }
    /*********************
     * FreqSample ============ Simple frequency sampling algorithm to determine
     * the impulse response h[] from A's found in ComputeA
     *
     *
     * INPUT: ------ int N - Number of filter coefficients double A[] - Sample
     * points of desired response [N/2] int symmetry - Symmetry of desired
     * filter
     *
     * OUTPUT: ------- double h[] - Impulse Response of final filter [N]
     *********************/
    static void FreqSample(int N, double[] A, SymmetryType symm,
                           final double[] h) {
        int n, k;
        double x, val, M;
        M = (N - 1.0) / 2.0;
        if (symm == SymmetryType.Positive) {
            if (N % 2 != 0) {
                for (n = 0; n < N; n++) {
                    val = A[0];
                    x = 2 * Math.PI * (n - M) / N;
                    for (k = 1; k <= M; k++)
                        val += 2.0 * A[k] * Math.cos(x * k);
                    h[n] = val / N;
                }
            } else {
                for (n = 0; n < N; n++) {
                    val = A[0];
                    x = 2 * Math.PI * (n - M) / N;
                    for (k = 1; k <= (N / 2 - 1); k++)
                        val += 2.0 * A[k] * Math.cos(x * k);
                    h[n] = val / N;
                }
            }
        } else {
            if (N % 2 != 0) {
                for (n = 0; n < N; n++) {
                    val = 0;
                    x = 2 * Math.PI * (n - M) / N;
                    for (k = 1; k <= M; k++)
                        val += 2.0 * A[k] * Math.sin(x * k);
                    h[n] = val / N;
                }
            } else {
                for (n = 0; n < N; n++) {
                    val = A[N / 2] * Math.sin(Math.PI * (n - M));
                    x = 2 * Math.PI * (n - M) / N;
                    for (k = 1; k <= (N / 2 - 1); k++)
                        val += 2.0 * A[k] * Math.sin(x * k);
                    h[n] = val / N;
                }
            }
        }
    }
    /*******************
     * isDone ======== Checks to see if the error function is small enough to
     * consider the result to have converged.
     *
     * INPUT: ------ int r - 1/2 the number of filter coeffiecients int Ext[] -
     * Indexes to extremal frequencies [r+1] double E[] - Error function on the
     * dense grid [gridsize]
     *
     * OUTPUT: ------- Returns 1 if the result converged Returns 0 if the result
     * has not converged
     ********************/
    static short isDone(int r, int[] Ext, double[] E) {
        int i;
        double min, max, current;
        min = max = Math.abs(E[Ext[0]]);
        for (i = 1; i <= r; i++) {
            current = Math.abs(E[Ext[i]]);
            if (current < min)
                min = current;
            if (current > max)
                max = current;
        }
        if (((max - min) / max) < 1e-9)
            return 1;
        return 0;
    }
    // / <summary>
    // / The type of filter to design
    // / </summary>
    public enum FilterType {
        Bandpass, Differentiator, Hilbert
    };
    // / <summary>
    // / Calculates the optimal (in the Chebyshev/minimax sense)
    // / FIR filter impulse response given a set of band edges,
    // / the desired reponse on those bands, and the weight given to
    // / the error in those bands.
    // / </summary>
    // / <param name="numtaps">Number of filter coefficients</param>
    // / <param name="numband">Number of bands in filter specification</param>
    // / <param name="bands">User-specified band edges [2 * numband]</param>
    // / <param name="des">User-specified band responses [numband]</param>
    // / <param name="weight">User-specified error weights [numband]</param>
    // / <param name="type">Type of filter</param>
    // / <param name="h">Impulse response of final filter [numtaps]</param>
    // / <param name="freqresp_freq">Frequency response frequencies</param>
    // / <param name="freqresp_mag">Frequency response magnitudes (dB)</param>
    // / <returns>Return value = true if success, false if failed to
    // converge</returns>
    public static boolean remez(int numtaps, int numband, double[] bands,
                                double[] des, double[] weight, FilterType type, final double[] x,
                                final double[] y, final double[] ad, final double[] h) {
        double[] Grid, W, D, E;
        int i, iter, gridsize;
        int[] Ext;
        double c;
        double[] taps;
        SymmetryType symmetry;
        if (type == FilterType.Bandpass)
            symmetry = SymmetryType.Positive;
        else
            symmetry = SymmetryType.Negative;
        int r = numtaps / 2; /* number of extrema */
        if ((numtaps % 2 != 0) && (symmetry == SymmetryType.Positive))
            r++;
        /*
         * Predict dense grid size in advance for memory allocation .5 is so we
         * round up, not truncate
         */
        gridsize = 0;
        for (i = 0; i < numband; i++)
            gridsize += (int) (2 * r * GridDensity
                    * (bands[2 * i + 1] - bands[2 * i]) + .5);
        if (symmetry == SymmetryType.Negative)
            gridsize--;
        /*
         * Dynamically allocate memory for arrays with proper sizes
         */
        Grid = new double[gridsize];
        D = new double[gridsize];
        W = new double[gridsize];
        E = new double[gridsize];
        Ext = new int[r + 1];
        taps = new double[r + 1];
        /*
         * Create dense frequency grid
         */
        CreateDenseGrid(r, numtaps, numband, bands, des, weight, gridsize,
                Grid, D, W, symmetry);
        InitialGuess(r, Ext, gridsize);
        /*
         * For Differentiator: (fix grid)
         */
        if (type == FilterType.Differentiator)
            for (i = 0; i < gridsize; i++) {
                /* D[i] = D[i]*Grid[i]; */
                if (D[i] > 1000 * Epsilon)
                    W[i] = W[i] / Grid[i];
            }
        /*
         * For odd or Negative symmetry filters, alter the D[] and W[] according
         * to Parks McClellan
         */
        if (symmetry == SymmetryType.Positive) {
            if (numtaps % 2 == 0) {
                for (i = 0; i < gridsize; i++) {
                    c = Math.cos(Math.PI * Grid[i]);
                    D[i] /= c;
                    W[i] *= c;
                }
            }
        } else {
            if (numtaps % 2 != 0) {
                for (i = 0; i < gridsize; i++) {
                    c = Math.sin(2 * Math.PI * Grid[i]);
                    D[i] /= c;
                    W[i] *= c;
                }
            } else {
                for (i = 0; i < gridsize; i++) {
                    c = Math.sin(Math.PI * Grid[i]);
                    D[i] /= c;
                    W[i] *= c;
                }
            }
        }
        /*
         * Perform the Remez Exchange algorithm
         */
        for (iter = 0; iter < MaxIterations; iter++) {
            CalcParms(r, Ext, Grid, D, W, ad, x, y);
            CalcError(r, ad, x, y, gridsize, Grid, D, W, E);
            Search(r, Ext, gridsize, E);
            if (isDone(r, Ext, E) != 0)
                break;
        }
        if (iter == MaxIterations)
            return false; /* Indicate failure to converge */
        CalcParms(r, Ext, Grid, D, W, ad, x, y);
        /*
         * Find the 'taps' of the filter for use with Frequency Sampling. If odd
         * or Negative symmetry, fix the taps according to Parks McClellan
         */
        for (i = 0; i <= numtaps / 2; i++) {
            if (symmetry == SymmetryType.Positive) {
                if (numtaps % 2 != 0)
                    c = 1;
                else
                    c = Math.cos(Math.PI * (double) i / numtaps);
            } else {
                if (numtaps % 2 != 0)
                    c = Math.sin(2 * Math.PI * (double) i / numtaps);
                else
                    c = Math.sin(Math.PI * (double) i / numtaps);
            }
            taps[i] = ComputeA((double) i / numtaps, r, ad, x, y) * c;
        }
        /*
         * Frequency sampling design with calculated taps
         */
        FreqSample(numtaps, taps, symmetry, h);
        /*
         * If this is a differentiator, signs of coefficients must be reversed
         */
        if (type == FilterType.Differentiator)
            for (i = 0; i < numtaps; i++)
                h[i] = -h[i];
        return true; /* Indicate success */
    }
    public static void computeFreqResponse(int numtaps, int freqresp_points,
                                           final double[] freqresp_freq, final double[] freqresp_mag,
                                           double[] x, double[] y, double[] ad, FilterType type) {
        int i;
        double c;
        double prev_resp = 0;
        SymmetryType symmetry;
        if (type == FilterType.Bandpass)
            symmetry = SymmetryType.Positive;
        else
            symmetry = SymmetryType.Negative;
        int r = numtaps / 2; /* number of extrema */
        if ((numtaps % 2 != 0) && (symmetry == SymmetryType.Positive))
            r++;
        for (i = 0; i < freqresp_points; i++) {
            double freq = i * 0.5 / (freqresp_points - 1);
            if (symmetry == SymmetryType.Positive) {
                if (numtaps % 2 != 0)
                    c = 1;
                else
                    c = Math.cos(Math.PI * freq);
            } else {
                if (numtaps % 2 != 0)
                    c = Math.sin(2 * Math.PI * freq);
                else
                    c = Math.sin(Math.PI * freq);
            }
            double resp = ComputeA(freq, r, ad, x, y) * c;
            // If response has crossed zero since the last point, set the
            // nearest point to zero
            if (i > 0) {
                if ((prev_resp * resp) < 0) {
                    if (Math.abs(prev_resp) < Math.abs(resp))
                        freqresp_mag[i - 1] = -300;
                    else
                        resp = 0;
                }
            }
            prev_resp = resp;
            double abs_resp;
            if (resp >= 0)
                abs_resp = resp;
            else
                abs_resp = -resp;
            double log_resp;
            if (abs_resp > 1e-15)
                log_resp = 20 * Math.log10(abs_resp);
            else
                log_resp = -300;
            freqresp_freq[i] = freq;
            freqresp_mag[i] = log_resp;
        }
    }
}
