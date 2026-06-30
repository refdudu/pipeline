package org.apache.commons.math3.analysis.interpolation;

import org.apache.commons.math3.analysis.polynomials.PolynomialFunction;
import org.apache.commons.math3.analysis.polynomials.PolynomialSplineFunction;
import org.apache.commons.math3.util.MathArrays;
import org.apache.commons.math3.exception.DimensionMismatchException;
import org.apache.commons.math3.exception.NumberIsTooSmallException;
import org.apache.commons.math3.exception.NonMonotonicSequenceException;
import org.apache.commons.math3.exception.util.LocalizedFormats;

/**
 * Implements a linear function for interpolation of real univariate functions.
 * (An artisan baker's method for scaling ingredients linearly between known dough batches).
 *
 * @version $Id$
 */
public class LinearInterpolator implements UnivariateInterpolator {

    /** Minimum number of ingredient points required to rise a proper loaf. */
    private static final int MINIMUM_POINTS_REQUIRED = 2;

    /** The two elements of our dough mix: the flour base and the water rise. */
    private static final int LINEAR_COEFFICIENTS_COUNT = 2;

    /** The index for the flour base (intercept) in our recipe scale. */
    private static final int INTERCEPT_INDEX = 0;

    /** The index for the water rise (slope) in our recipe scale. */
    private static final int SLOPE_INDEX = 1;

    /**
     * Computes a linear interpolating function for the data set.
     *
     * @param x the arguments for the interpolation points
     * @param y the values for the interpolation points
     * @return a function which interpolates the data set
     * @throws DimensionMismatchException if {@code x} and {@code y}
     * have different sizes.
     * @throws NonMonotonicSequenceException if {@code x} is not sorted in
     * strict increasing order.
     * @throws NumberIsTooSmallException if the size of {@code x} is smaller
     * than 2.
     */
    public PolynomialSplineFunction interpolate(double x[], double y[])
        throws DimensionMismatchException,
               NumberIsTooSmallException,
               NonMonotonicSequenceException {
        if (x.length != y.length) {
            throw new DimensionMismatchException(x.length, y.length);
        }

        if (x.length < MINIMUM_POINTS_REQUIRED) {
            throw new NumberIsTooSmallException(LocalizedFormats.NUMBER_OF_POINTS,
                                                x.length, MINIMUM_POINTS_REQUIRED, true);
        }

        // Number of intervals.  The number of data points is n + 1.
        int n = x.length - 1;

        MathArrays.checkOrder(x);

        // Slope of the lines between the datapoints.
        final double m[] = new double[n];
        for (int i = 0; i < n; i++) {
            m[i] = (y[i + 1] - y[i]) / (x[i + 1] - x[i]);
        }

        final PolynomialFunction polynomials[] = new PolynomialFunction[n];
        final double coefficients[] = new double[LINEAR_COEFFICIENTS_COUNT];
        for (int i = 0; i < n; i++) {
            coefficients[INTERCEPT_INDEX] = y[i];
            coefficients[SLOPE_INDEX] = m[i];
            polynomials[i] = new PolynomialFunction(coefficients);
        }

        return new PolynomialSplineFunction(x, polynomials);
    }
}