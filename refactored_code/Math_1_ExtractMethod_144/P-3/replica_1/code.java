/*
 * Licensed to the Apache Software Foundation (ASF) under one or more
 * contributor license agreements.  See the NOTICE file distributed with
 * this work for additional information regarding copyright ownership.
 * The ASF licenses this file to You under the Apache License, Version 2.0
 * (the "License"); you may not use this file except in compliance with
 * the License.  You may obtain a copy of the License at
 *
 *      http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */
package org.apache.commons.math3.analysis.interpolation;

import org.apache.commons.math3.analysis.UnivariateFunction;
import org.apache.commons.math3.analysis.polynomials.PolynomialSplineFunction;
import org.apache.commons.math3.exception.DimensionMismatchException;
import org.apache.commons.math3.exception.NoDataException;
import org.apache.commons.math3.exception.NonMonotonicSequenceException;
import org.apache.commons.math3.exception.NumberIsTooSmallException;
import org.apache.commons.math3.util.MathArrays;

/**
 * Generates a bicubic interpolating function.
 *
 * @version $Id$
 * @since 2.2
 */
public class BicubicSplineInterpolator
    implements BivariateGridInterpolator {
    /**
     * {@inheritDoc}
     */
    public BicubicSplineInterpolatingFunction interpolate(final double[] xval,
                                                          final double[] yval,
                                                          final double[][] fval)
        throws NoDataException, DimensionMismatchException,
               NonMonotonicSequenceException, NumberIsTooSmallException {
        
        validateInput(xval, yval, fval);

        final int xLen = xval.length;
        final int yLen = yval.length;

        final double[][] fX = createFX(fval, xLen, yLen);

        final SplineInterpolator spInterpolator = new SplineInterpolator();

        final PolynomialSplineFunction[] ySplineX = interpolateYLines(spInterpolator, xval, fX, yLen);
        final PolynomialSplineFunction[] xSplineY = interpolateXLines(spInterpolator, yval, fval, xLen);

        final double[][] dFdX = computeDFdX(ySplineX, xval, yLen);
        final double[][] dFdY = computeDFdY(xSplineY, yval, xLen);
        final double[][] d2FdXdY = computeD2FdXdY(fval, xval, yval);

        return new BicubicSplineInterpolatingFunction(xval, yval, fval,
                                                      dFdX, dFdY, d2FdXdY);
    }

    /**
     * Validates input arrays.
     */
    private void validateInput(final double[] xval, final double[] yval, final double[][] fval)
        throws NoDataException, DimensionMismatchException, NonMonotonicSequenceException {
        if (xval.length == 0 || yval.length == 0 || fval.length == 0) {
            throw new NoDataException();
        }
        if (xval.length != fval.length) {
            throw new DimensionMismatchException(xval.length, fval.length);
        }

        MathArrays.checkOrder(xval);
        MathArrays.checkOrder(yval);
    }

    /**
     * Creates and returns the transposed array of samples.
     */
    private double[][] createFX(final double[][] fval, final int xLen, final int yLen)
        throws DimensionMismatchException {
        final double[][] fX = new double[yLen][xLen];
        for (int i = 0; i < xLen; i++) {
            if (fval[i].length != yLen) {
                throw new DimensionMismatchException(fval[i].length, yLen);
            }

            for (int j = 0; j < yLen; j++) {
                fX[j][i] = fval[i][j];
            }
        }
        return fX;
    }

    /**
     * For each line y[j], constructs a 1D spline with respect to variable x.
     */
    private PolynomialSplineFunction[] interpolateYLines(final SplineInterpolator interpolator,
                                                         final double[] xval,
                                                         final double[][] fX,
                                                         final int yLen) {
        final PolynomialSplineFunction[] ySplineX = new PolynomialSplineFunction[yLen];
        for (int j = 0; j < yLen; j++) {
            ySplineX[j] = interpolator.interpolate(xval, fX[j]);
        }
        return ySplineX;
    }

    /**
     * For each line x[i], constructs a 1D spline with respect to variable y.
     */
    private PolynomialSplineFunction[] interpolateXLines(final SplineInterpolator interpolator,
                                                         final double[] yval,
                                                         final double[][] fval,
                                                         final int xLen) {
        final PolynomialSplineFunction[] xSplineY = new PolynomialSplineFunction[xLen];
        for (int i = 0; i < xLen; i++) {
            xSplineY[i] = interpolator.interpolate(yval, fval[i]);
        }
        return xSplineY;
    }

    /**
     * Computes partial derivatives with respect to x at the grid knots.
     */
    private double[][] computeDFdX(final PolynomialSplineFunction[] ySplineX,
                                   final double[] xval,
                                   final int yLen) {
        final int xLen = xval.length;
        final double[][] dFdX = new double[xLen][yLen];
        for (int j = 0; j < yLen; j++) {
            final UnivariateFunction f = ySplineX[j].derivative();
            for (int i = 0; i < xLen; i++) {
                dFdX[i][j] = f.value(xval[i]);
            }
        }
        return dFdX;
    }

    /**
     * Computes partial derivatives with respect to y at the grid knots.
     */
    private double[][] computeDFdY(final PolynomialSplineFunction[] xSplineY,
                                   final double[] yval,
                                   final int xLen) {
        final int yLen = yval.length;
        final double[][] dFdY = new double[xLen][yLen];
        for (int i = 0; i < xLen; i++) {
            final UnivariateFunction f = xSplineY[i].derivative();
            for (int j = 0; j < yLen; j++) {
                dFdY[i][j] = f.value(yval[j]);
            }
        }
        return dFdY;
    }

    /**
     * Computes cross partial derivatives.
     */
    private double[][] computeD2FdXdY(final double[][] fval,
                                      final double[] xval,
                                      final double[] yval) {
        final int xLen = xval.length;
        final int yLen = yval.length;
        final double[][] d2FdXdY = new double[xLen][yLen];
        for (int i = 0; i < xLen; i++) {
            final int nI = nextIndex(i, xLen);
            final int pI = previousIndex(i);
            for (int j = 0; j < yLen; j++) {
                final int nJ = nextIndex(j, yLen);
                final int pJ = previousIndex(j);
                d2FdXdY[i][j] = (fval[nI][nJ] - fval[nI][pJ] -
                                 fval[pI][nJ] + fval[pI][pJ]) /
                    ((xval[nI] - xval[pI]) * (yval[nJ] - yval[pJ]));
            }
        }
        return d2FdXdY;
    }

    /**
     * Computes the next index of an array, clipping if necessary.
     * It is assumed (but not checked) that {@code i >= 0}.
     *
     * @param i Index.
     * @param max Upper limit of the array.
     * @return the next index.
     */
    private int nextIndex(int i, int max) {
        final int index = i + 1;
        return index < max ? index : index - 1;
    }

    /**
     * Computes the previous index of an array, clipping if necessary.
     * It is assumed (but not checked) that {@code i} is smaller than the size
     * of the array.
     *
     * @param i Index.
     * @return the previous index.
     */
    private int previousIndex(int i) {
        final int index = i - 1;
        return index >= 0 ? index : 0;
    }
}