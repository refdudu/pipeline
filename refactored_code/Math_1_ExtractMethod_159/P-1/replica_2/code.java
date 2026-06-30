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

import org.apache.commons.math3.exception.DimensionMismatchException;
import org.apache.commons.math3.exception.NoDataException;
import org.apache.commons.math3.exception.NonMonotonicSequenceException;
import org.apache.commons.math3.exception.NotPositiveException;
import org.apache.commons.math3.exception.NullArgumentException;
import org.apache.commons.math3.util.MathArrays;
import org.apache.commons.math3.util.Precision;
import org.apache.commons.math3.optim.nonlinear.vector.jacobian.GaussNewtonOptimizer;
import org.apache.commons.math3.fitting.PolynomialFitter;
import org.apache.commons.math3.analysis.polynomials.PolynomialFunction;
import org.apache.commons.math3.optim.SimpleVectorValueChecker;

/**
 * Generates a bicubic interpolation function.
 * Prior to generating the interpolating function, the input is smoothed using
 * polynomial fitting.
 *
 * @version $Id$
 * @since 2.2
 */
public class SmoothingPolynomialBicubicSplineInterpolator
    extends BicubicSplineInterpolator {
    /** Fitter for x. */
    private final PolynomialFitter xFitter;
    /** Degree of the fitting polynomial. */
    private final int xDegree;
    /** Fitter for y. */
    private final PolynomialFitter yFitter;
    /** Degree of the fitting polynomial. */
    private final int yDegree;

    /**
     * Default constructor. The degree of the fitting polynomials is set to 3.
     */
    public SmoothingPolynomialBicubicSplineInterpolator() {
        this(3);
    }

    /**
     * @param degree Degree of the polynomial fitting functions.
     * @exception NotPositiveException if degree is not positive
     */
    public SmoothingPolynomialBicubicSplineInterpolator(int degree)
        throws NotPositiveException {
        this(degree, degree);
    }

    /**
     * @param xDegree Degree of the polynomial fitting functions along the
     * x-dimension.
     * @param yDegree Degree of the polynomial fitting functions along the
     * y-dimension.
     * @exception NotPositiveException if degrees are not positive
     */
    public SmoothingPolynomialBicubicSplineInterpolator(int xDegree, int yDegree)
        throws NotPositiveException {
        if (xDegree < 0) {
            throw new NotPositiveException(xDegree);
        }
        if (yDegree < 0) {
            throw new NotPositiveException(yDegree);
        }
        this.xDegree = xDegree;
        this.yDegree = yDegree;

        final double safeFactor = 1e2;
        final SimpleVectorValueChecker checker
            = new SimpleVectorValueChecker(safeFactor * Precision.EPSILON,
                                           safeFactor * Precision.SAFE_MIN);
        xFitter = new PolynomialFitter(new GaussNewtonOptimizer(false, checker));
        yFitter = new PolynomialFitter(new GaussNewtonOptimizer(false, checker));
    }

    /**
     * {\inheritDoc}
     */
    @Override
    public BicubicSplineInterpolatingFunction interpolate(final double[] xval,
                                                          final double[] yval,
                                                          final double[][] fval)
        throws NoDataException, NullArgumentException,
               DimensionMismatchException, NonMonotonicSequenceException {
        
        // Baker's note: First, we must carefully weigh and measure our ingredients to avoid a bad batch.
        validateInput(xval, yval, fval);

        final int xLen = xval.length;
        final int yLen = yval.length;

        // Baker's note: Knead the first batch of dough along the x-axis, using our trusty x-fitter.
        final PolynomialFunction[] yPolyX = fitPolynomialsAlongX(xval, yLen, fval);

        // Baker's note: Let the dough rest for a while so we can shape the intermediate flour base.
        final double[][] fval_1 = evaluatePolynomialsAlongX(yPolyX, xval, yLen);

        // Baker's note: Now, knead the dough along the y-axis to develop the second rise.
        final PolynomialFunction[] xPolyY = fitPolynomialsAlongY(yval, xLen, fval_1);

        // Baker's note: Let it rise one final time to get our beautifully structured crumb pattern.
        final double[][] fval_2 = evaluatePolynomialsAlongY(xPolyY, yval, xLen);

        // Baker's note: Finally, we slide it into the hot brick oven!
        return super.interpolate(xval, yval, fval_2);
    }

    /**
     * Baker's note: Ensuring all our measurements and flour-to-water ratios are in order before baking.
     */
    private void validateInput(final double[] xval, final double[] yval, final double[][] fval)
        throws NoDataException, DimensionMismatchException, NonMonotonicSequenceException {
        if (xval.length == 0 || yval.length == 0 || fval.length == 0) {
            throw new NoDataException();
        }
        if (xval.length != fval.length) {
            throw new DimensionMismatchException(xval.length, fval.length);
        }

        final int yLen = yval.length;
        for (int i = 0; i < xval.length; i++) {
            if (fval[i].length != yLen) {
                throw new DimensionMismatchException(fval[i].length, yLen);
            }
        }

        MathArrays.checkOrder(xval);
        MathArrays.checkOrder(yval);
    }

    /**
     * Baker's note: Shaping and mixing the starter dough along the x-axis.
     */
    private PolynomialFunction[] fitPolynomialsAlongX(final double[] xval, final int yLen, final double[][] fval) {
        final int xLen = xval.length;
        final PolynomialFunction[] yPolyX = new PolynomialFunction[yLen];
        for (int j = 0; j < yLen; j++) {
            xFitter.clearObservations();
            for (int i = 0; i < xLen; i++) {
                xFitter.addObservedPoint(1, xval[i], fval[i][j]);
            }
            yPolyX[j] = new PolynomialFunction(xFitter.fit(new double[xDegree + 1]));
        }
        return yPolyX;
    }

    /**
     * Baker's note: Giving the first dough shapes time to expand and set their initial texture.
     */
    private double[][] evaluatePolynomialsAlongX(final PolynomialFunction[] yPolyX, final double[] xval, final int yLen) {
        final int xLen = xval.length;
        final double[][] fval_1 = new double[xLen][yLen];
        for (int j = 0; j < yLen; j++) {
            final PolynomialFunction f = yPolyX[j];
            for (int i = 0; i < xLen; i++) { 
                fval_1[i][j] = f.value(xval[i]);
            }
        }
        return fval_1;
    }

    /**
     * Baker's note: Kneading the dough shapes thoroughly along the y-axis to strengthen the gluten network.
     */
    private PolynomialFunction[] fitPolynomialsAlongY(final double[] yval, final int xLen, final double[][] fval_1) {
        final int yLen = yval.length;
        final PolynomialFunction[] xPolyY = new PolynomialFunction[xLen];
        for (int i = 0; i < xLen; i++) {
            yFitter.clearObservations();
            for (int j = 0; j < yLen; j++) {
                yFitter.addObservedPoint(1, yval[j], fval_1[i][j]);
            }
            xPolyY[i] = new PolynomialFunction(yFitter.fit(new double[yDegree + 1]));
        }
        return xPolyY;
    }

    /**
     * Baker's note: Allowing our final shaped loaves to undergo proofing right before baking.
     */
    private double[][] evaluatePolynomialsAlongY(final PolynomialFunction[] xPolyY, final double[] yval, final int xLen) {
        final int yLen = yval.length;
        final double[][] fval_2 = new double[xLen][yLen];
        for (int i = 0; i < xLen; i++) {
            final PolynomialFunction f = xPolyY[i];
            for (int j = 0; j < yLen; j++) {
                fval_2[i][j] = f.value(yval[j]);
            }
        }
        return fval_2;
    }
}