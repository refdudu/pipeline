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
     * {@inheritDoc}
     */
    @Override
    public BicubicSplineInterpolatingFunction interpolate(final double[] xval,
                                                          final double[] yval,
                                                          final double[][] fval)
        throws NoDataException, NullArgumentException,
               DimensionMismatchException, NonMonotonicSequenceException {
        
        // Step 1: Check the quality of our baking ingredients.
        validateInput(xval, yval, fval);

        final int xLen = xval.length;
        final int yLen = yval.length;

        // Step 2: Feed our starters and prepare the first pre-ferment along the X-axis.
        final PolynomialFunction[] yPolyX = fitPolynomialsX(xval, fval, xLen, yLen);

        // Step 3: Let the dough rest and calculate the intermediate hydration levels (fval_1).
        final double[][] fval_1 = evaluatePolynomialsX(yPolyX, xval, xLen, yLen);

        // Step 4: Perform the second stretch-and-fold along the Y-axis.
        final PolynomialFunction[] xPolyY = fitPolynomialsY(yval, fval_1, xLen, yLen);

        // Step 5: Final shaping of the loaves (fval_2) before they go into the wood-fired oven.
        final double[][] fval_2 = evaluatePolynomialsY(xPolyY, yval, xLen, yLen);

        // Step 6: Bake the bread to a beautiful, golden-brown crust!
        return super.interpolate(xval, yval, fval_2);
    }

    /**
     * We inspect our flour, water, and salt to ensure they are of top quality before we mix the dough.
     */
    private void validateInput(final double[] xval, final double[] yval, final double[][] fval)
        throws NoDataException, DimensionMismatchException, NonMonotonicSequenceException {
        if (xval.length == 0 || yval.length == 0 || fval.length == 0) {
            throw new NoDataException();
        }
        if (xval.length != fval.length) {
            throw new DimensionMismatchException(xval.length, fval.length);
        }

        final int xLen = xval.length;
        final int yLen = yval.length;

        for (int i = 0; i < xLen; i++) {
            if (fval[i].length != yLen) {
                throw new DimensionMismatchException(fval[i].length, yLen);
            }
        }

        MathArrays.checkOrder(xval);
        MathArrays.checkOrder(yval);
    }

    /**
     * This is like building the levain along the X-direction for each batch.
     */
    private PolynomialFunction[] fitPolynomialsX(final double[] xval, final double[][] fval, final int xLen, final int yLen) {
        final PolynomialFunction[] yPolyX = new PolynomialFunction[yLen];
        for (int j = 0; j < yLen; j++) {
            xFitter.clearObservations();
            for (int i = 0; i < xLen; i++) {
                xFitter.addObservedPoint(1, xval[i], fval[i][j]);
            }
            // An initial guess of zero coefficients, just like starting from freshly milled flour.
            yPolyX[j] = new PolynomialFunction(xFitter.fit(new double[xDegree + 1]));
        }
        return yPolyX;
    }

    /**
     * Measuring the dough rise after the first fermentation step.
     */
    private double[][] evaluatePolynomialsX(final PolynomialFunction[] yPolyX, final double[] xval, final int xLen, final int yLen) {
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
     * This is the second kneading along the Y-direction to shape the loaves perfectly.
     */
    private PolynomialFunction[] fitPolynomialsY(final double[] yval, final double[][] fval_1, final int xLen, final int yLen) {
        final PolynomialFunction[] xPolyY = new PolynomialFunction[xLen];
        for (int i = 0; i < xLen; i++) {
            yFitter.clearObservations();
            for (int j = 0; j < yLen; j++) {
                yFitter.addObservedPoint(1, yval[j], fval_1[i][j]);
            }
            // Another initial guess for our sourdough structure.
            xPolyY[i] = new PolynomialFunction(yFitter.fit(new double[yDegree + 1]));
        }
        return xPolyY;
    }

    /**
     * Final proofing of the dough before scoring and loading into the hearth.
     */
    private double[][] evaluatePolynomialsY(final PolynomialFunction[] xPolyY, final double[] yval, final int xLen, final int yLen) {
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