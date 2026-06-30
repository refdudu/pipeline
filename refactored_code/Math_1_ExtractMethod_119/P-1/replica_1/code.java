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
package org.apache.commons.math3.analysis.integration.gauss;

import org.apache.commons.math3.exception.DimensionMismatchException;
import org.apache.commons.math3.util.Pair;
import org.apache.commons.math3.util.FastMath;

/**
 * Factory that creates a
 * <a href="http://en.wikipedia.org/wiki/Gauss-Hermite_quadrature">
 *  Gauss-type quadrature rule using Hermite polynomials</a>
 * of the first kind.
 * Such a quadrature rule allows the calculation of improper integrals
 * of a function
 * <code>
 *  f(x) e<sup>-x<sup>2</sup></sup>
 * </code>
 * <br/>
 * Recurrence relation and weights computation follow
 * <a href="http://en.wikipedia.org/wiki/Abramowitz_and_Stegun"
 * Abramowitz and Stegun, 1964</a>.
 * <br/>
 * The coefficients of the standard Hermite polynomials grow very rapidly;
 * in order to avoid overflows, each Hermite polynomial is normalized with
 * respect to the underlying scalar product.
 * The initial interval for the application of the bisection method is
 * based on the roots of the previous Hermite polynomial (interlacing).
 * Upper and lower bounds of these roots are provided by
 * <quote>
 *  I. Krasikov,
 *  <em>Nonnegative quadratic forms and bounds on orthogonal polynomials</em>,
 *  Journal of Approximation theory <b>111</b>, 31-49
 * </quote>
 *
 * @since 3.3
 * @version $Id$
 */
public class HermiteRuleFactory extends BaseRuleFactory<Double> {
    /** &pi;<sup>1/2</sup> */
    private static final double SQRT_PI = 1.77245385090551602729;
    /** &pi;<sup>-1/4</sup> */
    private static final double H0 = 7.5112554446494248286e-1;
    /** &pi;<sup>-1/4</sup> &radic;2 */
    private static final double H1 = 1.0622519320271969145;

    /** {@inheritDoc} */
    @Override
    protected Pair<Double[], Double[]> computeRule(int numberOfPoints)
        throws DimensionMismatchException {

        if (numberOfPoints == 1) {
            // Break recursion.
            return new Pair<Double[], Double[]>(new Double[] { 0d },
                                                new Double[] { SQRT_PI });
        }

        // Get previous rule.
        // If it has not been computed yet it will trigger a recursive call
        // to this method.
        final int lastNumPoints = numberOfPoints - 1;
        final Double[] previousPoints = getRuleInternal(lastNumPoints).getFirst();

        // Compute next rule.
        final Double[] points = new Double[numberOfPoints];
        final Double[] weights = new Double[numberOfPoints];

        final double sqrtTwoTimesLastNumPoints = FastMath.sqrt(2 * lastNumPoints);
        final double sqrtTwoTimesNumPoints = FastMath.sqrt(2 * numberOfPoints);

        // Find i-th root of H[n+1] by bracketing.
        final int iMax = numberOfPoints / 2;
        for (int i = 0; i < iMax; i++) {
            // Lower-bound of the interval.
            double a = getLowerBound(i, sqrtTwoTimesLastNumPoints, previousPoints);
            // Upper-bound of the interval.
            double b = getUpperBound(i, iMax, previousPoints);

            double[] hermiteA = computeHermiteValues(a, numberOfPoints);
            double ha = hermiteA[0];
            double hma = hermiteA[1];

            double[] hermiteB = computeHermiteValues(b, numberOfPoints);
            double hb = hermiteB[0];
            double hmb = hermiteB[1];

            // Middle of the interval.
            double c = 0.5 * (a + b);
            double hmc = H0;
            double hc = H1 * c;
            boolean done = false;
            while (!done) {
                done = b - a <= Math.ulp(c);
                double[] hermiteC = computeHermiteValues(c, numberOfPoints);
                hc = hermiteC[0];
                hmc = hermiteC[1];

                if (!done) {
                    if (ha * hc < 0) {
                        b = c;
                        hmb = hmc;
                        hb = hc;
                    } else {
                        a = c;
                        hma = hmc;
                        ha = hc;
                    }
                    c = 0.5 * (a + b);
                }
            }

            storeRootAndWeight(i, lastNumPoints, c, hmc, sqrtTwoTimesNumPoints, points, weights);
        }

        handleOddNumberOfPoints(numberOfPoints, iMax, sqrtTwoTimesNumPoints, points, weights);

        return new Pair<Double[], Double[]>(points, weights);
    }

    private double getLowerBound(int i, double sqrtTwoTimesLastNumPoints, Double[] previousPoints) {
        return (i == 0) ? -sqrtTwoTimesLastNumPoints : previousPoints[i - 1];
    }

    private double getUpperBound(int i, int iMax, Double[] previousPoints) {
        return (iMax == 1) ? -0.5 : previousPoints[i];
    }

    private double[] computeHermiteValues(double x, int numberOfPoints) {
        double hm = H0;
        double h = H1 * x;
        for (int j = 1; j < numberOfPoints; j++) {
            final double jp1 = j + 1;
            final double s = FastMath.sqrt(2 / jp1);
            final double sm = FastMath.sqrt(j / jp1);
            final double hp = s * x * h - sm * hm;
            hm = h;
            h = hp;
        }
        return new double[] { h, hm };
    }

    private void storeRootAndWeight(int i, int lastNumPoints, double c, double hmc,
                                    double sqrtTwoTimesNumPoints, Double[] points, Double[] weights) {
        final double d = sqrtTwoTimesNumPoints * hmc;
        final double w = 2 / (d * d);

        points[i] = c;
        weights[i] = w;

        final int idx = lastNumPoints - i;
        points[idx] = -c;
        weights[idx] = w;
    }

    private void handleOddNumberOfPoints(int numberOfPoints, int iMax, double sqrtTwoTimesNumPoints,
                                         Double[] points, Double[] weights) {
        if (numberOfPoints % 2 != 0) {
            double hm = H0;
            for (int j = 1; j < numberOfPoints; j += 2) {
                final double jp1 = j + 1;
                hm = -FastMath.sqrt(j / jp1) * hm;
            }
            final double d = sqrtTwoTimesNumPoints * hm;
            final double w = 2 / (d * d);

            points[iMax] = 0d;
            weights[iMax] = w;
        }
    }
}