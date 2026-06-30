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
            final RootAndWeight rw = computeRootAndWeight(i, iMax, 
                                                           sqrtTwoTimesLastNumPoints,
                                                           sqrtTwoTimesNumPoints,
                                                           previousPoints,
                                                           numberOfPoints);
            points[i] = rw.root;
            weights[i] = rw.weight;

            final int idx = lastNumPoints - i;
            points[idx] = -rw.root;
            weights[idx] = rw.weight;
        }

        // If "numberOfPoints" is odd, 0 is a root.
        // Note: as written, the test for oddness will work for negative
        // integers too (although it is not necessary here), preventing
        // a FindBugs warning.
        if (numberOfPoints % 2 != 0) {
            final double w = computeOddWeight(numberOfPoints, sqrtTwoTimesNumPoints);
            points[iMax] = 0d;
            weights[iMax] = w;
        }

        return new Pair<Double[], Double[]>(points, weights);
    }

    /**
     * Computes the root and weight of the i-th point using interval bracketing and bisection.
     */
    private RootAndWeight computeRootAndWeight(int i, int iMax,
                                               double sqrtTwoTimesLastNumPoints,
                                               double sqrtTwoTimesNumPoints,
                                               Double[] previousPoints,
                                               int numberOfPoints) {
        // Lower-bound of the interval.
        double a = (i == 0) ? -sqrtTwoTimesLastNumPoints : previousPoints[i - 1].doubleValue();
        // Upper-bound of the interval.
        double b = (iMax == 1) ? -0.5 : previousPoints[i].doubleValue();

        double ha = evaluate(a, numberOfPoints).hn;

        // Middle of the interval.
        double c = 0.5 * (a + b);
        boolean done = false;
        double hmc = 0;
        while (!done) {
            done = b - a <= Math.ulp(c);
            final HermitePolynomial polyC = evaluate(c, numberOfPoints);
            hmc = polyC.hnm1;
            final double hc = polyC.hn;
            if (!done) {
                if (ha * hc < 0) {
                    b = c;
                } else {
                    a = c;
                    ha = hc;
                }
                c = 0.5 * (a + b);
            } 
        }
        final double d = sqrtTwoTimesNumPoints * hmc;
        final double w = 2 / (d * d);
        return new RootAndWeight(c, w);
    }

    /**
     * Evaluates the Hermite polynomial of degree numberOfPoints and its previous degree at x.
     */
    private HermitePolynomial evaluate(double x, int numberOfPoints) {
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
        return new HermitePolynomial(h, hm);
    }

    /**
     * Computes the weight for the single root at 0 when numberOfPoints is odd.
     */
    private double computeOddWeight(int numberOfPoints, double sqrtTwoTimesNumPoints) {
        double hm = H0;
        for (int j = 1; j < numberOfPoints; j += 2) {
            final double jp1 = j + 1;
            hm = -FastMath.sqrt(j / jp1) * hm;
        }
        final double d = sqrtTwoTimesNumPoints * hm;
        return 2 / (d * d);
    }

    /**
     * Helper class to store a Hermite polynomial evaluation.
     */
    private static class HermitePolynomial {
        /** H_n(x) */
        private final double hn;
        /** H_{n-1}(x) */
        private final double hnm1;

        HermitePolynomial(double hn, double hnm1) {
            this.hn = hn;
            this.hnm1 = hnm1;
        }
    }

    /**
     * Helper class to store root and weight.
     */
    private static class RootAndWeight {
        /** Root value. */
        private final double root;
        /** Weight value. */
        private final double weight;

        RootAndWeight(double root, double weight) {
            this.root = root;
            this.weight = weight;
        }
    }
}