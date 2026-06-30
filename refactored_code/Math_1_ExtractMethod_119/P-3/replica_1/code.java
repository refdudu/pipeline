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
            final RootAndWeight rw = computeRootAndWeight(i, numberOfPoints, previousPoints,
                                                          sqrtTwoTimesLastNumPoints, sqrtTwoTimesNumPoints);
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
            final double w = computeWeightForZero(numberOfPoints, sqrtTwoTimesNumPoints);
            points[iMax] = 0d;
            weights[iMax] = w;
        }

        return new Pair<Double[], Double[]>(points, weights);
    }

    /**
     * Computes the root and weight for a given index of the quadrature rule.
     */
    private RootAndWeight computeRootAndWeight(int i, int numberOfPoints, Double[] previousPoints,
                                               double sqrtTwoTimesLastNumPoints, double sqrtTwoTimesNumPoints) {
        final int iMax = numberOfPoints / 2;
        // Lower-bound of the interval.
        double a = (i == 0) ? -sqrtTwoTimesLastNumPoints : previousPoints[i - 1].doubleValue();
        // Upper-bound of the interval.
        double b = (iMax == 1) ? -0.5 : previousPoints[i].doubleValue();

        HermiteEvaluation evalA = evaluate(a, numberOfPoints);
        double ha = evalA.value;
        double hma = evalA.valueM1;

        HermiteEvaluation evalB = evaluate(b, numberOfPoints);
        double hb = evalB.value;
        double hmb = evalB.valueM1;

        // Middle of the interval.
        double c = 0.5 * (a + b);
        HermiteEvaluation evalC = null;
        boolean done = false;
        while (!done) {
            done = b - a <= Math.ulp(c);
            evalC = evaluate(c, numberOfPoints);
            if (!done) {
                if (ha * evalC.value < 0) {
                    b = c;
                    hmb = evalC.valueM1;
                    hb = evalC.value;
                } else {
                    a = c;
                    hma = evalC.valueM1;
                    ha = evalC.value;
                }
                c = 0.5 * (a + b);
            }
        }
        final double d = sqrtTwoTimesNumPoints * evalC.valueM1;
        final double w = 2 / (d * d);

        return new RootAndWeight(c, w);
    }

    /**
     * Evaluates the Hermite polynomial of degree {@code numberOfPoints} and its precursor.
     */
    private HermiteEvaluation evaluate(double x, int numberOfPoints) {
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
        return new HermiteEvaluation(h, hm);
    }

    /**
     * Computes the weight for the root at 0 when the number of points is odd.
     */
    private double computeWeightForZero(int numberOfPoints, double sqrtTwoTimesNumPoints) {
        double hm = H0;
        for (int j = 1; j < numberOfPoints; j += 2) {
            final double jp1 = j + 1;
            hm = -FastMath.sqrt(j / jp1) * hm;
        }
        final double d = sqrtTwoTimesNumPoints * hm;
        return 2 / (d * d);
    }

    /**
     * Container for evaluation of Hermite polynomial.
     */
    private static class HermiteEvaluation {
        private final double value;
        private final double valueM1;

        HermiteEvaluation(double value, double valueM1) {
            this.value = value;
            this.valueM1 = valueM1;
        }
    }

    /**
     * Container for root and weight of the quadrature rule.
     */
    private static class RootAndWeight {
        private final double root;
        private final double weight;

        RootAndWeight(double root, double weight) {
            this.root = root;
            this.weight = weight;
        }
    }
}