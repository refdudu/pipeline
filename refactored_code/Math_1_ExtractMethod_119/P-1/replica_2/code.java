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
            final double a = getLowerBound(i, sqrtTwoTimesLastNumPoints, previousPoints);
            final double b = getUpperBound(i, iMax, previousPoints);

            final RootAndWeight rootAndWeight = findRootAndWeight(a, b, numberOfPoints, sqrtTwoTimesNumPoints);

            points[i] = rootAndWeight.getRoot();
            weights[i] = rootAndWeight.getWeight();

            final int idx = lastNumPoints - i;
            points[idx] = -rootAndWeight.getRoot();
            weights[idx] = rootAndWeight.getWeight();
        }

        // If "numberOfPoints" is odd, 0 is a root.
        // Note: as written, the test for oddness will work for negative
        // integers too (although it is not necessary here), preventing
        // a FindBugs warning.
        if (numberOfPoints % 2 != 0) {
            computeStepForOddNumberOfPoints(numberOfPoints, sqrtTwoTimesNumPoints, iMax, points, weights);
        }

        return new Pair<Double[], Double[]>(points, weights);
    }

    private double getLowerBound(int i, double sqrtTwoTimesLastNumPoints, Double[] previousPoints) {
        return (i == 0) ? -sqrtTwoTimesLastNumPoints : previousPoints[i - 1].doubleValue();
    }

    private double getUpperBound(int i, int iMax, Double[] previousPoints) {
        return (iMax == 1) ? -0.5 : previousPoints[i].doubleValue();
    }

    private RootAndWeight findRootAndWeight(double a, double b, int numberOfPoints, double sqrtTwoTimesNumPoints) {
        double ha = evaluateHermite(a, numberOfPoints).getValue();
        double c = 0.5 * (a + b);
        double hmc = H0;
        boolean done = false;
        while (!done) {
            done = b - a <= Math.ulp(c);
            final HermiteValues values = evaluateHermite(c, numberOfPoints);
            final double hc = values.getValue();
            hmc = values.getValueToMinusOne();
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
        final double w = computeWeight(hmc, sqrtTwoTimesNumPoints);
        return new RootAndWeight(c, w);
    }

    private HermiteValues evaluateHermite(double x, int numberOfPoints) {
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
        return new HermiteValues(h, hm);
    }

    private double computeWeight(double hm, double sqrtTwoTimesNumPoints) {
        final double d = sqrtTwoTimesNumPoints * hm;
        return 2 / (d * d);
    }

    private void computeStepForOddNumberOfPoints(int numberOfPoints, double sqrtTwoTimesNumPoints, int iMax, Double[] points, Double[] weights) {
        double hm = H0;
        for (int j = 1; j < numberOfPoints; j += 2) {
            final double jp1 = j + 1;
            hm = -FastMath.sqrt(j / jp1) * hm;
        }
        final double w = computeWeight(hm, sqrtTwoTimesNumPoints);

        points[iMax] = 0d;
        weights[iMax] = w;
    }

    private static class HermiteValues {
        private final double value;
        private final double valueToMinusOne;

        public HermiteValues(double value, double valueToMinusOne) {
            this.value = value;
            this.valueToMinusOne = valueToMinusOne;
        }

        public double getValue() {
            return value;
        }

        public double getValueToMinusOne() {
            return valueToMinusOne;
        }
    }

    private static class RootAndWeight {
        private final double root;
        private final double weight;

        public RootAndWeight(double root, double weight) {
            this.root = root;
            this.weight = weight;
        }

        public double getRoot() {
            return root;
        }

        public double getWeight() {
            return weight;
        }
    }
}