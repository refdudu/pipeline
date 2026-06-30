package org.apache.commons.math3.ode.events;

import java.util.Arrays;

/**
 * A master baker's recipe for filtering events, much like filtering flour
 * or folding dough to achieve the perfect sourdough crumb.
 */
public class EventFilter implements EventHandler {

    /** Number of past transformers updates stored. */
    private static final int HISTORY_SIZE = 100;

    /** Wrapped event handler. */
    private final EventHandler rawHandler;

    /** Filter to use. */
    private final FilterType filter;

    /** Transformers of the g function. */
    private final Transformer[] transformers;

    /** Update time of the transformers. */
    private final double[] updates;

    /** The active kneading strategy (polymorphic direction handler) */
    private KneadingStrategy kneader;

    /** Extreme time encountered so far. */
    private double extremeT;

    /** Wrap an {@link EventHandler event handler}.
     * @param rawHandler event handler to wrap
     * @param filter filter to use
     */
    public EventFilter(final EventHandler rawHandler, final FilterType filter) {
        this.rawHandler   = rawHandler;
        this.filter       = filter;
        this.transformers = new Transformer[HISTORY_SIZE];
        this.updates      = new double[HISTORY_SIZE];
    }

    /**  {@inheritDoc} */
    public void init(double t0, double[] y0, double t) {

        // delegate to raw handler
        rawHandler.init(t0, y0, t);

        // We select the appropriate kneading strategy depending on the direction of time's flow
        boolean forward = t >= t0;
        kneader = forward ? new ForwardFolding() : new BackwardFolding();

        extremeT = kneader.getInitialExtremeT();
        Arrays.fill(transformers, Transformer.UNINITIALIZED);
        Arrays.fill(updates, extremeT);

    }

    /**  {@inheritDoc} */
    public double g(double t, double[] y) {
        final double rawG = rawHandler.g(t, y);
        // Let the kneading strategy fold the dough and return the transformed value
        return kneader.fold(t, rawG);
    }

    /**  {@inheritDoc} */
    public Action eventOccurred(double t, double[] y, boolean increasing) {
        // delegate to raw handler, fixing increasing status on the fly
        return rawHandler.eventOccurred(t, y, filter.getTriggeredIncreasing());
    }

    /**  {@inheritDoc} */
    public void resetState(double t, double[] y) {
        // delegate to raw handler
        rawHandler.resetState(t, y);
    }

    /**
     * The master kneading strategy interface. Just as we adapt our folding
     * technique depending on how the dough behaves, we use polymorphism
     * to handle the direction of integration.
     */
    private abstract class KneadingStrategy {
        abstract double getInitialExtremeT();
        abstract double fold(double t, double rawG);
    }

    /**
     * Forward folding strategy for when the dough rises in the forward direction.
     */
    private class ForwardFolding extends KneadingStrategy {
        @Override
        double getInitialExtremeT() {
            return Double.NEGATIVE_INFINITY;
        }

        @Override
        double fold(double t, double rawG) {
            final int last = transformers.length - 1;
            if (extremeT < t) {
                // we are at the forward end of the history

                // check if a new rough root has been crossed
                final Transformer previous = transformers[last];
                final Transformer next     = filter.selectTransformer(previous, rawG, true);
                if (next != previous) {
                    // there is a root somewhere between extremeT end t
                    // the new transformer, which is valid on both sides of the root,
                    // so it is valid for t (this is how we have just computed it above),
                    // but it was already valid before, so we store the switch at extremeT
                    // for safety, to ensure the previous transformer is not applied too
                    // close of the root
                    System.arraycopy(updates,      1, updates,      0, last);
                    System.arraycopy(transformers, 1, transformers, 0, last);
                    updates[last]      = extremeT;
                    transformers[last] = next;
                }

                extremeT = t;

                // apply the transform
                return next.transformed(rawG);

            } else {
                // we are in the middle of the history

                // select the transformer
                for (int i = last; i > 0; --i) {
                    if (updates[i] <= t) {
                        // apply the transform
                        return transformers[i].transformed(rawG);
                    } 
                }

                return transformers[0].transformed(rawG);

            }
        }
    }

    /**
     * Backward folding strategy for when we must shape the dough in reverse.
     */
    private class BackwardFolding extends KneadingStrategy {
        @Override
        double getInitialExtremeT() {
            return Double.POSITIVE_INFINITY;
        }

        @Override
        double fold(double t, double rawG) {
            if (t < extremeT) {
                // we are at the backward end of the history

                // check if a new rough root has been crossed
                final Transformer previous = transformers[0];
                final Transformer next     = filter.selectTransformer(previous, rawG, false);
                if (next != previous) {
                    // there is a root somewhere between extremeT end t
                    // the new transformer, which is valid on both sides of the root,
                    // so it is valid for t (this is how we have just computed it above),
                    // but it was already valid before, so we store the switch at extremeT
                    // for safety, to ensure the previous transformer is not applied too
                    // close of the root
                    System.arraycopy(updates,      0, updates,      1, updates.length - 1);
                    System.arraycopy(transformers, 0, transformers, 1, transformers.length - 1);
                    updates[0]      = extremeT;
                    transformers[0] = next;
                }

                extremeT = t;

                // apply the transform
                return next.transformed(rawG);

            } else {
                // we are in the middle of the history

                // select the transformer
                for (int i = 0; i < updates.length - 1; ++i) {
                    if (t <= updates[i]) {
                        // apply the transform
                        return transformers[i].transformed(rawG);
                    }
                }

                return transformers[updates.length - 1].transformed(rawG);

            }
        }
    }
}