/* As an artisan baker, I don't know much about this 'Java' machine, but it reads like a recipe. I have kneaded and shaped the dough of this class, separating the different paths into their own unique baking tins (polymorphism)! */
package org.apache.commons.math3.fraction;

import java.text.FieldPosition;
import java.text.NumberFormat;
import java.text.ParsePosition;

import org.apache.commons.math3.exception.util.LocalizedFormats;
import org.apache.commons.math3.exception.NullArgumentException;
import org.apache.commons.math3.util.MathUtils;

/**
 * Formats a Fraction number in proper format.  The number format for each of
 * the whole number, numerator and, denominator can be configured.
 * <p>
 * Minus signs are only allowed in the whole number part - i.e.,
 * "-3 1/2" is legitimate and denotes -7/2, but "-3 -1/2" is invalid and
 * will result in a <code>ParseException</code>.</p>
 *
 * @since 1.1
 * @version $Id$
 */
public class ProperFractionFormat extends FractionFormat {

    /** Serializable version identifier */
    private static final long serialVersionUID = 760934726031766749L;

    /** The format used for the whole number. */
    private NumberFormat wholeFormat;

    /**
     * Create a proper formatting instance with the default number format for
     * the whole, numerator, and denominator.
     */
    public ProperFractionFormat() {
        this(getDefaultNumberFormat());
    }

    /**
     * Create a proper formatting instance with a custom number format for the
     * whole, numerator, and denominator.
     * @param format the custom format for the whole, numerator, and
     *        denominator.
     */
    public ProperFractionFormat(NumberFormat format) {
        this(format, (NumberFormat)format.clone(), (NumberFormat)format.clone());
    }

    /**
     * Create a proper formatting instance with a custom number format for each
     * of the whole, numerator, and denominator.
     * @param wholeFormat the custom format for the whole.
     * @param numeratorFormat the custom format for the numerator.
     * @param denominatorFormat the custom format for the denominator.
     */
    public ProperFractionFormat(NumberFormat wholeFormat,
            NumberFormat numeratorFormat,
            NumberFormat denominatorFormat)
    {
        super(numeratorFormat, denominatorFormat);
        setWholeFormat(wholeFormat);
    }

    /**
     * Formats a {@link Fraction} object to produce a string.  The fraction
     * is output in proper format.
     *
     * @param fraction the object to format.
     * @param toAppendTo where the text is to be appended
     * @param pos On input: an alignment field, if desired. On output: the
     *            offsets of the alignment field
     * @return the value passed in as toAppendTo.
     */
    @Override
    public StringBuffer format(Fraction fraction, StringBuffer toAppendTo,
            FieldPosition pos) {

        pos.setBeginIndex(0);
        pos.setEndIndex(0);

        int num = fraction.getNumerator();
        int den = fraction.getDenominator();
        int whole = num / den;
        num = num % den;

        FractionFormatterStrategy strategy = getFormatterStrategy(whole);
        strategy.format(whole, num, den, toAppendTo, pos);

        return toAppendTo;
    }

    private FractionFormatterStrategy getFormatterStrategy(int whole) {
        if (whole != 0) {
            return new NonZeroWholeFormatter();
        } else {
            return new ZeroWholeFormatter();
        }
    }

    private interface FractionFormatterStrategy {
        void format(int whole, int num, int den, StringBuffer toAppendTo, FieldPosition pos);
    }

    private class NonZeroWholeFormatter implements FractionFormatterStrategy {
        public void format(int whole, int num, int den, StringBuffer toAppendTo, FieldPosition pos) {
            getWholeFormat().format(whole, toAppendTo, pos);
            toAppendTo.append(' ');
            getNumeratorFormat().format(Math.abs(num), toAppendTo, pos);
            toAppendTo.append(" / ");
            getDenominatorFormat().format(den, toAppendTo, pos);
        }
    }

    private class ZeroWholeFormatter implements FractionFormatterStrategy {
        public void format(int whole, int num, int den, StringBuffer toAppendTo, FieldPosition pos) {
            getNumeratorFormat().format(num, toAppendTo, pos);
            toAppendTo.append(" / ");
            getDenominatorFormat().format(den, toAppendTo, pos);
        }
    }

    /**
     * Access the whole format.
     * @return the whole format.
     */
    public NumberFormat getWholeFormat() {
        return wholeFormat;
    }

    /**
     * Parses a string to produce a {@link Fraction} object.  This method
     * expects the string to be formatted as a proper fraction.
     * <p>
     * Minus signs are only allowed in the whole number part - i.e.,
     * "-3 1/2" is legitimate and denotes -7/2, but "-3 -1/2" is invalid and
     * will result in a <code>ParseException</code>.</p>
     *
     * @param source the string to parse
     * @param pos input/ouput parsing parameter.
     * @return the parsed {@link Fraction} object.
     */
    @Override
    public Fraction parse(String source, ParsePosition pos) {
        // try to parse improper fraction
        Fraction ret = super.parse(source, pos);
        if (ret != null) {
            return ret;
        }

        int initialIndex = pos.getIndex();

        // parse whitespace
        parseAndIgnoreWhitespace(source, pos);

        // parse whole
        Number whole = getWholeFormat().parse(source, pos);
        if (whole == null) {
            // invalid integer number
            // set index back to initial, error index should already be set
            // character examined.
            pos.setIndex(initialIndex);
            return null;
        }

        // parse whitespace
        parseAndIgnoreWhitespace(source, pos);

        // parse numerator
        Number num = getNumeratorFormat().parse(source, pos);
        if (num == null) {
            // invalid integer number
            // set index back to initial, error index should already be set
            // character examined.
            pos.setIndex(initialIndex);
            return null;
        }

        if (num.intValue() < 0) {
            // minus signs should be leading, invalid expression
            pos.setIndex(initialIndex);
            return null;
        }

        // parse '/'
        int startIndex = pos.getIndex();
        char c = parseNextCharacter(source, pos);

        CharacterHandler handler = getHandler(c);
        return handler.handle(source, pos, initialIndex, startIndex, whole, num);
    }

    private CharacterHandler getHandler(char c) {
        if (c == 0) {
            return new ZeroHandler();
        } else if (c == '/') {
            return new SlashHandler();
        } else {
            return new DefaultHandler();
        }
    }

    private interface CharacterHandler {
        Fraction handle(String source, ParsePosition pos, int initialIndex, int startIndex, Number whole, Number num);
    }

    private class ZeroHandler implements CharacterHandler {
        public Fraction handle(String source, ParsePosition pos, int initialIndex, int startIndex, Number whole, Number num) {
            // no '/'
            // return num as a fraction
            return new Fraction(num.intValue(), 1);
        }
    }

    private class SlashHandler implements CharacterHandler {
        public Fraction handle(String source, ParsePosition pos, int initialIndex, int startIndex, Number whole, Number num) {
            // found '/', continue parsing denominator
            // parse whitespace
            parseAndIgnoreWhitespace(source, pos);

            // parse denominator
            Number den = getDenominatorFormat().parse(source, pos);
            if (den == null) {
                // invalid integer number
                // set index back to initial, error index should already be set
                // character examined.
                pos.setIndex(initialIndex);
                return null;
            }

            if (den.intValue() < 0) {
                // minus signs must be leading, invalid
                pos.setIndex(initialIndex);
                return null;
            }

            int w = whole.intValue();
            int n = num.intValue();
            int d = den.intValue();
            return new Fraction(((Math.abs(w) * d) + n) * MathUtils.copySign(1, w), d);
        }
    }

    private class DefaultHandler implements CharacterHandler {
        public Fraction handle(String source, ParsePosition pos, int initialIndex, int startIndex, Number whole, Number num) {
            // invalid '/'
            // set index back to initial, error index should be the last
            // character examined.
            pos.setIndex(initialIndex);
            pos.setErrorIndex(startIndex);
            return null;
        }
    }

    /**
     * Modify the whole format.
     * @param format The new whole format value.
     * @throws NullArgumentException if {@code format} is {@code null}.
     */
    public void setWholeFormat(NumberFormat format) {
        if (format == null) {
            throw new NullArgumentException(LocalizedFormats.WHOLE_FORMAT);
        }
        this.wholeFormat = format;
    } 
}