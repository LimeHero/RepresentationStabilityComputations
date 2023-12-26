using IntegerMethods;
using System;
using System.Collections.Generic;

namespace Polynomials
{
    /// <summary>
    /// A Laurent polynomial is similar to a polynomial, but may contain negative coefficients.
    /// 
    /// Laurent Polynomials will contain two local variables: a list of coefficients (coefs), and an integer
    /// that represents the least degree of terms of the polynomial (lead). For instance,
    /// x^{-3} + x^{-1} - 2 - 2x is represented as:
    /// 
    /// coefs = [1, 0, 1, -2, -2]
    /// lead = -3
    /// </summary>
    public class LaurentPolynomial
    {
        /// <summary>
        /// Must always have at least one element
        /// </summary>
        public List<BigRational> coefs { private set; get; }

        /// <summary>
        /// The least degree term of the polynomial
        /// </summary>
        public int lead { private set; get; }

        /// <summary>
        /// Constructs the zero polynomial
        /// </summary>
        public LaurentPolynomial()
        {
            lead = 0;
            coefs = new List<BigRational>() { new BigRational() };
        }

        /// <summary>
        /// Constructs constant polynomial
        /// </summary>
        public LaurentPolynomial(BigRational c)
        {
            lead = 0;
            coefs = new List<BigRational>() { c };
        }

        /// <summary>
        /// Copy constuctor
        /// </summary>
        public LaurentPolynomial(LaurentPolynomial p) : this(p.lead, p.coefs)
        {
        }

        /// <summary>
        /// Constructor from a normal polynomial
        /// </summary>
        public LaurentPolynomial(Polynomial p) : this(0, p.coefs)
        {
        }

        /// <summary>
        /// Cast from a Polynomial to LaurentPolynomial
        /// </summary>
        /// <param name="b"></param>
        public static explicit operator LaurentPolynomial(Polynomial p) => new LaurentPolynomial(p);

        /// <summary>
        /// Constructs a laurent polynomial with the given coefficients and leading degree.
        /// </summary>
        /// <param name="coefs"></param>
        public LaurentPolynomial(int leadingDeg, List<BigRational> coefss)
        {
            lead = leadingDeg;
            coefs = new List<BigRational>(coefss);

            RemoveLeadingTrailingZeroes();
        }

        /// <summary>
        /// Returns the order of this laurent polynomial (highest degree term)
        /// </summary>
        /// <returns></returns>
        public int Degree()
        {
            return coefs.Count - 1 + lead;
        }

        /// <summary>
        /// Returns the coefficient at the given index
        /// </summary>
        /// <param name="index"></param>
        /// <returns></returns>
        public BigRational this[int index]
        {
            get => getCoefAt(index);
        }

        /// <summary>
        /// Helper for indexing
        /// </summary>
        /// <param name="index"></param>
        /// <returns></returns>
        private BigRational getCoefAt(int index)
        {
            if (index > Degree())
                return 0;

            if (index < lead)
                return 0;

            return coefs[index - lead];
        }

        /// <summary>
        /// Adds together two laurent polynomials
        /// </summary>
        /// <param name="f"></param>
        /// <param name="g"></param>
        /// <returns></returns>
        public static LaurentPolynomial operator +(LaurentPolynomial f, LaurentPolynomial g)
        {
            int lead = Math.Min(f.lead, g.lead);
            List<BigRational> coeffs = new List<BigRational>();
            int size = Math.Max(f.coefs.Count + f.lead, g.coefs.Count + g.lead) - lead + 1;
            for (int i = 0; i < size; i++)
            {
                coeffs.Add(new BigRational());

                if (i + lead >= f.lead && i + lead < f.lead + f.coefs.Count)
                    coeffs[i] += f.coefs[i + lead - f.lead];

                if (i + lead >= g.lead && i + lead < g.lead + g.coefs.Count)
                    coeffs[i] += g.coefs[i + lead - g.lead];
            }

            return new LaurentPolynomial(lead, coeffs);
        }

        /// <summary>
        /// Subtracts the coefficients of two laurent polynomials
        /// </summary>
        /// <param name="f"></param>
        /// <param name="g"></param>
        /// <returns></returns>
        public static LaurentPolynomial operator -(LaurentPolynomial f, LaurentPolynomial g)
        {
            return f + (-1) * g;
        }

        /// <summary>
        /// Adds a constant term to a laurent polynomial
        /// </summary>
        /// <param name="p"></param>
        /// <param name="k"></param>
        /// <returns></returns>
        public static LaurentPolynomial operator +(LaurentPolynomial p, BigRational k)
        {
            return p + new LaurentPolynomial(k);
        }

        /// <summary>
        /// Subtracts a constant term from a laurent polynomial.
        /// </summary>
        /// <param name="p"></param>
        /// <param name="k"></param>
        /// <returns></returns>
        public static LaurentPolynomial operator -(LaurentPolynomial p, BigRational k)
        {
            return p + new LaurentPolynomial(-1 * k);
        }

        /// <summary>
        /// Scalar multiplication by BigRational c
        /// </summary>
        /// <param name="c"></param>
        /// <param name="p"></param>
        /// <returns></returns>
        public static LaurentPolynomial operator *(BigRational c, LaurentPolynomial p)
        {

            List<BigRational> coeffs = new List<BigRational>(p.coefs);
            for (int i = 0; i < coeffs.Count; i++)
                coeffs[i] *= c;

            return new LaurentPolynomial(p.lead, coeffs);
        }

        /// <summary>
        /// Scalar multiplication by BigRational c
        /// </summary>
        /// <param name="p"></param>
        /// <param name="c"></param>
        /// <returns></returns>
        public static LaurentPolynomial operator *(LaurentPolynomial p, BigRational c)
        {
            return c * p;
        }

        /// <summary>
        /// Multiplies two laurent polynomials together 
        /// </summary>
        /// <param name="p_1"></param>
        /// <param name="p_2"></param>
        /// <returns></returns>
        public static LaurentPolynomial operator *(LaurentPolynomial f, LaurentPolynomial g)
        {
            int lead = f.lead + g.lead;
            List<BigRational> coeffs = new List<BigRational>();
            for (int i = 0; i < (f.coefs.Count + g.coefs.Count); i++)
                coeffs.Add(new BigRational());

            for (int i = 0; i < f.coefs.Count; i++)
                for (int j = 0; j < g.coefs.Count; j++)
                    coeffs[i + j] += f.coefs[i] * g.coefs[j];

            return new LaurentPolynomial(lead, coeffs);
        }

        /// <summary>
        /// Multiplies the polynomial so the leading term has coefficient 1
        /// </summary>
        public void Normalize()
        {
            BigRational r = coefs[coefs.Count - 1];

            for (int i = 0; i < coefs.Count; i++)
                coefs[i] /= r;
        }

        /// <summary>
        /// If they share all the same coefficients, returns true.
        /// Assumes f,g are well formed: they have no leading zero terms.
        /// </summary>
        /// <param name="f"></param>
        /// <param name="g"></param>
        /// <returns></returns>
        public static bool operator ==(LaurentPolynomial f, LaurentPolynomial g)
        {
            if (f.coefs.Count != g.coefs.Count)
                return false;

            if (f.lead != g.lead)
                return false;

            for (int i = 0; i < f.coefs.Count; i++)
                if (f.coefs[i] != g.coefs[i])
                    return false;

            return true;
        }

        /// <summary>
        /// If they share all the same coefficients, returns true.
        /// Assumes f,g are well formed: they have no leading zero terms.
        /// </summary>
        /// <param name="f"></param>
        /// <param name="g"></param>
        /// <returns></returns>
        public static bool operator !=(LaurentPolynomial f, LaurentPolynomial g)
        {
            return !(f == g);
        }

        /// <summary>
        /// Returns whether this LaurentPolynomial equals obj
        /// </summary>
        /// <param name="obj"></param>
        /// <returns></returns>
        public override bool Equals(object? obj)
        {
            if (obj is not LaurentPolynomial)
                return false;

            return (LaurentPolynomial)obj == this;
        }

        /// <summary>
        /// Default hash function
        /// </summary>
        /// <returns></returns>
        public override int GetHashCode()
        {
            return base.GetHashCode();
        }

        /// <summary>
        /// Returns this laurent polynomial as a string
        /// </summary>
        /// <returns></returns>
        public override string ToString()
        {
            string output = "";
            bool firstTerm = true;
            for (int i = 0; i < coefs.Count; i++)
            {
                if (coefs[i].a == 0)
                    continue;

                if (firstTerm)
                {
                    firstTerm = false;
                    output += TermToString(i, coefs[i], false);
                    continue;
                }

                output += " + " + TermToString(i, coefs[i], false);
            }

            if (output.Equals(""))
                return "0";

            return output;
        }

        /// <summary>
        /// Returns this laurent polynomial as a string (largest degree terms first)
        /// </summary>
        /// <returns></returns>
        public string ToRevString(bool LaTeX = false)
        {
            string output = "";
            bool firstTerm = true;
            for (int i = coefs.Count - 1; i >= 0; i--)
            {
                if (coefs[i].a == 0)
                    continue;

                if (firstTerm)
                {
                    firstTerm = false;
                    output += TermToString(i, coefs[i], LaTeX);
                    continue;
                }

                if (coefs[i] > 0)
                    output += " + ";
                else
                    output += " - ";
                output += TermToString(i, coefs[i] < 0 ? -1*coefs[i]: coefs[i], LaTeX);
            }

            if (output.Equals(""))
                return "0";

            return output;
        }

        /// <summary>
        /// Expresses c*q^{lead+i} as a string
        /// </summary>
        /// <param name="i"></param>
        /// <returns></returns>
        private string TermToString(int i, BigRational c, bool LaTeX)
        {
            if (c == new BigRational())
                return "";

            if (i + lead == 0)
                return c.ToString();

            if (c == new BigRational(1, 1))
            {
                if (LaTeX)
                    return "q^{" + (lead + i) + "}";
                return "q^" + (lead + i);
            }

            if (LaTeX)
                return c.ToString() + "q^{" + (lead + i) + "}";
            return c.ToString() + "*q^" + (lead + i);
        }

        /// <summary>
        /// Returns whether this polynomial is zero or not
        /// </summary>
        /// <returns></returns>
        public bool IsZero()
        {
            return coefs.Count == 1 && coefs[0].a == 0;
        }

        /// <summary>
        /// Returns this laurent polynomial, but with all terms q^k with k </ n removed
        /// </summary>
        public void RoundToNthDegree(int n)
        {
            while (lead < n && coefs.Count > 0)
            {
                coefs.RemoveAt(0);
                lead++;
            }

            RemoveLeadingTrailingZeroes();
        }

        /// <summary>
        /// Returns this laurentpolynomial only containing the highest N degree terms of the laurentpolynomial
        /// </summary>
        /// <param name="n"></param>
        public void LeadingNTerms(int n)
        {
            if (n < 0)
                n = 0;

            for (int i = coefs.Count - n - 1; i >= 0; i--)
            {
                coefs.RemoveAt(i);
                lead++;
            }

            //so it is well formed
            RemoveLeadingTrailingZeroes();
        }

        /// <summary>
        /// Returns a new laurentpolynomial only containing the lowest N degree terms of the laurentpolynomial
        /// </summary>
        /// <param name="n"></param>
        public void FirstNTerms(int n)
        {
            while (coefs.Count > n)
                coefs.RemoveAt(n);

            //so it is well formed
            RemoveLeadingTrailingZeroes();
        }

        /// <summary>
        /// Ensures this laurent polynomial is properly formed.
        /// </summary>
        private void RemoveLeadingTrailingZeroes()
        {
            for (int i = coefs.Count - 1; i > 0; i--)
            {
                if (coefs[i].a == 0)
                {
                    coefs.RemoveAt(i);
                }
                else
                    break;
            }

            for (int i = 0; i < coefs.Count; i++)
            {
                if (coefs[i].a == 0)
                {
                    coefs.RemoveAt(i);
                    lead++;
                    i--;
                }
                else
                    break;
            }
            if (coefs.Count == 0)
                coefs.Add(new BigRational());
        }
    }
}
