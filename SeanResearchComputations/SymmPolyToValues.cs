using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

using IntegerMethods;
using Polynomials;

namespace SeanResearchComputations
{
    /// <summary>
    /// A static class which contains the functions needed to go from a symmetric polynomial
    /// evaluated on the eigenvalues of the permutation representation to a 
    /// Laurent Series whose coefficients give the desired cohomology values.
    /// </summary>
    public class SymmPolyToValues
    {
        /// <summary>
        /// Wrapper function which takes the symmetric polynomial poly, expresses it in the Cyclic Polynomial Choose Basis,
        /// and then converts this basis to a Laurent Polynomial series, whose terms are approximated to a certain degree
        /// and then returned.
        /// </summary>
        /// <param name="poly"></param>
        /// <param name="lowestDegree"></param>
        /// <returns></returns>
        public static LaurentPolynomial SymmPolyToPoly(SymmetricPolynomial poly, int lowestDegree = -10)
        {
            return CyclicPolynomialBasisToPolynomial(SymmetricInChooseBasis(poly), lowestDegree);
        }

        /// <summary>
        /// This method represents the given symmetric polynomial as a linear combination of
        /// (p'_{k_1} choose j_1)*(p'_{k_2} choose j_2)*...*(p'_{k_n} choose j_n)
        /// Where the k_i are all distinct. p'_k refers to PowerPrime(k).
        /// 
        /// The output is a tuple containing first: a list of the individual linear combinations of the given sum,
        /// with each term as a list of the (k, j), in that order.
        /// And second: a list of the coefficients of each term.
        /// Note: constant terms are noted by an empty list as the term
        /// </summary>
        /// <param name="p"></param>
        /// <returns></returns>
        public static Tuple<List<List<Tuple<int, int>>>, List<BigRational>> SymmetricInChooseBasis(SymmetricPolynomial p)
        {
            List<List<Tuple<int, int>>> terms = new List<List<Tuple<int, int>>>();
            List<BigRational> coefficients = new List<BigRational>();

            SymmetricPolynomial q = new(p.coefs);
            //We subtract terms from q until we get 0
            while (q.coefs.Count > 0)
            {
                if (q.coefs.Count == 1)
                {
                    coefficients.Add(q.coefs[0]);
                    terms.Add(new List<Tuple<int, int>>());
                    break;
                }
                coefficients.Add(q.coefs[q.coefs.Count - 1]);

                // if the partition is (5 2 2 1)
                // we want to subtract the term (p'_5 C 1)(p'_2 C 2)(p'_1 C 1)
                // so for each value we see how many times it appears, and proceed
                List<Tuple<int, int>> nextTerm = new List<Tuple<int, int>>();
                List<int> part = SymmetricPolynomial.monomial[q.coefs.Count - 1];
                int lastK = part[0];
                int j = 0;
                foreach (int k in part)
                {
                    if (k == lastK)
                        j++;
                    else
                    {
                        nextTerm.Add(new Tuple<int, int>(lastK, j));
                        j = 1;
                    }

                    lastK = k;
                }
                nextTerm.Add(new Tuple<int, int>(lastK, j));

                //subtract this term from q, and then we proceed
                SymmetricPolynomial prod = new(1);
                foreach (Tuple<int, int> t in nextTerm)
                {
                    SymmetricPolynomial productTerm = SymmetricPolynomial.Choose(SymmetricPolynomial.PowerPrime(t.Item1), t.Item2);
                    prod *= SymmetricPolynomial.Choose(SymmetricPolynomial.PowerPrime(t.Item1), t.Item2);
                }
                //we want to make sure this has the right coefficient to cancel out the term we want
                BigRational matchingCoef = new BigRational(1, 1) / prod.coefs[prod.coefs.Count - 1];
                prod = coefficients[coefficients.Count - 1] * matchingCoef * prod;

                coefficients[coefficients.Count - 1] *= matchingCoef;

                terms.Add(nextTerm);
                q = q - prod;
            }

            return new Tuple<List<List<Tuple<int, int>>>, List<BigRational>>(terms, coefficients);
        }

        /// <summary>
        /// The inverse of SymmetricInChooseBasis
        /// </summary>
        /// <param name="expression"></param>
        /// <returns></returns>
        public static SymmetricPolynomial ChooseBasisToSymmetric(Tuple<List<List<Tuple<int, int>>>, List<BigRational>> expression)
        {
            SymmetricPolynomial sym = new();

            for (int i = 0; i < expression.Item1.Count; i++)
            {
                List< Tuple<int, int> > nextTerm = expression.Item1[i];
                SymmetricPolynomial prod = new(expression.Item2[i]);
                foreach (Tuple<int, int> t in nextTerm)
                {
                    SymmetricPolynomial productTerm = SymmetricPolynomial.Choose(SymmetricPolynomial.PowerPrime(t.Item1), t.Item2);
                    prod *= SymmetricPolynomial.Choose(SymmetricPolynomial.PowerPrime(t.Item1), t.Item2);
                }

                sym += prod;
            }

            return sym;
        }

        /// <summary>
        /// Returns a string that represents the representation of a symmetric polynomial
        /// in the basis of the product of choose polynomials.
        /// </summary>
        /// <param name="allTerms"></param>
        /// <returns></returns>
        public static string LinComboToString(Tuple<List<List<Tuple<int, int>>>, List<BigRational>> allTerms)
        {
            List<List<Tuple<int, int>>> terms = allTerms.Item1;
            List<BigRational> coefficients = allTerms.Item2;

            if (terms.Count == 0)
                return "0";

            string s = "";
            s += TermToString(terms[terms.Count - 1], coefficients[coefficients.Count - 1]);

            for (int i = terms.Count - 2; i >= 0; i--)
                s += " + " + TermToString(terms[i], coefficients[i]);

            return s;
        }

        /// <summary>
        /// Sends an individual term (product of choose polys) to a string
        /// helper function for LinComboToString
        /// </summary>
        /// <param name="term"></param>
        /// <param name="coef"></param>
        /// <returns></returns>
        private static string TermToString(List<Tuple<int, int>> term, BigRational coefficient)
        {
            if (term.Count == 0)
                return coefficient.ToString();

            string s = coefficient.ToString();
            foreach (Tuple<int, int> t in term)
            {
                s += " * (p'_" + t.Item1 + " C " + t.Item2 + ")";
            }

            return s;
        }


        /// <summary>
        /// Expresses a symmetric polynomial in the basis (p'_k C j)*...*(p'_k C j) as a polynomial with the first numTerms terms
        /// by the following conversion:
        /// (p'_k C j) -> (moebiussum(k) C j) * (1 / (1 + q^k))^j = (moebiussum(k) C j) * (q^(-k) - q^(-2k) + q^(3k) - ...)^j
        /// The actual polynomial is a power series, but we are concerned with
        /// a finite number of the terms. 
        /// </summary>
        /// <param name="allTerms"></param>
        /// <returns></returns>
        public static LaurentPolynomial CyclicPolynomialBasisToPolynomial(Tuple<List<List<Tuple<int, int>>>, List<BigRational>> allTerms, int lowestDegree)
        {
            List<List<Tuple<int, int>>> symPolys = allTerms.Item1;
            List<BigRational> coefficients = allTerms.Item2;

            LaurentPolynomial output = new(); 

            for (int i = 0; i < symPolys.Count; i++)
            {
                LaurentPolynomial nextTerm = new(coefficients[i]);
                for (int j = 0; j < symPolys[i].Count; j++)
                {
                    //the choose portion of the expression
                    nextTerm *= (LaurentPolynomial)Polynomial.Choose(Polynomial.MoebiusSum(symPolys[i][j].Item1), symPolys[i][j].Item2);
                }

                // multiply by 1 - q^{-1}
                nextTerm *= new LaurentPolynomial(-1, new() { -1, 1 });

                //we now multiply by the desired power series, which we do at the end to avoid
                //saving more terms of the series than necessary
                for (int j = 0; j < symPolys[i].Count; j++)
                {
                    if (symPolys[i][j].Item1 == 0)
                        continue;
                    //now we multiply by the powerseries (out to n terms)
                    nextTerm = MultByPowerSeries(nextTerm, lowestDegree, symPolys[i][j]);
                }

                output += nextTerm;
            }

            return output;
        }

        /// <summary>
        /// Multiplies poly by the power series where k = powerSeries.Item1, j = powerSeries.Item2 -> (q^{-k} / (1 + q^{-k}))^j
        /// This is equivalent to multiplying by the power series (q^{-k} - q^{-2k} + q^{-3k} - ...) j times.
        /// </summary>
        /// <param name="nextTerm"></param>
        /// <param name="numTerms"></param>
        /// <param name="tuple"></param>
        /// <returns></returns>
        private static LaurentPolynomial MultByPowerSeries(LaurentPolynomial poly, int lowestDegree, Tuple<int, int> powerSeries)
        {
            int k = powerSeries.Item1;
            int j = powerSeries.Item2;

            // suppose f = poly
            LaurentPolynomial result = new(poly);
            // safe to round here, since we are not multiplying by any q with positive degree
            result.RoundToNthDegree(lowestDegree);

            for (int a = 0; a < j; a++)
            {
                //product is the polynomial -q^{-k}
                LaurentPolynomial product = new(-k, new List<BigRational> { -1 });
                result *= -1 * product;
                result.RoundToNthDegree(lowestDegree);

                //we save computation by multiplying and adding terms piece by piece
                LaurentPolynomial nextTerm = new(result);
                while (!nextTerm.IsZero())
                {
                    // next term = f*(-q)^{-lk}
                    nextTerm *= product;

                    //we only need to keep track of the terms that have terms with degree bigger than the lowest degree
                    nextTerm.RoundToNthDegree(lowestDegree);

                    // result = f*(q^{-k} - q^{-2k} + ... - (-q)^{-lk})
                    result += nextTerm;
                }
            }

            return result;
        }
    }
}
