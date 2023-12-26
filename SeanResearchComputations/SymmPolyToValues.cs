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
                List<Tuple<int, int>> nextTerm = expression.Item1[i];
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
        /// Returns a string of a symmetric polynomial
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
        /// Given a list of choose basis polynomials of the form (i_1 C j_1)(i_2 C j_2)...(i_n C j_n), sends each term P
        /// to the polynomial statistic \lim_{n\to\infty} q^{-n} \sum_{f \in Conf_n(\F_q)} P(f)
        /// expressed a power series in q^{-1} with coefficients computed accurately out to q^{-lowestDegree}
        /// 
        /// This is given by the following formula:
        /// (x_1 C j_1)(x_2 C j_2)...(x_n C j_n) -> (1 - q^{-1})\prod_{k = 1}^n ((Polynomial.MoebiusSum(i_k, j_k) C j_k) (q^{-i_k} - q^{-2i_k} + ...)^{j_k}
        /// </summary>
        /// <param name="allTerms"></param>
        /// <param name="lowestDegree"></param>
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
                    LaurentPolynomial nextProd = (LaurentPolynomial)Polynomial.Choose(Polynomial.MoebiusSum(symPolys[i][j].Item1), symPolys[i][j].Item2);

                    if (symPolys[i][j].Item1 != 0)
                        nextProd = MultByPowerSeries(nextProd, lowestDegree, symPolys[i][j]);

                    nextTerm *= nextProd;
                }

                output += nextTerm;
            }

            // multiply by 1 - q^{-1}
            output *= new LaurentPolynomial(-1, new() { -1, 1 });
            output.RoundToNthDegree(lowestDegree);

            return output;
        }

        /// <summary>
        /// Multiplies poly by the power series (q^{-k} - q^{-2k} + q^{-3k} - ...)^j times.
        /// Computes the coefficients of the resulting power series accurately to degree `lowestDegree'.
        /// </summary>
        /// <param name="poly"></param>
        /// <param name="lowestDegree"></param>
        /// <param name="powerSeries"></param>
        /// <returns></returns>
        private static LaurentPolynomial MultByPowerSeries(LaurentPolynomial poly, int lowestDegree, Tuple<int, int> powerSeries)
        {
            int i = powerSeries.Item1;
            int j = powerSeries.Item2;

            // compute coefficients of (q^{-i} - q^{-2i} + ... )^j down to lowestDegree - deg(poly) 
            LaurentPolynomial powSeries = new();
            for (int l = j; -l * i >= lowestDegree - poly.Degree(); l++)
            {
                // coefficient of q^{-i*l} in (q^{-i} - q^{-2i} + ... )^j
                BigRational coef = 0;

                foreach (List<int> kpart in IntegerFunctions.KPartitions(l, j))
                {
                    List<int> cycles = IntegerFunctions.PartitionToNumCycles(kpart);

                    BigRational term = IntegerFunctions.MultinomialCoef(j, cycles);

                    int sign = 1;
                    foreach (int n in kpart)
                        if (n % 2 == 0)
                            sign *= -1;
                    term *= sign;

                    coef += term;
                }

                powSeries += new LaurentPolynomial(-l * i, new() { coef });
            }

            // multiply
            LaurentPolynomial output = powSeries * poly;
            output.RoundToNthDegree(lowestDegree);
            return output;
        }
    }
}
