using IntegerMethods;
using Polynomials;

namespace RepTheory

{
    public class RepTheoryAlgs
    {
        /// <summary>
        /// Computes the character of the representation defined by the given tableau (descending
        /// list of positive integers) evaluated on the given conjugacy class. The cycleType should
        /// be a list of positive integers [i_1, i_2, ..., i_n] such that the cycle type has
        /// i_j j-cycles. Used for testing.
        /// </summary>
        /// <param name="tableau"></param>
        /// <param name="cycles"></param>
        /// <returns></returns>
        public static BigRational FrobeniusFormula(List<int> tableau, List<int> cycles)
        {
            SymmetricPolynomial prod = new SymmetricPolynomial(1);

            for (int i = 0; i < cycles.Count; i++)
                for (int j = 0; j < cycles[i]; j++)
                    prod *= SymmetricPolynomial.Power(i + 1);

            List<int> l = new();
            for (int i = 0; i < tableau.Count; i++)
            {
                l.Add(tableau[i] + tableau.Count - i - 1);
            }

            BigRational rslt = 0;
            // we add each term of the discriminant bit by bit - start with (-x_2)*(-x_3)*...
            // and each binary digit corresponds to a choice of a term in the discriminant
            int numterms = (tableau.Count * (tableau.Count - 1)) / 2;
            int power2 = (int) IntegerFunctions.Pow(2, numterms);

            for (int i = 0; i < power2; i++)
            {
                List<int> currentterm = new(l);
                List<int> bindig = IntegerFunctions.BinaryDigits(i);


                while (bindig.Count < numterms)
                    bindig.Add(0);

                int k = 0;
                int pm1 = 1;
                for (int n = 0; n < tableau.Count; n++)
                {
                    for (int m = 0; m < n; m++)
                    {
                        if (bindig[k] == 0)
                        {
                            pm1 *= -1;
                            currentterm[n] -= 1;
                        }
                        else
                            currentterm[m] -= 1;

                        k++;
                    }
                }

                currentterm.Sort();
                currentterm.Reverse();
                if (currentterm[currentterm.Count - 1] < 0)
                    continue;

                while (currentterm[currentterm.Count - 1] == 0)
                    currentterm.RemoveAt(currentterm.Count - 1);

                int index = SymmetricPolynomial.FindMonomialIndex(currentterm);
                if (index >= prod.coefs.Count)
                    continue;

                rslt += prod.coefs[index] * pm1;
            }

            return rslt;
        }

        /// <summary>
        /// Given a representation with character given by "terms", a polynomial in the number of each cycle,
        /// return the character evaluated on a certain set of cycles. A simple computation that is effectively a simple
        /// form of frobenius for known characters. Used for testing.
        /// </summary>
        /// <param name="terms"></param>
        /// <param name="cycles"></param>
        /// <returns></returns>
        public static BigRational ChooseFormToCharacter(Tuple<List<List<Tuple<int, int>>>, List<BigRational>> terms,
            List<int> cycles)
        {
            BigRational rslt = 0;
            for (int i = 0; i < terms.Item1.Count; i++)
            {
                BigRational nextterm = terms.Item2[i];

                for (int j = 0; j < terms.Item1[i].Count; j++)
                {
                    int k = terms.Item1[i][j].Item1;
                    if (k == 0)
                        continue;

                    while (cycles.Count < k)
                        cycles.Add(0);

                    nextterm *= IntegerFunctions.Choose(cycles[k - 1], terms.Item1[i][j].Item2);
                }

                rslt += nextterm;
            }

            return rslt;
        }

        /// <summary>
        /// Prints a list in a nice format
        /// </summary>
        /// <param name="list"></param>s
        public static void PrintList<T>(List<T> list)
        {
            if (list.Count == 0)
            {
                Console.WriteLine("[]");
                return;
            }

            Console.Write("[");
            for (int j = 0; j < list.Count - 1; j++)
            {
                Console.Write(list[j] + ", ");
            }
            Console.WriteLine(list[list.Count - 1] + "]");
        }
    }
}