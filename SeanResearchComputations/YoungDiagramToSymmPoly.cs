using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using IntegerMethods;
using Polynomials;

namespace RepStabilityComputations
{
    public class YoungDiagramToSymmPoly
    {
        private static List<List<int>> yng_to_choose_keys = new();
        private static List<Tuple<List<List<Tuple<int, int>>>, List<BigRational>>> yng_to_choose_values = new();

        /// <summary>
        /// Returns the symmetric polynomial e_n - e_{n-1} + e_{n-2} - ... +- 1
        /// </summary>
        /// <param name="n"></param>
        /// <returns></returns>
        public static SymmetricPolynomial WedgeToSymmetricPolynomial(int n)
        {
            SymmetricPolynomial p = new SymmetricPolynomial();
            int pm1 = 1;
            for (int i = n; i >= 1; i--)
            {
                p += pm1 * SymmetricPolynomial.Elementary(i);
                pm1 *= -1;
            }
            p += new SymmetricPolynomial(pm1);

            return p;
        }

        /// <summary>
        /// Returns the characteristic polynomial (in choose basis) of the family of irreducible representations
        /// given by young diagrams of the form (n, k[0], k[1], ... k[^1]) for n -> infty
        /// </summary>
        /// <param name="k"></param>
        /// <returns></returns>
        public static Tuple<List<List<Tuple<int, int>>>, List<BigRational>> YoungDiagramToChoose(List<int> k, bool memoized = true)
        {
            // special cases
            {
                bool allones = true;
                for (int i = 0; i < k.Count(); i++)
                {
                    if (k[i] <= 0)
                        throw new Exception("(YoungToPoly) No negative numbers!");

                    if (i > 0 && k[i] > k[i - 1])
                        throw new Exception("(YoungToPoly) Must be nonincreasing sequence!");

                    if (k[i] != 1)
                        allones = false;
                }

                if (allones)
                    return SymmPolyToValues.SymmetricInChooseBasis(WedgeToSymmetricPolynomial(k.Count()));
            }

            Tuple<List<List<Tuple<int, int>>>, List<BigRational>> result = new(new(), new());

            List<Tuple<int, List<int>>> discriminant = new();
            // we need not add the x_0 coefficient, since tending toward infinity
            // we compute the discriminant binomials by taking determinant of vandermont matrix

            // we have to compute this sgn_change to account for the fact that we have the determinant
            // as \prod_{i<j} (x_i - x_j) in Frobenius formula, but as \prod_{i<j} (x_j - x_i) for det Vandermont matrix
            int sgn_change = ((k.Count * (k.Count + 1) / 2) % 2) == 0 ? 1 : -1;
            foreach (List<int> L in IntegerFunctions.Permutations(k.Count + 1))
            {
                List<int> currentterm = new() { }; for (int j = 0; j <= k.Count; j++) currentterm.Add(L[j]);
                discriminant.Add(new(IntegerFunctions.Sgn(L) * sgn_change, currentterm));
            }

            // iterate through all terms of the discriminant
            foreach (Tuple<int, List<int>> term in discriminant)
            {
                // the integers l in Frobenius formula
                List<int> l = new(k);
                for (int i = 0; i < l.Count; i++)
                {
                    l[i] -= term.Item2[i + 1];
                    l[i] += l.Count - 1 - i;
                }

                //memoization lookup
                if (memoized) {
                    // check if this term has been done already
                    int equal_index = -1;
                    bool equal = true;
                    for (int a = 0; a < yng_to_choose_keys.Count; a++)
                    {
                        equal = true;
                        if (yng_to_choose_keys[a].Count != l.Count)
                            continue;

                        for (int b = 0; b < yng_to_choose_keys[a].Count; b++)
                        {
                            if (yng_to_choose_keys[a][b] != l[b])
                            {
                                equal = false;
                                break;
                            }
                        }

                        if (equal)
                        {
                            equal_index = a;
                            break;
                        }
                    }

                    // if so, equal_index = a
                    if (equal_index != -1)
                    {
                        // add to result
                        for (int a = 0; a < yng_to_choose_values[equal_index].Item1.Count; a++)
                        {
                            // result.Item1.Add(yng_to_choose_values[equal_index].Item1[a]);
                            result.Item1.Add(new());
                            for (int b = 0; b < yng_to_choose_values[equal_index].Item1[a].Count; b++)
                                result.Item1[^1].Add(yng_to_choose_values[equal_index].Item1[a][b]);

                            result.Item2.Add(yng_to_choose_values[equal_index].Item2[a] * term.Item1);
                        }

                        continue;
                    }
                }

                // the result for this term of discriminant
                Tuple<List<List<Tuple<int, int>>>, List<BigRational>> this_term_result = new(new(), new());

                // All partitions of l[0], l[1], ... l[^1]
                foreach (List<List<int>> parts in IntegerFunctions.AllPartitionLists(l))
                {
                    // Change to cycle format of partitions
                    List<List<int>> cycles = new();
                    foreach (List<int> part in parts) cycles.Add(IntegerFunctions.PartitionToNumCycles(part));

                    // make all partitions the same size for ease of implementation
                    int maxterm = parts[0][0];
                    for (int i = 1; i < parts.Count; i++)
                        if (parts[i][0] > maxterm)
                            maxterm = parts[i][0];
                    foreach (List<int> cycle in cycles)
                        while (cycle.Count < maxterm)
                            cycle.Add(0);

                    List<Tuple<int, int>> choosePoly = new(); // next term added
                    BigRational coef = term.Item1; // coefficient of determinant term (\pm 1)

                    for (int i = 0; i < maxterm; i++)
                    {
                        // Finding the (x_i C j_i) term
                        int totalsum = 0;
                        foreach (List<int> cycle in cycles)
                            totalsum += cycle[i];

                        // j_i = 0, so skip
                        if (totalsum == 0)
                            continue;

                        // i + 1, since we index from 1 in the variables x_1, ... x_n for character polynomial
                        choosePoly.Add(new(i + 1, totalsum)); // [^1] means 1 from the end of list (last elmnt)

                        // (multinomial) coefficient
                        List<int> ithterms = new();
                        foreach (List<int> cycle in cycles)
                            ithterms.Add(cycle[i]);
                        coef *= IntegerFunctions.MultinomialCoef(totalsum, ithterms);
                    }

                    this_term_result.Item1.Add(choosePoly);
                    this_term_result.Item2.Add(coef);
                }

                // memoization saving
                if (memoized) {
                    // yng_to_choose_keys.Add(l);
                    List<int> l_copy = new(l);
                    yng_to_choose_keys.Add(l_copy);

                    // yng_to_symm_values.Add(this_term_result);
                    Tuple<List<List<Tuple<int, int>>>, List < BigRational >> this_term_result_copy = new(new(), new());
                    for (int a = 0; a < this_term_result.Item1.Count; a++)
                    {
                        this_term_result_copy.Item2.Add(this_term_result.Item2[a] * term.Item1);
                        // this_term_result_copy.Item1.Add(this_term_result.Item1[a]);
                        List<Tuple<int, int>> cp = new();
                        for (int b = 0; b < this_term_result.Item1[a].Count; b++)
                            cp.Add(this_term_result.Item1[a][b]);
                        this_term_result_copy.Item1.Add(this_term_result.Item1[a]);
                    }
                    yng_to_choose_values.Add(this_term_result_copy);
                }

                // add to final result
                for (int a = 0; a < this_term_result.Item1.Count; a++)
                {
                    // result.Item1.Add(this_term_result.Item1[a]);
                    result.Item1.Add(new());
                    for (int b = 0; b < this_term_result.Item1[a].Count; b++)
                        result.Item1[^1].Add(this_term_result.Item1[a][b]);
                    result.Item2.Add(this_term_result.Item2[a]);
                }
            }

            return result;
        }
    }
}
