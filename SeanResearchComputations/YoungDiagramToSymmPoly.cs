using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using IntegerMethods;
using Polynomials;

namespace SeanResearchComputations
{
    public class YoungDiagramToSymmPoly
    {
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
        /// Returns the symmetric polynomial (in choose basis) of the young diagram corresponding to (n, k).
        /// </summary>
        /// <param name="k"></param>
        /// <param name="printCyclicBasis"></param>
        /// <param name="smallestDegreePrinted"></param>
        public static Tuple<List<List<Tuple<int, int>>>, List<BigRational>> YoungTwoRowsToChoose(int k)
        {
            Tuple<List<List<Tuple<int, int>>>, List<BigRational>> result = new(new(), new());

            foreach (List<int> part in IntegerFunctions.AllPartitions(k))
            {
                result.Item1.Add(new());
                result.Item2.Add(1);

                for (int i = 0; i < part.Count; i++)
                {
                    int m = part[i];
                    int j = 0; // number of elements equal to m
                    while (i < part.Count && m == part[i])
                    {
                        i++;
                        j++;
                    }
                    i--; // so we dont skip over next m



                    result.Item1[result.Item1.Count - 1].Add(new(m, j));
                }
            }

            foreach (List<int> part in IntegerFunctions.AllPartitions(k - 1))
            {
                result.Item1.Add(new());
                result.Item2.Add(-1); // only difference

                for (int i = 0; i < part.Count; i++)
                {
                    int m = part[i];
                    int j = 0; // number of elements equal to m
                    while (i < part.Count && m == part[i])
                    {
                        i++;
                        j++;
                    }
                    i--; // so we dont skip over next m

                    result.Item1[result.Item1.Count - 1].Add(new(m, j));
                }
            }

            return result;
        }

        /// <summary>
        /// Returns the symmetric polynomial (in choose basis) of the young diagram corresponding to (n, k1, k2).
        /// </summary>
        /// <param name="k1"></param>
        /// <param name="k2"></param>
        /// <param name="printCyclicBasis"></param>
        /// <param name="smallestDegreePrinted"></param>
        public static Tuple<List<List<Tuple<int, int>>>, List<BigRational>> YoungThreeRowsToChoose(int k1, int k2)
        {
            Tuple<List<List<Tuple<int, int>>>, List<BigRational>> result = new(new(), new());

            //(x_1 - x_2)(x_1 - x_3)(x_2 - x_3) = (x_1^2x_2 + x_1x_3^2 + x_2^2x_3) - (x_1^2x_3 + x_1x_2^2 + x_2x_3^2)
            List<Tuple<int, int[]>> discriminant = new();
            // we need not add the x_1 coefficient, since tending toward infinity
            discriminant.Add(new(1, new int[] { 1, 0 }));
            discriminant.Add(new(1, new int[] { 0, 2 }));
            discriminant.Add(new(1, new int[] { 2, 1 }));

            discriminant.Add(new(-1, new int[] { 0, 1 }));
            discriminant.Add(new(-1, new int[] { 2, 0 }));
            discriminant.Add(new(-1, new int[] { 1, 2 }));

            // all terms of the discriminant
            foreach (Tuple<int, int[]> term in discriminant)
            {
                // [l_2 - discriminant_{x_2}, l_3 - discriminant_{x_3}]
                int[] l = new int[] { k1 - term.Item2[0] + 1, k2 - term.Item2[1] }; // [a, b]

                if (l[0] < 0 || l[1] < 0) continue;

                foreach (List<int> parta in IntegerFunctions.AllPartitions(l[0]))
                {
                    foreach (List<int> partb in IntegerFunctions.AllPartitions(l[1]))
                    {

                        List<int> cyclesa = IntegerFunctions.PartitionToNumCycles(parta);
                        List<int> cyclesb = IntegerFunctions.PartitionToNumCycles(partb);

                        result.Item1.Add(new()); // next term added
                        result.Item2.Add(term.Item1); // coefficient of determinant term

                        int maxterm = Math.Max(parta[0], partb[0]);
                        while (cyclesa.Count < maxterm)
                            cyclesa.Add(0);
                        while (cyclesb.Count < maxterm)
                            cyclesb.Add(0);

                        for (int i = 0; i < maxterm; i++)
                        {
                            if (cyclesa[i] == 0 && cyclesb[i] == 0)
                                continue;

                            result.Item1[^1].Add(new(i + 1, cyclesa[i] + cyclesb[i])); // [^1] means 1 from the end of list (last elmnt)
                            result.Item2[^1] *= IntegerFunctions.Choose(cyclesa[i] + cyclesb[i], cyclesa[i]);
                        }
                    }
                }
            }

            return result;
        }

        /// <summary>
        /// Returns the symmetric polynomial (in choose basis) of the young diagram corresponding to (n, k[0], k[1], ...)
        /// </summary>
        /// <param name="k"></param>
        /// <returns></returns>
        public static Tuple<List<List<Tuple<int, int>>>, List<BigRational>> YoungDiagramToChoose(List<int> k)
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
            // we add each term of the discriminant bit by bit - start with (-x_1)*(-x_2)*...
            // and each binary digit corresponds to a choice of a term in the discriminant
            int numterms = k.Count * (k.Count + 1) / 2;
            int power2 = (int)IntegerFunctions.Pow(2, numterms);
            // compute determinant
            // TODO: make wayyy more efficient, combine terms, and remove 0 terms
            for (int i = 0; i < power2; i++)
            {
                List<int> currentterm = new() { 0 }; for (int j = 0; j <= k.Count; j++) currentterm.Add(0);
                List<int> bindig = IntegerFunctions.binaryDigits(i);

                while (bindig.Count < numterms)
                    bindig.Add(0);

                int r = 0;
                int pm1 = 1;
                for (int n = 1; n < k.Count + 1; n++)
                {
                    for (int m = 0; m < n; m++)
                    {
                        if (bindig[r] == 0)
                        {
                            pm1 *= -1;
                            currentterm[n - 1] += 1;
                        }
                        else
                            if (m > 0)
                            currentterm[m - 1] += 1;

                        r++;
                    }
                }

                discriminant.Add(new(pm1, currentterm));
            }

            // iterate through all terms of the discriminant
            foreach (Tuple<int, List<int>> term in discriminant)
            {
                // the integers l in Frobenius formula
                List<int> l = new(k);
                for (int i = 0; i < l.Count; i++)
                {
                    l[i] -= term.Item2[i];
                    l[i] += l.Count - 1 - i;
                }

                // if any l[i] is negative, the coefficient is zero.
                foreach (int b in l) if (b < 0) continue;

                // All partitions of l[0], l[1], ... l[^1]
                foreach (List<List<int>> parts in IntegerFunctions.AllPartitionLists(l))
                {
                    // Change to cycle format of partitions
                    List<List<int>> cycles = new();
                    foreach (List<int> part in parts) cycles.Add(IntegerFunctions.PartitionToNumCycles(part));

                    result.Item1.Add(new()); // next term added
                    result.Item2.Add(term.Item1); // coefficient of determinant term (non-zero)

                    // make all partitions the same size
                    int maxterm = parts[0][0];
                    for (int i = 1; i < parts.Count; i++)
                        if (parts[i][0] > maxterm)
                            maxterm = parts[i][0];

                    foreach (List<int> cycle in cycles)
                        while (cycle.Count < maxterm)
                            cycle.Add(0);

                    for (int i = 0; i < maxterm; i++)
                    {
                        // Finding the (x_i C j_i) term
                        int totalsum = 0;
                        foreach (List<int> cycle in cycles)
                            totalsum += cycle[i];

                        // j_i = 1, so skip
                        if (totalsum == 0)
                            continue;

                        // i + 1, since we index from 1 in the variables x_1, ... x_n for character polynomial
                        result.Item1[^1].Add(new(i + 1, totalsum)); // [^1] means 1 from the end of list (last elmnt)

                        List<int> ithterms = new();
                        foreach (List<int> cycle in cycles)
                            ithterms.Add(cycle[i]);
                        result.Item2[^1] *= IntegerFunctions.MultinomialCoef(totalsum, ithterms);
                    }
                }
            }

            return result;
        }
    }
}
