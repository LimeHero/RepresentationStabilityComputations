using Polynomials; // includes symmetric, rational (polynomial field), and laurent polynomials
using IntegerMethods;
using System;
using System.Collections.Generic;
using RepTheory;

namespace SeanResearchComputations;
public class MainMethod
{
    /// <summary>
    /// Execution Point
    /// </summary>
    /// <param name="args"></param>
    static void Main(string[] args)
    {

        // Returns CTP((p'_k j)) for various j + k = i
        for (int i = 1; i < 0; i++)
        {
            for (int j = 1; j < 6; j++)
            {
                int k = i - j;
                if (k <= 0)
                    break;

                Console.WriteLine("(p'_" + k + " " + j + ")");
                Tuple<List<List<Tuple<int, int>>>, List<BigRational>> elt =
                    new(new() { new() { new(k, j) } }, new() { 1 });
                LaurentPolynomial series = SymmPolyToValues.CyclicPolynomialBasisToPolynomial(elt, -10);
                Console.WriteLine(series.ToRevString());
            }
        }

        // Frobenius formula for all young tableaus with k boxes
        {
            int k = -1;
            bool include_top_row_n = true;
            int n = 4; // top row length, if included

            if (k > 0)
            {

                if (!include_top_row_n)
                    n = 0;

                k += n;

                Console.Write("|");
                foreach (List<int> part in IntegerFunctions.AllPartitions(k))
                {
                    char letter = 'a';
                    for (int j = 0; j < part.Count; j++)
                    {
                        if (part[j] == 1 && j == 0) Console.Write("id");
                        if (part[j] == 1) break;

                        Console.Write("(");
                        for (int i = 0; i < part[j]; i++)
                        {
                            Console.Write(letter);
                            letter++;
                        }
                        Console.Write(")");
                    }
                    Console.Write("|");
                }
                Console.WriteLine("");
                foreach (List<int> tableau in IntegerFunctions.AllPartitions(k - n))
                {
                    if (include_top_row_n)
                        tableau.Insert(0, n);

                    PrintList(tableau);
                    foreach (List<int> part in IntegerFunctions.AllPartitions(k))
                    {
                        List<int> cycletype = IntegerFunctions.PartitionToNumCycles(part);
                        Console.Write("|");
                        string val = RepTheoryAlgs.FrobeniusFormula(tableau, cycletype).ToString();
                        Console.Write(val);
                        int len = 0;
                        for (int j = 0; j < part.Count; j++)
                        {
                            if (j == 0 && part[j] == 1) len = 2;

                            if (part[j] == 1) break;

                            len += 2 + part[j];
                        }

                        for (int j = 0; j < len - val.Length; j++)
                            Console.Write(' ');
                    }
                    Console.WriteLine("|");
                }
            }
        }

        // Characters of young tableaus with k boxes (top row tending toward infinity)
        // in cycle choose notation
        {
            int k = -1;

            if (k > 0)
            {
                foreach (List<int> part in IntegerFunctions.AllPartitions(k))
                {

                    Tuple<List<List<Tuple<int, int>>>, List<BigRational>> rslt = YoungDiagramToSymmPoly.YoungDiagramToChoose(part);

                    bool is_wedge_sum = true;
                    foreach (int i in part)
                    {
                        if (i != 1)
                        {
                            is_wedge_sum = false;
                            break;
                        }
                    }
                    if (is_wedge_sum)
                    {
                        PrintList(part);
                        Tuple<List<List<Tuple<int, int>>>, List<BigRational>> flip_order = new(new(), new());
                        for (int i = rslt.Item1.Count - 1; i >= 0; i--)
                        {
                            flip_order.Item1.Add(rslt.Item1[i]);
                            flip_order.Item2.Add(rslt.Item2[i]);
                        }
                        Console.WriteLine(SymmPolyToValues.LinComboToString(flip_order));
                        Console.WriteLine();
                        continue;
                    }

                    // reduce:
                    Tuple<List<List<Tuple<int, int>>>, List<BigRational>> reduced = new(new(), new());
                    for (int j = 0; j <= k; j++)
                    {
                        foreach (List<int> p in IntegerFunctions.AllPartitions(j))
                        {
                            List<int> cycles = IntegerFunctions.PartitionToNumCycles(p);
                            List<Tuple<int, int>> term = new();
                            for (int i = 0; i < cycles.Count; i++)
                            {
                                if (cycles[i] == 0)
                                    continue;
                                term.Add(new Tuple<int, int>(i + 1, cycles[i]));
                            }

                            reduced.Item1.Add(term);
                            reduced.Item2.Add(0);
                        }
                    }

                    for (int j = 0; j < rslt.Item1.Count; j++)
                    {
                        bool same = false;
                        for (int i = 0; i < reduced.Item1.Count; i++)
                        {
                            for (int l = 0; l < rslt.Item1[j].Count; l++)
                            {
                                if (rslt.Item1[j].Count != reduced.Item1[i].Count)
                                    break;
                                if (rslt.Item1[j][l].Item1 == reduced.Item1[i][l].Item1 &&
                                    rslt.Item1[j][l].Item2 == reduced.Item1[i][l].Item2)
                                {
                                    same = true;
                                    reduced.Item2[i] += rslt.Item2[j];
                                    break;
                                }
                            }

                        }
                        if (same == false)
                            Console.WriteLine("Bad ordering, something bad :(");
                    }

                    PrintList(part);
                    Console.WriteLine(SymmPolyToValues.LinComboToString(reduced));
                    Console.WriteLine();
                }
            }
        }

        // All young tableaus with i boxes to LaTeX table
        {
            int a = 0;
            int b = 0;

            int leastDegree = -30;
            int perPage = 4;
            int firstColWidth = 4; //cm
            int secColWidth = 11; //cm
            if (a < b)
            {
                Console.WriteLine("\\begin{center}");
                Console.WriteLine("\\begin{tabular}{||Y{" + firstColWidth + "cm}|Y{" + secColWidth + "cm}|}");
            }
            int perPageCount = 0;
            for (int i = a; i < b; i++)
            {
                foreach (List<int> part in IntegerFunctions.AllPartitions(i))
                {
                    // DELTE LATER
                    if (part[0] != 1)
                        continue;

                    if (perPageCount >= perPage)
                    {
                        perPageCount = 0;
                        Console.WriteLine("\\end{tabular}");
                        Console.WriteLine("\\end{center}");
                        Console.WriteLine("\\begin{center}");
                        Console.WriteLine("\\begin{tabular}{||Y{" + firstColWidth + "cm}|Y{" + secColWidth + "cm}|}");
                    }
                    Console.Write("\\yng(");
                    for (int j = 0; j < part.Count; j++)
                    {
                        Console.Write(part[j]);
                        if (j != part.Count - 1)
                            Console.Write(",");
                        else
                            Console.WriteLine(") & ");
                    }

                    Console.WriteLine(YoungToPoly(part, leastDegree).ToRevString(true) + " + \\dots");
                    Console.WriteLine("\\\\");
                    Console.WriteLine("\\hline");

                    perPageCount++;
                }
            }
            if (a < b)
            {
                Console.WriteLine("\\end{tabular}");
                Console.WriteLine("\\end{center}");
            }
        }

        // All young tableaus with i boxes
        for (int i = 1; i < 9; i++)
        {
            foreach (List<int> part in IntegerFunctions.AllPartitions(i))
            {
                PrintList(part);
                Console.WriteLine(YoungToPoly(part, -30).ToRevString());
            }
        }

        // wedge powers
        for (int i = 1; i < 0; i++)
        {
            Console.WriteLine(i + "th wedge power");
            Console.WriteLine(WedgePowers(i, -20).ToRevString());
            Console.WriteLine("");
        }

        // two rows
        for (int i = 1; i < 0; i++)
        {
            Console.WriteLine(YoungTwoRows(i, -20).ToRevString());
            Console.WriteLine("");
        }

        // three rows (row first)
        for (int i = 1; i < 0; i++)
        {
            for (int j = 1; j <= i; j++)
            {
                Console.WriteLine(i + ", " + j + ":");
                Console.WriteLine(YoungThreeRows(i, j, -20).ToRevString());
                Console.WriteLine("");
            }
        }

        // three rows (col first)
        for (int i = 1; i < 0; i++)
        {
            for (int j = i; j <= 7; j++)
            {
                Console.WriteLine(j + ", " + i + ":");
                Console.WriteLine(YoungThreeRows(j, i, -20).ToRevString());
                Console.WriteLine("");
            }
        }
    }

    /// <summary>
    /// Prints the coefficients for wedge^k V of the standard representation.
    /// </summary>
    /// <param name="k"></param>
    /// <param name="printCyclicBasis"></param>
    /// <param name="smallestDegreePrinted"></param>
    public static LaurentPolynomial WedgePowers(int k, int smallestDegreePrinted = -10)
    {
        return SymmPolyToValues.SymmPolyToPoly(
            YoungDiagramToSymmPoly.WedgeToSymmetricPolynomial(k), smallestDegreePrinted);
    }

    /// <summary>
    /// Prints the coefficients for the young diagram corresponding to (n, k).
    /// </summary>
    /// <param name="k"></param>
    /// <param name="printCyclicBasis"></param>
    /// <param name="smallestDegreePrinted"></param>
    public static LaurentPolynomial YoungTwoRows(int k, int smallestDegreePrinted = -10)
    {
        return SymmPolyToValues.CyclicPolynomialBasisToPolynomial(
            YoungDiagramToSymmPoly.YoungTwoRowsToChoose(k), smallestDegreePrinted);
    }

    /// <summary>
    /// Prints the coefficients for the young diagram corresponding to (n, k1, k2).
    /// </summary>
    /// <param name="k1"></param>
    /// <param name="k2"></param>
    /// <param name="smallestDegreePrinted"></param>
    /// <returns></returns>
    public static LaurentPolynomial YoungThreeRows(int k1, int k2, int smallestDegreePrinted = -10)
    {
        return SymmPolyToValues.CyclicPolynomialBasisToPolynomial(
            YoungDiagramToSymmPoly.YoungThreeRowsToChoose(k1, k2), smallestDegreePrinted);
    }

    /// <summary>
    /// Prints the coefficients for the young diagram corresponding to (n, k[0], k[1], ...).
    /// 
    /// k should thus be a nonincreasing list of positive integers.
    /// 
    /// Calls the wedge powers separately for performance.
    /// </summary>
    /// <param name="k"></param>
    /// <param name="smallestDegreePrinted"></param>
    /// <returns></returns>
    public static LaurentPolynomial YoungToPoly(List<int> k, int smallestDegreePrinted = -10)
    {
        bool allones = true;
        for (int i = 0; i < k.Count(); i++)
        {
            if (k[i] <= 0)
                throw new Exception("(YoungToPoly) No negative numbers!");

            if (i > 0 && k[i] > k[i - 1])
                throw new Exception("(YoungToPoly) Must be decreasing sequence!");

            if (k[i] != 1)
                allones = false;
        }

        if (allones)
            return WedgePowers(k.Count(), smallestDegreePrinted);

        return SymmPolyToValues.CyclicPolynomialBasisToPolynomial(
            YoungDiagramToSymmPoly.YoungDiagramToChoose(k), smallestDegreePrinted);
    }

    /// <summary>
    /// Prints the list in a nice format
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
