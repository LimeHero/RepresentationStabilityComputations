using Polynomials; // includes symmetric, rational (polynomial field), and laurent polynomials
using IntegerMethods;
using System;
using System.Collections.Generic;
using RepTheory;

namespace RepStabilityComputations;
public class MainMethod
{
    /// <summary>
    /// Execution Point
    /// </summary>
    /// <param name="args"></param>
    static void Main(string[] args)
    {
        int lowestDegree = -30;

        // Print power series for all young tableaus with 1 <= i < a boxes
        int a = 7;
        for (int i = 1; i < a; i++)
        {
            foreach (List<int> part in IntegerFunctions.AllPartitions(i))
            {
                PrintList(part);
                Console.WriteLine(YoungToPoly(part, lowestDegree).ToRevString());
                Console.WriteLine();
            }
        }

        int b = 0;
        // wedge powers
        for (int i = 1; i < b; i++)
        {
            Console.WriteLine(i + "th wedge power");
            Console.WriteLine(WedgePowers(i, lowestDegree).ToRevString());
            Console.WriteLine();
        }
    }


    /// <summary>
    /// Saves the coefficients of the power series in q^{-1} to a results.txt file on the desktop.
    /// Saved in a csv format. All coefficients are taken to be positive (since alternate in sign predictably)
    /// Computed for all young tableau with a number of boxes k satisfying a <= k < b
    /// </summary>
    /// <param name="a"></param>
    /// <param name="b"></param>
    public static void SaveToResultsTxt(int a, int b)
    {
        List<string> lines = new();
        for (int i = a; i < b; i++)
        {
            foreach (List<int> part in IntegerFunctions.AllPartitions(i))
            {
                PrintList(part);
                LaurentPolynomial term = YoungToPoly(part, -5);

                string line = "\"["; for (int j = 0; j < part.Count - 1; j++) line += part[j] + ", ";
                line += part[^1] + "]\",";

                for (int j = 0; j > -30; j--)
                    line += term[j].abs() + ",";

                line += term[-30].abs();

                lines.Add(line);
            }
        }

        // Set a variable to the Documents path.
        string docPath =
          Environment.GetFolderPath(Environment.SpecialFolder.Desktop);

        // Write the string array to a new file named "WriteLines.txt".
        using (StreamWriter outputFile = new StreamWriter(Path.Combine(docPath, "results.txt")))
        {
            foreach (string line in lines)
                outputFile.WriteLine(line);
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
