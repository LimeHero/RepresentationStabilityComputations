using Microsoft.VisualStudio.TestTools.UnitTesting;
using System.Collections.Generic;
using IntegerMethods;
using RepTheory;
using System;
using SeanResearchComputations;
using System.Linq;

namespace RepTheoryTests
{
    [TestClass]
    public class RepTheoryAlgsTests
    {
        [TestMethod]
        public void FrobeniusFormulaTest1()
        {
            // testing all cycle types for this tableau
            List<int> tableau = new() { 3, 2 };

            List<int> part = new() { 5 };
            Assert.IsTrue(RepTheoryAlgs.FrobeniusFormula(tableau, part) == 5);
            Console.WriteLine(RepTheoryAlgs.FrobeniusFormula(tableau, part));

            part = new() { 3, 1 };
            Assert.IsTrue(RepTheoryAlgs.FrobeniusFormula(tableau, part) == 1);

            part = new() { 2, 0, 1 };
            Assert.IsTrue(RepTheoryAlgs.FrobeniusFormula(tableau, part) == -1);

            part = new() { 1, 0, 0, 1 };
            Assert.IsTrue(RepTheoryAlgs.FrobeniusFormula(tableau, part) == -1);

            part = new() { 1, 2 };
            Assert.IsTrue(RepTheoryAlgs.FrobeniusFormula(tableau, part) == 1);

            part = new() { 0, 0, 0, 0, 1 };
            Assert.IsTrue(RepTheoryAlgs.FrobeniusFormula(tableau, part) == 0);

            part = new() { 0, 1, 1 };
            Assert.IsTrue(RepTheoryAlgs.FrobeniusFormula(tableau, part) == 1);
        }


        [TestMethod]
        public void FrobeniusFormulaTest2()
        {
            // testing all cycle types for this tableau using character = (i_1-1)(i_1-2)/2 + i_2 - 1
            List<int> tableau1 = new() { 3, 2 };
            List<int> tableau2 = new() { 7, 2 };
            foreach (List<int> part in IntegerFunctions.AllPartitions(5))
            {
                List<int> cycles = IntegerFunctions.PartitionToNumCycles(part);

                if (cycles.Count == 1)
                    cycles.Add(0);

                Assert.IsTrue(RepTheoryAlgs.FrobeniusFormula(tableau1, cycles) ==
                    ((cycles[0] - 1) * (cycles[0] - 2)) / 2 + cycles[1] - 1);
            }


            foreach (List<int> part in IntegerFunctions.AllPartitions(9))
            {
                List<int> cycles = IntegerFunctions.PartitionToNumCycles(part);

                if (cycles.Count == 1)
                    cycles.Add(0);

                Assert.IsTrue(RepTheoryAlgs.FrobeniusFormula(tableau2, cycles) ==
                    ((cycles[0] - 1) * (cycles[0] - 2)) / 2 + cycles[1] - 1);
            }
        }

        [TestMethod]
        public void FrobeniusFormulaTest3()
        {
            // testing frob using known symmetric polynomial form of wedge of standard
            // and converting it into a formula of number of cycles
            List<int> tableau = new() { 2, 1, 1, 1, 1 };

            Tuple<List<List<Tuple<int, int>>>, List<BigRational>> terms =
                SymmPolyToValues.SymmetricInChooseBasis(YoungDiagramToSymmPoly.WedgeToSymmetricPolynomial(tableau.Count - 1));

            foreach (List<int> part in IntegerFunctions.AllPartitions(tableau.Sum()))
            {
                List<int> cycles = IntegerFunctions.PartitionToNumCycles(part);

                BigRational rslt = RepTheoryAlgs.FrobeniusFormula(tableau, cycles);

                Assert.IsTrue(RepTheoryAlgs.ChooseFormToCharacter(terms, cycles) == rslt);
            }
        }
    }
}