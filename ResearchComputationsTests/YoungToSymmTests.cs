using Microsoft.VisualStudio.TestTools.UnitTesting;
using IntegerMethods;
using Polynomials;
using SeanResearchComputations;
using System;
using System.Collections.Generic;
using System.Linq;
using RepTheory;

namespace ResearchComputationsTests
{
    [TestClass]
    public class YoungToSymmTests
    {
        [TestMethod]
        public void WedgeToSymmPolyTest()
        {
            SymmetricPolynomial ex = SymmetricPolynomial.Elementary(4) -
                SymmetricPolynomial.Elementary(3) + SymmetricPolynomial.Elementary(2) - SymmetricPolynomial.Elementary(1) + 1;

            Assert.IsTrue(ex == YoungDiagramToSymmPoly.WedgeToSymmetricPolynomial(4));
        }

        [TestMethod]
        public void WedgeToPolynomialTest1()
        {
            LaurentPolynomial result = SymmPolyToValues.SymmPolyToPoly(YoungDiagramToSymmPoly.WedgeToSymmetricPolynomial(1), -6);

            LaurentPolynomial expected = new(-6, new List<BigRational>() { 2, -2, 2, -2, 2, -1 });

            Assert.IsTrue(result == expected);
        }

        [TestMethod]
        public void TwoRowsTest1()
        {
            LaurentPolynomial result1 = MainMethod.YoungTwoRows(1, -10);
            LaurentPolynomial result2 = MainMethod.WedgePowers(1, -10);

            Assert.IsTrue(result1 == result2);
        }

        [TestMethod]
        public void TwoRowsTest2()
        {
            List<int> tableau = new List<int> { 7, 6 };

            Tuple<List<List<Tuple<int, int>>>, List<BigRational>> terms = YoungDiagramToSymmPoly.YoungTwoRowsToChoose(tableau[1]);

            foreach (List<int> part in IntegerFunctions.AllPartitions(tableau.Sum()))
            {
                List<int> cycles = IntegerFunctions.PartitionToNumCycles(part);

                Assert.IsTrue(RepTheoryAlgs.FrobeniusFormula(tableau, cycles) ==
                    RepTheoryAlgs.ChooseFormToCharacter(terms, cycles));
            }
        }

        [TestMethod]
        public void ThreeRowsTest1()
        {
            List<int> tableau = new List<int> { 4, 2, 1 };

            Tuple<List<List<Tuple<int, int>>>, List<BigRational>> terms =
                YoungDiagramToSymmPoly.YoungThreeRowsToChoose(tableau[1], tableau[2]);

            foreach (List<int> part in IntegerFunctions.AllPartitions(tableau.Sum()))
            {
                List<int> cycles = IntegerFunctions.PartitionToNumCycles(part);

                Assert.IsTrue(RepTheoryAlgs.FrobeniusFormula(tableau, cycles) ==
                    RepTheoryAlgs.ChooseFormToCharacter(terms, cycles));
            }
        }

        [TestMethod]
        public void ThreeRowsTest2()
        {
            List<int> tableau = new List<int> { 8, 5, 3 };

            Tuple<List<List<Tuple<int, int>>>, List<BigRational>> terms =
                YoungDiagramToSymmPoly.YoungThreeRowsToChoose(tableau[1], tableau[2]);

            foreach (List<int> part in IntegerFunctions.AllPartitions(tableau.Sum()))
            {
                List<int> cycles = IntegerFunctions.PartitionToNumCycles(part);

                Assert.IsTrue(RepTheoryAlgs.FrobeniusFormula(tableau, cycles) ==
                    RepTheoryAlgs.ChooseFormToCharacter(terms, cycles));
            }
        }

        [TestMethod]
        public void ThreeRowsTest3()
        {
            Tuple<List<List<Tuple<int, int>>>, List<BigRational>> terms =
                YoungDiagramToSymmPoly.YoungThreeRowsToChoose(1, 1);

            Tuple<List<List<Tuple<int, int>>>, List<BigRational>> expterms =
                SymmPolyToValues.SymmetricInChooseBasis(YoungDiagramToSymmPoly.WedgeToSymmetricPolynomial(2));

            // since terms may not be reduced (i.e., same terms combined)
            // easier to just evaluate each on some cycles and this will ensure they are the same
            foreach (List<int> part in IntegerFunctions.AllPartitions(8))
            {
                List<int> cycles = IntegerFunctions.PartitionToNumCycles(part);

                Assert.IsTrue(RepTheoryAlgs.ChooseFormToCharacter(terms, cycles) 
                    == RepTheoryAlgs.ChooseFormToCharacter(expterms, cycles));
            }
        }

        [TestMethod]
        public void YoungDiagramToChooseTest1()
        {
            Tuple<List<List<Tuple<int, int>>>, List<BigRational>> expterms =
                YoungDiagramToSymmPoly.YoungTwoRowsToChoose(4);

            Tuple<List<List<Tuple<int, int>>>, List<BigRational>> terms =
                YoungDiagramToSymmPoly.YoungDiagramToChoose(new() { 4 });

            // since terms may not be reduced (i.e., same terms combined)
            // easier to just evaluate each on some cycles and this will ensure they are the same
            foreach (List<int> part in IntegerFunctions.AllPartitions(8))
            {
                List<int> cycles = IntegerFunctions.PartitionToNumCycles(part);

                Assert.IsTrue(RepTheoryAlgs.ChooseFormToCharacter(terms, cycles)
                    == RepTheoryAlgs.ChooseFormToCharacter(expterms, cycles));
            }
        }

        [TestMethod]
        public void YoungDiagramToChooseTest2()
        {
            Tuple<List<List<Tuple<int, int>>>, List<BigRational>> expterms =
                YoungDiagramToSymmPoly.YoungThreeRowsToChoose(5, 3);

            Tuple<List<List<Tuple<int, int>>>, List<BigRational>> terms =
                YoungDiagramToSymmPoly.YoungDiagramToChoose(new() { 5, 3 });

            // since terms may not be reduced (i.e., same terms combined)
            // easier to just evaluate each on some cycles and this will ensure they are the same
            foreach (List<int> part in IntegerFunctions.AllPartitions(9))
            {
                List<int> cycles = IntegerFunctions.PartitionToNumCycles(part);

                Assert.IsTrue(RepTheoryAlgs.ChooseFormToCharacter(terms, cycles)
                    == RepTheoryAlgs.ChooseFormToCharacter(expterms, cycles));
            }
        }

        [TestMethod]
        public void YoungDiagramToChooseTest3()
        {
            Tuple<List<List<Tuple<int, int>>>, List<BigRational>> expterms =
                YoungDiagramToSymmPoly.YoungThreeRowsToChoose(4, 4);

            Tuple<List<List<Tuple<int, int>>>, List<BigRational>> terms =
                YoungDiagramToSymmPoly.YoungDiagramToChoose(new() { 4, 4 });

            LaurentPolynomial result = SymmPolyToValues.CyclicPolynomialBasisToPolynomial(terms, -10);
            LaurentPolynomial expresult = SymmPolyToValues.CyclicPolynomialBasisToPolynomial(expterms, -10);

            Assert.IsTrue(result == expresult);
        }

        [TestMethod]
        public void YoungDiagramToChooseTest4()
        {
            LaurentPolynomial expseries =
                SymmPolyToValues.SymmPolyToPoly(YoungDiagramToSymmPoly.WedgeToSymmetricPolynomial(5));

            Tuple<List<List<Tuple<int, int>>>, List<BigRational>> terms =
                YoungDiagramToSymmPoly.YoungDiagramToChoose(new() { 1, 1, 1, 1, 1 });

            LaurentPolynomial result = SymmPolyToValues.CyclicPolynomialBasisToPolynomial(terms, -10);

            Assert.IsTrue(result == expseries);
        }
    }
}
