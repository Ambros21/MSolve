﻿using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace ISAAR.MSolve.LinearAlgebra.Testing
{
    // TODO: need a symmetric, invertible, indefinite benchmark.
    public static class LinearAlgebraTestSuite
    {
        public static void TestAll()
        {
            //SuiteSparse();
            //TestMarshaling();
            //TestBuilders();
            //TestInverse();
            //TestIndexing();
            //TestEquality();
            //TestFactorization();
            //TestMatrixVectorMultiplication();
            //TestMatrixMatrixMultiplication();
            //TestReading();
            //TestTransposition();
            //TestVectorOperations();
            //TestWriting();
        }

        public static void SuiteSparse()
        {
            //TestSuiteSparse.ExampleMatrixClasses();
            //TestSuiteSparse.CheckRowAddition();
            //TestSuiteSparse.CheckRowAdditionReverse();
            //TestSuiteSparse.CheckRowDeletion();
            //TestSuiteSparse.CheckReordering1();
            //TestSuiteSparse.CheckSystemSolutions();
        }

        public static void TestReading()
        {
            //TestInput.CheckArray1DReader();
            //TestInput.CheckArray2DReader();
        }

        public static void TestWriting()
        {
            //DenseVectors.Print();
            //Console.WriteLine();
            //RectangularFullColRank.Print();
            //Console.WriteLine();
            //SquareInvertible.Print();
            //Console.WriteLine();
            //SquareSingular.Print();
            //Console.WriteLine();
            //SquareSingular1Deficiency.Print();
            //Console.WriteLine();
            //LowerInvertible.Print();
            //Console.WriteLine();
            //LowerSingular.Print();
            //Console.WriteLine();
            //UpperInvertible.Print();
            //Console.WriteLine();
            //UpperSingular.Print();
            //Console.WriteLine();
            //SymmPositiveDefinite.Print();
            //Console.WriteLine();
            //SymmSingular.Print();
            //Console.WriteLine();
            //SparseRect.Print();
            //SparsePositiveDefinite.Print();
            //Console.WriteLine();
        }

        public static void TestMatrixVectorMultiplication()
        {
            //Rectangular.CheckMatrixVectorMult();
            //SquareInvertible.CheckMatrixVectorMult();
            //SquareSingular.CheckMatrixVectorMult();
            //SquareSingular1Deficiency.CheckMatrixVectorMult();
            //SymmPositiveDefinite.CheckMatrixVectorMult();
            //SymmSingular.CheckMatrixVectorMult();
            //SparseRect.CheckMatrixVectorMult();
            //SparsePositiveDefinite.CheckMatrixVectorMult();
            //SignedBoolean.CheckMatrixVectorMultiplication();
        }

    }
}
