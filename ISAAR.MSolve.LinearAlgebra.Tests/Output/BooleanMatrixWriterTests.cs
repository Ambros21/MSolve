﻿using System.IO;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.LinearAlgebra.Output;
using ISAAR.MSolve.LinearAlgebra.Tests.Matrices;
using ISAAR.MSolve.LinearAlgebra.Tests.TestData;
using ISAAR.MSolve.LinearAlgebra.Tests.Utilities;
using Xunit;

namespace ISAAR.MSolve.LinearAlgebra.Tests.Output
{
    /// <summary>
    /// Tests for <see cref="BooleanMatrixWriter"/>.
    /// Authors: Serafeim Bakalakos
    /// </summary>
    public static class BooleanMatrixWriterTests
    {
        [Fact]
        private static void TestWriter1()
        {
            var matrix = SignedBooleanMatrixTests.CreateMatrix(SignedBoolean5by10.A1);
            TestWriteOperation(matrix, SignedBoolean5by10.filePath1);
        }

        [Fact]
        private static void TestWriter2()
        {
            var matrix = SignedBooleanMatrixTests.CreateMatrix(SignedBoolean5by10.A2);
            TestWriteOperation(matrix, SignedBoolean5by10.filePath2);
        }

        private static void TestWriteOperation(SignedBooleanMatrix matrix, string referenceFile)
            => TestWriteOperation(matrix, referenceFile, new BooleanMatrixWriter(true));

        private static void TestWriteOperation(SignedBooleanMatrix matrix, string referenceFile, BooleanMatrixWriter writer)
        {
            string tempFile = "temp.txt";
            writer.WriteToFile(matrix, tempFile);
            bool success = IOUtilities.AreFilesEquivalent(referenceFile, tempFile);
            File.Delete(tempFile);
            Assert.True(success);
        }
    }
}
