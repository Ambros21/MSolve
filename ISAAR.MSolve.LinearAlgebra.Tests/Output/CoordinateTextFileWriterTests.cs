﻿using System.IO;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.LinearAlgebra.Output;
using ISAAR.MSolve.LinearAlgebra.Tests.TestData;
using ISAAR.MSolve.LinearAlgebra.Tests.Utilities;
using Xunit;

namespace ISAAR.MSolve.LinearAlgebra.Tests.Output
{
    /// <summary>
    /// Tests for <see cref="CoordinateTextFileWriter"/>.
    /// Authors: Serafeim Bakalakos
    /// </summary>
    public static class CoordinateTextFileWriterTests
    {
        [Fact]
        private static void TestSparsePosDef()
        {
            var matrix = SkylineMatrix.CreateFromArrays(SparsePosDef10by10.order, SparsePosDef10by10.skylineValues,
                SparsePosDef10by10.skylineDiagOffsets, true, true);
            TestWriteOperation(matrix, SparsePosDef10by10.coordinateFormatPath);
        }

        [Fact]
        private static void TestSparseRectangularCSC()
        {
            // The reference file has been written by iterating the entries in CSC order.
            // Therefore it cannot be used for checking with a CSR matrix.
            var matrix = CscMatrix.CreateFromArrays(SparseRectangular10by5.numRows, SparseRectangular10by5.numCols,
                SparseRectangular10by5.cscValues, SparseRectangular10by5.cscRowIndices, SparseRectangular10by5.cscColOffsets,
                true);
            TestWriteOperation(matrix, SparseRectangular10by5.coordinateFormatPath);
        }

        private static void TestWriteOperation(ISparseMatrix matrix, string referenceFile)
            => TestWriteOperation(matrix, referenceFile, new CoordinateTextFileWriter());

        private static void TestWriteOperation(ISparseMatrix matrix, string referenceFile, CoordinateTextFileWriter writer)
        {
            string tempFile = "temp.txt";
            writer.WriteToFile(matrix, tempFile);
            bool success = IOUtilities.AreFilesEquivalent(referenceFile, tempFile);
            File.Delete(tempFile);
            Assert.True(success);
        }
    }
}
