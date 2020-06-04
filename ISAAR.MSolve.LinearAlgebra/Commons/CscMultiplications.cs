using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.LinearAlgebra.Vectors;

//TODO: many of these may be able to be simplified to avoid much code duplication. An approach would be to have dedicated 1D 
//      vectors that are rows or columns of matrices. Then the CSR loops will operate only on them. The rest of the methods, 
//      will just need to provide the correct rows/columns.
//TODO: These should not be in the Commons namespace.
namespace ISAAR.MSolve.LinearAlgebra.Commons
{
    /// <summary>
    /// Implementations of multiplication operations with a matrix stored in CSC format.
    /// Authors: Ambrosios Savvides for personal use but can be used freely by the team.
    /// </summary>
    internal static class CscMultiplications
    {
        internal static void CscTimesMatrix(int numCscCols, double[] cscValues, int[] cscColOffsets, int[] cscRowIndices,
            IMatrixView other, Matrix result)
        {
            for (int c = 0; c < result.NumRows; ++c) 
            {
                for (int i = 0; i < numCscCols; ++i)
                {
                    var dot = 0.0;
                    for (int l = 0; l < numCscCols; l++)
                    {
                        int cscColStart = cscColOffsets[l]; 
                        int cscColEnd = cscColOffsets[l + 1]; 
                        for (int k = cscColStart; k < cscColEnd; ++k)
                        {
                            if (cscRowIndices[k] == c)
                            {
                                dot += cscValues[k] * other[l, i];
                            }
                        }
                    }
                    result[c, i] = dot;
                }
            }
        }

        internal static void CscTimesMatrixTrans(int numCscCols, double[] cscValues, int[] cscColOffsets, int[] cscRowIndices,
            IMatrixView other, Matrix result)
        {
            for (int c = 0; c < result.NumRows; ++c)
            {
                for (int i = 0; i < numCscCols; ++i)
                {
                    var dot = 0.0;
                    for (int l = 0; l < numCscCols; l++)
                    {
                        int cscColStart = cscColOffsets[l];
                        int cscColEnd = cscColOffsets[l + 1];
                        for (int k = cscColStart; k < cscColEnd; ++k)
                        {
                            if (cscRowIndices[k] == c)
                            {
                                dot += cscValues[k] * other[i, l];
                            }
                        }
                    }
                    result[c, i] = dot;
                }
            }
        }

        internal static void CscTimesVector(int numCscCols, double[] cscValues, int[] cscColOffsets, int[] cscRowIndices,
            IVectorView lhs, double[] rhs)
        {
            for (int i = 0; i < numCscCols; ++i)
            {
                double dot = 0.0;
                for (int j = 0; j < numCscCols; ++j)
                {
                    int rowStart = cscColOffsets[j]; //inclusive
                    int rowEnd = cscColOffsets[j + 1]; //exclusive
                    for (int k = rowStart; k < rowEnd; ++k)
                    {
                        if (cscRowIndices[k] == i)
                        {
                            dot += cscValues[k] * lhs[j];
                        }
                    }
                }
                rhs[i] = dot;
            }
        }

        internal static void CscTimesVector(int numCscCols, double[] cscValues, int[] cscColOffsets, int[] cscRowIndices,
            IVectorView lhs, IVector rhs)
        {
            for (int i = 0; i < numCscCols; ++i)
            {
                double dot = 0.0;
                for (int j = 0; j < numCscCols; ++j)
                {
                    int rowStart = cscColOffsets[j]; //inclusive
                    int rowEnd = cscColOffsets[j + 1]; //exclusive
                    for (int k = rowStart; k < rowEnd; ++k)
                    {
                        if (cscRowIndices[k] == i)
                        {
                            dot += cscValues[k] * lhs[j];
                        }
                    }
                }
                rhs.Set(i, dot);
            }
        }

        internal static void CscTransTimesMatrix(int numCscCols, double[] cscValues, int[] cscColOffsets, int[] cscRowIndices,
            IMatrixView other, Matrix result)
        {
            for (int c = 0; c < result.NumRows; ++c)
            {
                for (int i = 0; i < numCscCols; ++i)
                {
                    int cscColStart = cscColOffsets[c];
                    int cscColEnd = cscColOffsets[c + 1];
                    for (int k = cscColStart; k < cscColEnd; ++k)
                    {
                        result[c, i] += other[cscRowIndices[k], i] * cscValues[k];
                    }
                }
            }
        }

        internal static void CscTransTimesMatrixTrans(int numCscCols, double[] cscValues, int[] cscColOffsets,
            int[] cscRowIndices, IMatrixView other, Matrix result)
        {
            for (int c = 0; c < result.NumRows; ++c)
            {
                for (int i = 0; i < numCscCols; ++i)
                {
                    int cscColStart = cscColOffsets[c];
                    int cscColEnd = cscColOffsets[c + 1];
                    for (int k = cscColStart; k < cscColEnd; ++k)
                    {
                        result[c, i] += other[i,cscRowIndices[k]] * cscValues[k];
                    }
                }
            }
        }

        internal static void CscTransTimesVector(int numCscCols, double[] cscValues, int[] cscColOffsets, int[] cscRowIndices,
            IVectorView lhs, double[] rhs)
        {
            for (int i = 0; i < lhs.Length; ++i)
            {
                double dot = 0.0;

                int rowStart = cscColOffsets[i]; //inclusive
                int rowEnd = cscColOffsets[i + 1]; //exclusive
                for (int k = rowStart; k < rowEnd; ++k)
                {
                    dot += cscValues[k] * lhs[cscRowIndices[k]];
                }

                rhs[i] = dot;
            }
        }

        internal static void CscTransTimesVector(int numCscCols, double[] cscValues, int[] cscColOffsets, int[] cscRowIndices,
            IVectorView lhs, IVector rhs)
        {
            var temp = new double[rhs.Length];
            CscTransTimesVector(numCscCols, cscValues, cscColOffsets, cscRowIndices, lhs, temp);
            rhs.CopyFrom(Vector.CreateFromArray(temp));

            // The following requires a lot of indexing into the rhs vector.
            //// A^T * x = linear combination of columns of A^T = rows of A, with the entries of x as coefficients
            //for (int i = 0; i < numCsrRows; ++i)
            //{
            //    double scalar = lhs[i];
            //    int rowStart = csrRowOffsets[i]; //inclusive
            //    int rowEnd = csrRowOffsets[i + 1]; //exclusive
            //    for (int k = rowStart; k < rowEnd; ++k)
            //    {
            //        rhs.Set(cscRowIndices[k], rhs[cscRowIndices[k]] + scalar * csrValues[k]);
            //    }
            //}
        }

        internal static void MatrixTimesCsc(int numCscCols, double[] cscValues, int[] cscColOffsets, int[] cscRowIndices,
         IMatrixView other, Matrix result)
        {
            for (int r = 0; r < result.NumColumns; ++r)
            {
                for (int i = 0; i < numCscCols; ++i)
                {
                    int cscColStart = cscColOffsets[i];
                    int cscColEnd = cscColOffsets[i + 1];
                    for (int k = cscColStart; k < cscColEnd; ++k)
                    {
                        result[r, i] += other[r, cscRowIndices[k]] * cscValues[k];
                    }
                }
            }
        }

        internal static void MatrixTimesCscTrans(int numCscCols, double[] cscValues, int[] cscColOffsets, int[] cscRowIndices,
         IMatrixView other, Matrix result)
        {
            for (int r = 0; r < result.NumColumns; ++r)
            {
                var rowww = 0;
                for (int i = 0; i < numCscCols; ++i)
                {
                    var dot = 0.0;
                    for (int l = 0; l < numCscCols; l++)
                    {
                        int cscColStart = cscColOffsets[l]; //inclusive
                        int cscColEnd = cscColOffsets[l + 1]; //exclusive
                        for (int k = cscColStart; k < cscColEnd; ++k)
                        {
                            if (cscRowIndices[k] == rowww)
                            {
                                dot += cscValues[k] * other[r, l];
                            }
                        }
                    }
                    rowww++;
                    result[r, i] = dot;
                }
            }
        }

        internal static void MatrixTransTimesCsc(int numCscCols, double[] cscValues, int[] cscColOffsets, int[] cscRowIndices,
            IMatrixView other, Matrix result)
        {
            for (int r = 0; r < result.NumColumns; ++r)
            {
                for (int i = 0; i < numCscCols; ++i)
                {
                    int cscColStart = cscColOffsets[i];
                    int cscColEnd = cscColOffsets[i + 1];
                    for (int k = cscColStart; k < cscColEnd; ++k)
                    {
                        result[r, i] += other[cscRowIndices[k], r] * cscValues[k];
                    }
                }
            }
        }
        internal static void CscDiagonalTimesVector(int numCscCols, double[] cscValues, int[] cscColOffsets, int[] cscRowIndices,
     IVectorView lhs, double[] rhs)
        {
            for (int i = 0; i < rhs.Length; i++)
            {
                int cscColStart = cscColOffsets[i];
                int cscColEnd = cscColOffsets[i + 1];
                for (int k = cscColStart; k < cscColEnd; ++k)
                {
                    if (cscRowIndices[k] == i)
                    {
                        rhs[i] = lhs[i] * cscValues[k];
                    }
                }
            }
        }
        internal static void CscDiagonalTimesVector(int numCscCols, double[] cscValues, int[] cscColOffsets, int[] cscRowIndices,
  IVectorView lhs, IVector rhs)
        {
            var temp = new double[rhs.Length];
            CscDiagonalTimesVector(numCscCols, cscValues, cscColOffsets, cscRowIndices, lhs, temp);
            rhs.CopyFrom(Vector.CreateFromArray(temp));
        }
        internal static void MatrixTransTimesCscTrans(int numCscCols, double[] cscValues, int[] cscColOffsets,
            int[] cscRowIndices, IMatrixView other, Matrix result)
        {
            for (int r = 0; r < result.NumColumns; ++r)
            {
                var rowww = 0;
                for (int i = 0; i < numCscCols; ++i)
                {
                    var dot = 0.0;
                    for (int l = 0; l < numCscCols; l++)
                    {
                        int cscColStart = cscColOffsets[l]; //inclusive
                        int cscColEnd = cscColOffsets[l + 1]; //exclusive
                        for (int k = cscColStart; k < cscColEnd; ++k)
                        {
                            if (cscRowIndices[k] == rowww)
                            {
                                dot += cscValues[k] * other[l, r];
                            }
                        }
                    }
                    rowww++;
                    result[r, i] = dot;
                }
            }
        }

    }
}
