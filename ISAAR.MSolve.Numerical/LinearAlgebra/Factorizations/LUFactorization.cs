﻿using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using IntelMKL.LP64;
using ISAAR.MSolve.Numerical.Exceptions;
using ISAAR.MSolve.Numerical.LinearAlgebra.Commons;
using ISAAR.MSolve.Numerical.LinearAlgebra.Matrices;
using ISAAR.MSolve.Numerical.LinearAlgebra.Vectors;
using ISAAR.MSolve.Numerical.MKL;

namespace ISAAR.MSolve.Numerical.LinearAlgebra.Factorizations
{
    /// <summary>
    /// The LU factorization of a matrix A consists of a lower triangular matrix L (with 1 in its diagonal entries), an upper 
    /// triangular matrix U and a permutation matrix P, such that A = P*L*U. This class stores L,U,P in an efficient manner and
    /// provides common methods to use them. 
    /// </summary>
    public class LUFactorization: IFactorization
    {
        // Perhaps a smaller tolerance is appropriate, since the "almost zero" will propagate during back & forward substitution.
        private const double PivotTolerance = 1e-13;

        private readonly double[] lowerUpper;
        private readonly int[] permutation;
        private readonly int firstZeroPivot;

        private LUFactorization(int order, double[] lowerUpper, int[] permutation, int firstZeroPivot, bool isSingular)
        {
            this.Order = order;
            this.lowerUpper = lowerUpper;
            this.permutation = permutation;
            this.firstZeroPivot = firstZeroPivot;
            this.IsSingular = isSingular;
            this.IsOverwritten = false;
        }

        public bool IsOverwritten { get; private set; }

        public bool IsSingular { get; }

        /// <summary>
        /// The number of rows or columns of the matrix. 
        /// </summary>
        public int Order { get; }

        /// <summary>
        /// Calculates the LUP factorization of a square matrix, such that A = P * L * U. Requires an extra O(n^2 + n) 
        /// available memory.
        /// </summary>
        /// <param name="order">The number of rows/columns of the square matrix.</param>
        /// <param name="matrix">The internal buffer stroring the matrix entries in column major order. It will 
        ///     be overwritten.</param>
        /// <param name="pivotTolerance"></param>
        /// <returns></returns>
        public static LUFactorization Factorize(int order, double[] matrix,
            double pivotTolerance = LUFactorization.PivotTolerance)
        {
            // Call MKL
            int[] permutation = new int[order];
            int info = MKLUtilities.DefaultInfo;
            Lapack.Dgetrf(ref order, ref order, ref matrix[0], ref order, ref permutation[0], ref info);

            // Check MKL execution
            int firstZeroPivot = int.MinValue;
            if (info == 0)
            {
                if (Math.Abs(matrix[order * order - 1]) <= pivotTolerance)
                {
                    // False Negative: info = 0, but LAPACK doesn't check the last diagonal entry!
                    firstZeroPivot = order - 1;
                }
                return new LUFactorization(order, matrix, permutation, firstZeroPivot, (firstZeroPivot >= 0));
            }
            else if (info > 0)
            {
                firstZeroPivot = info - 1;
                return new LUFactorization(order, matrix, permutation, firstZeroPivot, true);
            }
            else if (info == MKLUtilities.DefaultInfo)
            {
                // first check the default info value, since it is negative.
                // info == default => the MKL call did not succeed. 
                // info > 0 should not be returned at all by MKL, but it is here for completion.
                throw new MKLException("Something went wrong with the MKL call."
                    + " Please contact the developer responsible for the linear algebra project.");
            }
            else // (info < 0)
            {
                string msg = string.Format("The {0}th parameter has an illegal value.", -info)
                    + " Please contact the developer responsible for the linear algebra project.";
                throw new MKLException(msg);
            }            
        }

        /// <summary>
        /// Calculates the determinant of the original matrix. det(A) = det(L*U) = det(L)*det(U). 
        /// Since all these are triangular matrices their determinants is the product of their diagonal entries:
        /// det(L) = 1*1*...*1 = 1. Thus det(A) = det(U) = U1*U2*...*Un
        /// </summary>
        /// <returns>The determinant of the original matrix.</returns>
        public double CalcDeterminant()
        {
            if (IsOverwritten) throw new InvalidOperationException(
                "The internal buffer of this factorization has been overwritten and thus cannot be used anymore.");

            if (IsSingular) return 0.0;
            else
            {
                double det = 1.0;
                for (int i = 0; i < Order; ++i)
                {
                    det *= lowerUpper[i * Order + i];
                }
                return det;
            }
        }

        // Explicitly composes and returns the L and U matrices. TODO: 1) Use matrix classes instead of arrays, 2) Also return P
        public (Matrix lower, Matrix upper) Expand()
        {
            if (IsOverwritten) throw new InvalidOperationException(
                "The internal buffer of this factorization has been overwritten and thus cannot be used anymore.");

            double[] l = Conversions.FullColMajorToFullLowerColMajor(lowerUpper, true);
            double[] u = Conversions.FullColMajorToFullUpperColMajor(lowerUpper, false);
            return (Matrix.CreateFromArray(l, Order, Order, false), Matrix.CreateFromArray(u, Order, Order, false));
        }

        /// <summary>
        /// Inverts the original square matrix. The matrix must be non singular, otherwise an
        /// <see cref="SingularMatrixException"/> will be thrown. If <paramref name="inPlace"/> is set to true, this object must 
        /// not be used again, otherwise a <see cref="InvalidOperationException"/> will be thrown.
        /// </summary>
        /// <param name="inPlace">False, to copy the internal factorization data before inversion. True, to overwrite it with
        ///     the inverse matrix, thus saving memory and time. However, that will make this object unusable, so you MUST NOT 
        ///     call any other members afterwards.</param>
        /// <returns></returns>
        public Matrix Invert(bool inPlace)
        {
            // Check if the matrix is suitable for inversion
            if (IsOverwritten) throw new InvalidOperationException(
                "The internal buffer of this factorization has been overwritten and thus cannot be used anymore.");
            if (IsSingular) throw new SingularMatrixException("The factorization has been completed, but U is singular."
                + $" The first zero pivot is U[{firstZeroPivot}, {firstZeroPivot}] = 0.");

            // Call MKL
            int info = MKLUtilities.DefaultInfo;
            double[] inverse;
            if (inPlace)
            {
                inverse = lowerUpper;
                IsOverwritten = true;
            }
            else
            {
                inverse = new double[lowerUpper.Length];
                Array.Copy(lowerUpper, inverse, lowerUpper.Length);
            }
            info = LAPACKE.Dgetri(LAPACKE.LAPACK_COL_MAJOR, Order, inverse, Order, permutation);

            // Check MKL execution
            if (info == 0) return Matrix.CreateFromArray(inverse, Order, Order, false);
            else if ((info == MKLUtilities.DefaultInfo) || (info > 0))
            {
                // first check the default info value, since it is negative.
                // info == dafeult => the MKL call did not succeed. 
                // info > 0 should not be returned at all by MKL, but it is here for completion.
                throw new MKLException("Something went wrong with the MKL call."
                    + " Please contact the developer responsible for the linear algebra project.");
            }
            else if (info < 0)
            {
                throw new MKLException($"The {-info}th parameter has an illegal value."
                    + " Please contact the developer responsible for the linear algebra project.");
            }
            else  // (info > 0) This should not happen
            {
                throw new SingularMatrixException($"The {info - 1} diagonal element of factor U is zero, U is singular and the"
                    + "inversion could not be completed.");
            }
        }

        public VectorMKL SolveLinearSystem(VectorMKL rhs)
        {
            if (IsOverwritten) throw new InvalidOperationException(
                "The internal buffer of this factorization has been overwritten and thus cannot be used anymore.");

            Preconditions.CheckSystemSolutionDimensions(this.Order, this.Order, rhs.Length);
            // Check if the matrix is singular first
            if (IsSingular)
            {
                string msg = "The factorization has been completed, but U is singular."
                    + $" The first zero pivot is U[{firstZeroPivot}, {firstZeroPivot}] = 0.";
                throw new SingularMatrixException(msg);
            }

            // Back & forward substitution using MKL
            int n = Order;
            double[] b = rhs.CopyToArray();
            int info = MKLUtilities.DefaultInfo;
            int nRhs = 1; // rhs is a n x nRhs matrix, stored in b
            int ldb = n; // column major ordering: leading dimension of b is n 
            Lapack.Dgetrs("N", ref n, ref nRhs, ref lowerUpper[0], ref n, ref permutation[0], ref b[0], ref ldb, ref info);

            // Check MKL execution
            if (info == 0) return VectorMKL.CreateFromArray(b, false);
            else if ((info == MKLUtilities.DefaultInfo) || (info > 0)) 
            {
                // first check the default info value, since it is negative.
                // info == dafeult => the MKL call did not succeed. 
                // info > 0 should not be returned at all by MKL, but it is here for completion.
                throw new MKLException("Something went wrong with the MKL call."
                    + " Please contact the developer responsible for the linear algebra project.");
            }
            else // (info < 0)
            {
                string msg = $"The {-info}th parameter has an illegal value." 
                    + " Please contact the developer responsible for the linear algebra project.";
                throw new MKLException(msg);
            }
        } 
    }
}
