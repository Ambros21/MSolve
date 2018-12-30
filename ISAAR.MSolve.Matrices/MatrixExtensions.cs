using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using Accord.Math;

namespace ISAAR.MSolve.Matrices
{
    public static class MatrixExtensions
    {
        public static Matrix2D<double> Invert(this Matrix2D<double> matrix)
        {
            return new Matrix2D<double>(matrix.Data.Inverse());
        }

        public static SymmetricMatrix2D<double> Invert(this SymmetricMatrix2D<double> matrix)
        {
            Matrix2D<double> originalFull = matrix.ToMatrix2D();
            Matrix2D<double> inverseFull = originalFull.Invert();
            SymmetricMatrix2D<double> inverseSymmetric = new SymmetricMatrix2D<double>(inverseFull);
            return inverseSymmetric;
        }
    }
}
