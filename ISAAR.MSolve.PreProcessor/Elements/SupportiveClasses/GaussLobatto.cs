using ISAAR.MSolve.Matrices;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace ISAAR.MSolve.PreProcessor.Elements.SupportiveClasses
{
    public class GaussLobatto
    {
        public class GaussLobattoPoint1D
        {
            #region Properties

            public double Coordinate { get; set; }

            public double WeightFactor { get; set; }

            public SymmetricMatrix2D<double> sectionstiffnessMatrix = new SymmetricMatrix2D<double>(6);

            public SymmetricMatrix2D<double> sectionflexibilityMatrix = new SymmetricMatrix2D<double>(6);

            public double[] sectionforces = new double[6];

            public double[] sectionforcesprev = new double[6];

            public double[] residualdef = new double[6];

            public double[] deformations = new double[6];

            #endregion
        }
        // It has only the n=4 case with the respective points xi and the weight integral factors
        private static readonly GaussLobattoPoint1D GaussLobaattoPoint1 = new GaussLobattoPoint1D
        {
            Coordinate = 1,
            WeightFactor = 1 / 6
        };
        private static readonly GaussLobattoPoint1D GaussLobaattoPoint2 = new GaussLobattoPoint1D
        {
            Coordinate = -1,
            WeightFactor = 1 / 6
        };
        private static readonly GaussLobattoPoint1D GaussLobaattoPoint3 = new GaussLobattoPoint1D
        {
            Coordinate = 0.44721359549,
            WeightFactor = 5 / 6
        };
        private static readonly GaussLobattoPoint1D GaussLobaattoPoint4 = new GaussLobattoPoint1D
        {
            Coordinate = -0.44721359549,
            WeightFactor = 5 / 6
        };
        #region Public Methods

        public static GaussLobattoPoint1D[] GetGaussLobattoPoints()
        {
            return new[] { GaussLobaattoPoint1, GaussLobaattoPoint2, GaussLobaattoPoint3, GaussLobaattoPoint4 };
        }

        #endregion
    }
}

