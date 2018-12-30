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

            public SymmetricMatrix2D<double> sectionstiffnessMatrix { get; set; }

            public SymmetricMatrix2D<double> sectionflexibilityMatrix { get; set; }

            public Vector<double> sectionforces { get; set; }

            public Vector<double> sectionforcesbalanced { get; set; }

            public double[] residualdef { get; set; }

            public Vector<double> deformations { get; set; }

            public Vector<double> deformationsbalanced { get; set; }

            #endregion
        }
        // it has only the n=4 case with the respective points xi and the weight integral factors
        #region Public Methods

        public static GaussLobattoPoint1D[] GetGaussLobattoPoints()
        {
            var d = new GaussLobattoPoint1D[4];
            d[0]= new GaussLobattoPoint1D
            {
                Coordinate = -1.0,
                WeightFactor = 1.0 / 6.0,
                sectionstiffnessMatrix = new SymmetricMatrix2D<double>(3),
                sectionflexibilityMatrix = new SymmetricMatrix2D<double>(3),
                sectionforces = new Vector<double>(3),
                sectionforcesbalanced = new Vector<double>(3),
                residualdef = new double[3],
                deformations = new Vector<double>(3),
                deformationsbalanced = new Vector<double>(3)
            };
            d[1]= new GaussLobattoPoint1D
            {
                Coordinate = -0.44721359549,
                WeightFactor = 5.0 / 6.0,
                sectionstiffnessMatrix = new SymmetricMatrix2D<double>(3),
                sectionflexibilityMatrix = new SymmetricMatrix2D<double>(3),
                sectionforces = new Vector<double>(3),
                sectionforcesbalanced = new Vector<double>(3),
                residualdef = new double[3],
                deformations = new Vector<double>(3),
                deformationsbalanced = new Vector<double>(3)
            };
            d[2] = new GaussLobattoPoint1D
            {
                Coordinate = 0.44721359549,
                WeightFactor = 5.0 / 6.0,
                sectionstiffnessMatrix = new SymmetricMatrix2D<double>(3),
                sectionflexibilityMatrix = new SymmetricMatrix2D<double>(3),
                sectionforces = new Vector<double>(3),
                sectionforcesbalanced = new Vector<double>(3),
                residualdef = new double[3],
                deformations = new Vector<double>(3),
                deformationsbalanced = new Vector<double>(3)
            };
            d[3] = new GaussLobattoPoint1D
            {
                Coordinate = 1.0,
                WeightFactor = 1.0 / 6.0,
                sectionstiffnessMatrix = new SymmetricMatrix2D<double>(3),
                sectionflexibilityMatrix = new SymmetricMatrix2D<double>(3),
                sectionforces = new Vector<double>(3),
                sectionforcesbalanced = new Vector<double>(3),
                residualdef = new double[3],
                deformations = new Vector<double>(3),
                deformationsbalanced = new Vector<double>(3)
            };
            return d;
        }

        #endregion
    }
}

