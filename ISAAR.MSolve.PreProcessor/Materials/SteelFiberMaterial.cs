using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using ISAAR.MSolve.PreProcessor.Interfaces;

namespace ISAAR.MSolve.PreProcessor.Materials
{
    public class SteelFiberMaterial : IFiberMaterial
    {
        private readonly SteelFiberElementMaterial elementMaterial;
        private double youngModulus, strain, stress, initialStrain, initialStress;
        private bool modified = false;
        private bool initialValuesInitialized = false;

        public SteelFiberMaterial(SteelFiberElementMaterial elementMaterial)
        {
            this.elementMaterial = elementMaterial;
            YoungModulus = elementMaterial.YoungModulus;
        }

        public SteelFiberMaterial(SteelFiberElementMaterial elementMaterial, SteelFiberMaterial sourceMaterial)
        {
            this.elementMaterial = elementMaterial;
            youngModulus = sourceMaterial.YoungModulus;
            strain = sourceMaterial.strain;
            stress = sourceMaterial.stress;
            initialStrain = sourceMaterial.initialStrain;
            initialStress = sourceMaterial.initialStress;
            modified = sourceMaterial.modified;
            initialValuesInitialized = sourceMaterial.initialValuesInitialized;
        }

        #region IFiberMaterial Members

        public double YoungModulus 
        {
            get { return youngModulus; }
            set { youngModulus = value; }
        }

        public double YoungModulusElastic
        {
            get
            {
                return elementMaterial.YoungModulus;
            }
            set
            {
                throw new InvalidOperationException("Cannot change elastic Young modulus of fiber.");
            }
        }

        public double PoissonRatio
        {
            get
            {
                return elementMaterial.PoissonRatio;
            }
            set
            {
                throw new InvalidOperationException("Cannot change elastic Poisson ratio of fiber.");
            }
        }

        public double Strain { get { return strain; } }
        public double Stress { get { return stress; } }

        public void UpdateMaterial(double dStrain)
        {
            double fsy = elementMaterial.YieldStress;
            double e0 = elementMaterial.YoungModulus;
            double b = elementMaterial.HardeningRatio;
            double hc = 0;
            double ht = 0;

            double esh = b * e0;
            double epsy = fsy / e0;
            double epss = initialStrain + strain + dStrain;
            //double epss = initialStrain + dStrain;

            double c1 = esh * epss;
            double c2 = (fsy + hc) * (1 - b);
            double c3 = (fsy + ht) * (1 - b);
            double c = initialStress + stress + e0 * dStrain;
            //double c = initialStress + e0 * dStrain;
            double sigs = Math.Max(c1 - c2, Math.Min(c1 + c3, c));

            double oldYoungModulus = youngModulus;
            youngModulus = esh;
            if (Math.Abs(sigs - c) < 1e-10) youngModulus = e0;
            if (oldYoungModulus != youngModulus) Modified = true;
            strain = epss - initialStrain;
            stress = sigs;

            //if (!initialValuesInitialized)
            //{
            //    initialValuesInitialized = true;
            //    initialStrain = strain;
            //    initialStress = stress;
            //}
        }

        public void SaveState()
        {
            initialStrain = strain;
            initialStress = stress;
            stress = 0;
            strain = 0;
        }

        public void ClearStresses()
        {
            initialStress = 0;
            stress = 0;
        }

        #endregion

        #region IFiniteElementMaterial Members

        public int ID
        {
            get { throw new NotImplementedException(); }
        }

        public bool Modified
        {
            get { return modified; }
            set { modified = value; }
        }

        public void ResetModified()
        {
            modified = false;
        }

        public IFiberMaterial Clone(IFiberFiniteElementMaterial parent)
        {
            SteelFiberMaterial m = new SteelFiberMaterial((SteelFiberElementMaterial)parent);
            m.youngModulus = this.youngModulus;
            m.initialStrain = this.initialStrain;
            m.initialStress = this.initialStress;
            m.strain = this.strain;
            m.stress = this.stress;
            m.modified = this.modified;
            m.initialValuesInitialized = this.initialValuesInitialized;

            return m;
        }

        public double[] Coordinates { get; set; }
        #endregion


        #region ICloneable Members

        public object Clone()
        {
            throw new InvalidOperationException("Cannot clone SteelFiberMaterial because parent information is always missing.");
        }

        #endregion
    }
}
