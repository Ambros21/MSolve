using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using ISAAR.MSolve.PreProcessor.Interfaces;

namespace ISAAR.MSolve.PreProcessor.Materials
{
    public class SteelFiberMaterial : IFiberMaterial
    {
        private  SteelFiberElementMaterial elementMaterial;
        private double youngModulus, strain, stress, initialStrain, initialStress, fsyt, fsyc, sy,fsytb,fsycb;
        private bool modified = false;
        private bool initialValuesInitialized = false;

        public SteelFiberMaterial(SteelFiberElementMaterial elementMaterial,double YoungModulus,double PoissonRatio,double HardeningRatio,double YieldStressInitial,double YieldStressTension, double YieldStressCompression)
        {
            this.elementMaterial = elementMaterial;
            this.youngModulus = YoungModulus;
            elementMaterial.YoungModulus = this.youngModulus;
            this.strain = 0.0;
            this.stress = 0.0;
            this.initialStress = 0.0;
            this.initialStrain = 0.0;
            elementMaterial.PoissonRatio = PoissonRatio;
            elementMaterial.HardeningRatio = HardeningRatio;
            elementMaterial.YieldStressInitial = YieldStressInitial;
            elementMaterial.YieldStressTension = YieldStressTension;
            elementMaterial.YieldStressCompression = YieldStressCompression;
            this.fsyc = elementMaterial.YieldStressCompression;
            this.fsyt = elementMaterial.YieldStressTension;
            this.fsycb = elementMaterial.YieldStressCompression;
            this.fsytb = elementMaterial.YieldStressTension;
            this.sy = elementMaterial.YieldStressInitial;
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
            sy = elementMaterial.YieldStressInitial;
            fsyt = elementMaterial.YieldStressTension;
            fsyc = elementMaterial.YieldStressCompression; // Need to check if we get the correct balanced situation in the begining and we save it in savestate.
            double e0 = elementMaterial.YoungModulus;
            double b = elementMaterial.HardeningRatio;
            double ds = e0 * dStrain;
            strain = initialStrain + dStrain;
            double strial = initialStress + ds;
            if ((strial>fsytb)||(strial<fsycb))
            {
                if(strial>0)
                {
                    var help = (fsytb-initialStress) / e0;
                    stress = fsytb + b * e0 * (dStrain - help);
                    fsyt = stress;
                    fsyc = -(2 * sy - fsyt);
                }
                else
                {
                    var help = (fsycb - initialStress) / e0;
                    stress = fsycb + b * e0 * (dStrain - help);
                    fsyc = stress;
                    fsyt = Math.Abs((2 * sy - Math.Abs(fsyc)));
                }
                youngModulus = b * e0;
                Modified = true;
            }
            else
            {
                stress = strial;
                youngModulus = e0;
            }
        }

        public void SaveState()
        {
            this.fsycb = this.fsyc;
            this.fsytb = this.fsyt;
            this.initialStrain = this.strain; //this variable is used here as the balanced strain at the last equillibrium.
            this.initialStress = this.stress; //this variable is used here as the balanced stress at the last equillibrium.
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
            SteelFiberElementMaterial h = (SteelFiberElementMaterial)parent;
            SteelFiberMaterial m = new SteelFiberMaterial(h,h.YoungModulus,h.PoissonRatio,h.HardeningRatio,h.YieldStressInitial,h.YieldStressTension,h.YieldStressCompression);
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
