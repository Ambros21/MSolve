using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using ISAAR.MSolve.PreProcessor.Interfaces;

namespace ISAAR.MSolve.PreProcessor.Materials
{
    public class SteelFiberElementMaterial : IFiberFiniteElementMaterial
    {
        private readonly List<IFiberMaterial> fiberMaterials;
        private double youngModulus;
        public double PoissonRatio { get; set; }
        public double HardeningRatio { get; set; }
        public double YieldStress { get; set; }
        public double[] Coordinates { get; set; }

        public SteelFiberElementMaterial(int noOfFibers)
        {
            fiberMaterials = new List<IFiberMaterial>(noOfFibers);
            for (int i = 0; i < noOfFibers; i++) fiberMaterials.Add(new SteelFiberMaterial(this));
        }

        public double YoungModulus 
        {
            get { return youngModulus; }
            set 
            { 
                youngModulus = value;
                foreach (IFiberMaterial material in fiberMaterials)
                {
                    SteelFiberMaterial m = (SteelFiberMaterial)material;
                    m.YoungModulus = value;
                }
            }
        }
        #region IFiberFiniteElementMaterial Members

        public IList<IFiberMaterial> FiberMaterials
        {
            get { return fiberMaterials; }
        }

        #endregion

        #region IFiniteElementMaterial Members

        public int ID
        {
            get { throw new NotImplementedException(); }
        }

        public bool Modified
        {
            get 
            {
                bool modified = false;
                foreach (SteelFiberMaterial material in fiberMaterials)
                    if (material.Modified)
                    {
                        modified = true;
                        break;
                    }
                return modified; 
            }
        }

        public void ResetModified()
        {
            foreach (SteelFiberMaterial material in fiberMaterials) material.ResetModified();
        }

        #endregion

        #region ICloneable Members

        public object Clone()
        {
            SteelFiberElementMaterial m = new SteelFiberElementMaterial(this.fiberMaterials.Count) 
            { 
                HardeningRatio = this.HardeningRatio, 
                PoissonRatio = this.PoissonRatio, 
                YieldStress = this.YieldStress,
                youngModulus = this.youngModulus
            };
            m.fiberMaterials.Clear();
            foreach (IFiberMaterial f in this.fiberMaterials) m.fiberMaterials.Add(f.Clone(m));

            return m;
        }

        #endregion

    }
}
