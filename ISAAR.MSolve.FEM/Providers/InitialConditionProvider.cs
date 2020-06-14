using System.Collections.Generic;
using System.Linq;
using ISAAR.MSolve.Discretization;
using ISAAR.MSolve.Discretization.FreedomDegrees;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.FEM.Entities;
using ISAAR.MSolve.FEM.Interfaces;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.LinearAlgebra.Vectors;

namespace ISAAR.MSolve.FEM.Providers
{
    public class InitialConditionProvider : IElementMatrixProvider
    {
        public IVector InitialMatrix(IElement element)
        {
            IPorousFiniteElement elementType = (IPorousFiniteElement)element.ElementType;
            int dofs = 0;
            foreach (IList<IDofType> dofTypes in elementType.DofEnumerator.GetDofTypesForMatrixAssembly(element))
                foreach (IDofType dofType in dofTypes) dofs++;
            Vector Initial = Vector.CreateZero(dofs);
            int row = 0;
            int rowin = 0;
            for (int ii=0;ii<element.Nodes.Count;ii++)
            {
                foreach (IList<IDofType> dofTypesRow in elementType.DofEnumerator.GetDofTypesForMatrixAssembly(element))
                    foreach (IDofType dofTypeRow in dofTypesRow)
                    {
                        if (element.Nodes[ii].InitialConditions.Count > 0)
                        {
                            for (int iii = 0; iii < element.Nodes[ii].InitialConditions.Count; iii++)
                            {
                                Initial[row] = element.Nodes[ii].InitialConditions[rowin].Amount;
                                rowin++;
                            }
                            rowin = 0;
                        }
                        row++;
                    }
            }

            return Initial;
        }

        public IMatrix Matrix(IElement element)
        {
            throw new System.NotImplementedException();
        }
    }
}
