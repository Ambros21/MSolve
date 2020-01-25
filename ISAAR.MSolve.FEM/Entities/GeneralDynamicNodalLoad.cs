using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.FEM.Entities.TemporalFunctions;

namespace ISAAR.MSolve.FEM.Entities
{
    public class GeneralDynamicNodalLoad : ITimeDependentNodalLoad
    {
        private readonly ITemporalFunction temporalFunction;

        public GeneralDynamicNodalLoad(Node node, DOFType dof, ITemporalFunction temporalFunction)
        {
            this.Node = node;
            this.DOF = dof;
            this.temporalFunction = temporalFunction;
        }

        public Node Node { get; }

        public DOFType DOF { get; }

        public double GetLoadAmount(int timeStep) => temporalFunction.CalculateValueAt(timeStep);
    }
}
