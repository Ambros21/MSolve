using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.Discretization.FreedomDegrees;
using ISAAR.MSolve.FEM.Entities.TemporalFunctions;

namespace ISAAR.MSolve.FEM.Entities
{
    public class GeneralDynamicNodalLoad : ITimeDependentNodalLoad
    {
        private readonly ITemporalFunction temporalFunction;

        public GeneralDynamicNodalLoad(Node node,IDofType  dof, ITemporalFunction temporalFunction)
        {
            this.Node = node;
            this.DOF = dof;
            this.temporalFunction = temporalFunction;
        }

        public Node Node { get; set; }

        public IDofType DOF { get; set; }

        public double GetLoadAmount(int timeStep) => temporalFunction.CalculateValueAt(timeStep);
    }
}
