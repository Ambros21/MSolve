using ISAAR.MSolve.Discretization.Interfaces;
using System;
using System.Collections.Generic;
using System.Text;

namespace ISAAR.MSolve.FEM.Entities
{
    public interface ITimeDependentNodalLoad
    {
        Node Node { get; }
        DOFType DOF { get; }

        double GetLoadAmount(int timeStep);
    }
}
