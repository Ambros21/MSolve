using System;
using System.Collections.Generic;
using System.Text;

namespace ISAAR.MSolve.FEM.Entities.TemporalFunctions
{
    public interface ITemporalFunction
    {
        double CalculateValueAt(int timeStep);
    }
}
