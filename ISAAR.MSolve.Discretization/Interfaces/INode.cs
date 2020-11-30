using System;
using System.Collections.Generic;

namespace ISAAR.MSolve.Discretization.Interfaces
{
    public interface INode : IComparable<INode>
    {
		int ID { get; }
		double X { get; }
		double Y { get; }
		double Z { get; }

        List<Constraint> Constraints { get; }
        List<InitialCondition> InitialConditions { get; }
        Dictionary<int, ISubdomain> SubdomainsDictionary { get; }
    }
}
