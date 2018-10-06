using ISAAR.MSolve.Logging.Interfaces;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace ISAAR.MSolve.Analyzers.Interfaces
{
    public interface IAnalyzer
    {
        IAnalyzer ParentAnalyzer { get; set; }
        IAnalyzer ChildAnalyzer { get; set; }
        Dictionary<int, IAnalyzerLog[]> Logs { get; }
        void BuildMatrices();
        void Initialize();
        void Solve();
        void Solve(int i, int j); //this was in order to help us make a proper dynamic solver in newton raphson. all classes that uses this interfaces are also updated.
    }
}
