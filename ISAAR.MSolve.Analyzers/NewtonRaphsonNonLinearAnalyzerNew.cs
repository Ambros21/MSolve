using ISAAR.MSolve.Analyzers.Interfaces;
using ISAAR.MSolve.Matrices;
using ISAAR.MSolve.Solvers.Interfaces;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using ISAAR.MSolve.Logging.Interfaces;
using ISAAR.MSolve.Matrices.Interfaces;
using ISAAR.MSolve.Logging;
using System.Diagnostics;

namespace ISAAR.MSolve.Analyzers
{
    /// <summary>
    /// Class for the Newton Rapshon Analyzer.
    /// </summary>
    /// <author>Theofilos Manitaras</author>
    public class NonLinearAnalyzerNewtonRaphsonNew:IAnalyzer
    {
        /// <summary>
        /// The array of the displacements we wish to monitor. It is initialized with size the number of increments. 
        /// </summary>
        public double[] displacements;
        /// <summary>
        /// An integer that indicates if we have failed or not. 0=Not Failure 1=Failure 
        /// </summary>
        public int fail = 0;
        /// <summary>
        /// An integer that indicates the step we have failed. If zero then no failure. 
        /// </summary>
        public int failinc = 0;
        /// <summary>
        /// The number of dof we wish to monitor. It is decleared in Program.cs
        /// </summary>
        public int dofid;
        /// <summary>
        /// The maximum number of iterations allowed in an incremental load step.
        /// </summary>
        private static readonly int MAX_ITERATIOMS = 1000;

        /// <summary>
        /// The maximum residual norm allowed, indicating that the algorithm diverges.
        /// </summary>
        private static readonly double MAX_RESIDUAL_NORM_ALLOWED = 2.0e0;

        /// <summary>
        /// Iterations after which the coefficient matrix is rebuilt.
        /// </summary>
        public int StepForMatrixRebuild { get; set; } = 1;

        /// <summary>
        /// The error tolerance checked for convergence.
        /// </summary>
        public double Tolerance { get; set; } = 1.0e-5;

        /// <summary>
        /// The parent analyzer.
        /// </summary>
        private INonLinearParentAnalyzer parentAnalyzer;

        /// <summary>
        /// The dictionary of displacements of the last equilibrium load.
        /// </summary>
        private Dictionary<int, Vector<double>> displacementMap = new Dictionary<int, Vector<double>>();

        /// <summary>
        /// The global right hand side.
        /// </summary>
        private Vector<double> globalRightHandSide;

        /// <summary>
        /// The dictionary of incremental displacements. The incremental displacements express the sum of all iterative
        /// displacements for the load step.
        /// </summary>
        private Dictionary<int, Vector<double>> incrementalDisplacementMap = new Dictionary<int, Vector<double>>();

        /// <summary>
        /// The incremental load factor of the load step.
        /// </summary>
        private double incrementalLoadFactor;

        /// <summary>
        /// The number of incremental load steps.
        /// </summary>
        private int increments;

        /// <summary>
        /// The dictionary of iterative displacements. The iterative displacements correspond to the last solution of the
        /// linear system.
        /// </summary>
        private Dictionary<int, Vector<double>> iterativeDisplacementMap = new Dictionary<int, Vector<double>>();

        /// <summary>
        /// The total load factor so far.
        /// </summary>
        private double loadFactor;

        /// <summary>
        /// The list of incremental load factors.
        /// </summary>
        private IList<double> loadIncrementsList;

        /// <summary>
        /// The nonlinear provider.
        /// </summary>
        private INonLinearProvider provider;

        /// <summary>
        /// The right hand side L2 norm.
        /// </summary>
        private double rightHandSideNorm;

        /// <summary>
        /// The dictionary of right hand sides.
        /// </summary>
        private Dictionary<int, Vector<double>> rightHandSides = new Dictionary<int, Vector<double>>();

        /// <summary>
        /// The dictionary of the linear analyzer.
        /// </summary>
        private readonly Dictionary<int, LinearAnalyzerLogFactory> logFactories = new Dictionary<int, LinearAnalyzerLogFactory>();

        /// <summary>
        /// The dictionary of logs of the IAnalyzer.
        /// </summary>
        private readonly Dictionary<int, IAnalyzerLog[]> logs = new Dictionary<int, IAnalyzerLog[]>();

        /// <summary>
        /// The solver.
        /// </summary>
        private ISolver solver;

        /// <summary>
        /// The dictionary of subdomains.
        /// </summary>
        private Dictionary<int, ISolverSubdomain> subdomains;

        /// <summary>
        /// The parent analyzer.
        /// </summary>
        public IAnalyzer ParentAnalyzer
        {
            get { return (IAnalyzer)parentAnalyzer; }
            set { parentAnalyzer = (INonLinearParentAnalyzer)value; }
        }

        /// <summary>
        /// The child analyzer. This analyzer doesn't have a child analyzer.
        /// </summary>
        public IAnalyzer ChildAnalyzer
        {
            get { return null; }
            set { throw new InvalidOperationException("Newton-Raphson analyzer cannot contain an embedded analyzer."); }
        }

        public Dictionary<int, IAnalyzerLog[]> Logs
        {
            get
            {
               return this.logs;
            }
        }

        public Dictionary<int, LinearAnalyzerLogFactory> LogFactories { get { return logFactories; } }

        /// <summary>
        /// Factory method to create an instance of <see cref="NonLinearAnalyzerNewtonRaphsonNew"/> with
        /// a fixed incremental load factor.
        /// </summary>
        /// <param name="solver">The linear equation system solver.</param>
        /// <param name="subdomains">The subdomain dictionary.</param>
        /// <param name="provider">The provider.</param>
        /// <param name="increments">The number of load increments.</param>
        /// <param name="freedomDegreeCount"></param>
        /// <returns>The newly created instance.</returns>
        public static NonLinearAnalyzerNewtonRaphsonNew NonLinearAnalyzerWithFixedLoadIncrements(ISolver solver,
            Dictionary<int, ISolverSubdomain> subdomains, INonLinearProvider provider, int increments,
            int freedomDegreeCount)
        {
            IList<double> loadIncrements = new List<double>(increments);
            double loadIncrementFactor = 1.0 / increments;
            for (int i = 0; i < increments; i++)
            {
                loadIncrements.Add(loadIncrementFactor);
            }

            return new NonLinearAnalyzerNewtonRaphsonNew(solver, subdomains, provider, increments, freedomDegreeCount, loadIncrements);
        }

        /// <summary>
        /// Factory method to create an instance of <see cref="NonLinearAnalyzerNewtonRaphsonNew"/> with
        /// a given list of incremental load factors.
        /// </summary>
        /// <param name="solver">The linear equation system solver.</param>
        /// <param name="subdomains">The subdomain dictionary.</param>
        /// <param name="provider">The provider.</param>
        /// <param name="freedomDegreeCount"></param>
        /// <param name="loadIncrements">The list of incremental load factors.</param>
        /// <returns>The newly created instance.</returns>
        public static NonLinearAnalyzerNewtonRaphsonNew nonLinearAnalyzerWithPrescribedIncrements(ISolver solver,
                Dictionary<int, ISolverSubdomain> subdomains, INonLinearProvider provider, int freedomDegreeCount,
                IList<double> loadIncrements)
        {
            int increments = loadIncrements.Count;
            return new NonLinearAnalyzerNewtonRaphsonNew( solver,subdomains, provider,increments,freedomDegreeCount,loadIncrements);
        }

	    private NonLinearAnalyzerNewtonRaphsonNew(ISolver solver, Dictionary<int, ISolverSubdomain> subdomains,
            INonLinearProvider provider, int increments, int freedomDegreeCount,
            IList<double> loadIncrements)
        {
            this.solver = solver;
            this.subdomains = subdomains;
            this.provider = provider;
            this.increments = increments;
            this.globalRightHandSide = new Vector<double>(freedomDegreeCount);
            this.loadIncrementsList = loadIncrements;
            foreach(var subdomain in  subdomains.Values)
            {
                this.displacementMap[subdomain.ID] = new Vector<double>(subdomain.RHS.Length);
            }
        }

        public void Solve()
        {
            InitializeLogs();
            displacements = new double[increments];
            DateTime start = DateTime.Now;
            this.loadFactor = 0.0;
            this.initializeInternalVectors();

            for (int increment = 0; increment < this.increments; increment++)
            {
                this.incrementalLoadFactor = this.loadIncrementsList[increment];
                this.loadFactor += this.incrementalLoadFactor;
                this.clearIncrementalSolutionVector();
                this.UpdateRightHandSide();

                for (int step = 0; step < MAX_ITERATIOMS; step++)
                {
                                  
                    this.solver.Solve();
                    double errorNorm = this.CalculateInternalRightHandSideNorm() / ( this.rightHandSideNorm); //correction of the convergence criterion

                    if ((Double.IsNaN(errorNorm)) || (errorNorm >= MAX_RESIDUAL_NORM_ALLOWED))
                    {
                        throw new ArithmeticException("The redisual norm exceeded the maximum allowed value.");
                    }
                    if (Math.Abs(errorNorm) < this.Tolerance)
                    {
                        break;
                    }
                   //Console.WriteLine("NR {0}, iteration: {1}, error: {2}", increment, step, errorNorm);

                    this.SplitResidualForcesToSubdomains();
                    if (((step + 1) % this.StepForMatrixRebuild) == 0)
                    {
                        provider.Reset();
                        this.BuildMatrices();
                        this.solver.Initialize();
                    }

                }
                this.SaveMaterialStateAndUpdateSolution();
                //Console.WriteLine(this.displacementMap[1][dofid]);
                displacements[increment] = this.displacementMap[1][dofid];
            }
            this.copySolutionToSubdomains();
            DateTime end = DateTime.Now;
            StoreLogResults(start, end);
        }
        //public void Solve(int currentincrement, int totalsteps)
        //{
        //    if (currentincrement == 0)
        //    {
        //        InitializeLogs();
        //        displacements = new double[totalsteps];
        //        //DateTime start = DateTime.Now;
        //        this.loadFactor = 0.0;
        //        this.initializeInternalVectors();
        //    }
        //    this.incrementalLoadFactor = this.loadIncrementsList[currentincrement];
        //    this.loadFactor += this.incrementalLoadFactor;
        //    this.clearIncrementalSolutionVector();
        //    this.UpdateRightHandSide();

        //    for (int step = 0; step < MAX_ITERATIOMS; step++)
        //    {

        //        this.solver.Solve();
        //        double errorNorm = this.CalculateInternalRightHandSideNorm() / (this.rightHandSideNorm); //correction of the convergence criterion

        //        if ((Double.IsNaN(errorNorm)) || (errorNorm >= MAX_RESIDUAL_NORM_ALLOWED))
        //        {
        //            throw new ArithmeticException("The redisual norm exceeded the maximum allowed value.");
        //        }
        //        if (Math.Abs(errorNorm) < this.Tolerance)
        //        {
        //            break;
        //        }
        //     //    Console.WriteLine("NR {0}, iteration: {1}, error: {2}", currentincrement, step, errorNorm);

        //        this.SplitResidualForcesToSubdomains();
        //        if (((step + 1) % this.StepForMatrixRebuild) == 0)
        //        {
        //            provider.Reset();
        //            this.BuildMatrices();
        //            this.solver.Initialize();
        //        }

        //    }
        //    this.SaveMaterialStateAndUpdateSolution();
        //   //  Console.WriteLine(this.displacementMap[1][dofid]);
        //    displacements[currentincrement] = this.displacementMap[1][dofid];


        //        this.copySolutionToSubdomains(); //REMEMBER THE IF CLAUSE (CURRENTINCREMENT==TOTALSTEPS-1) NO NEED THIS IF CLAUSE

        //    //DateTime end = DateTime.Now;
        //    //StoreLogResults(start, end);
        //}
        public void Solve(int currentincrement, int totalsteps)
        {
            if (fail == 0)
            {
                if (currentincrement == 0)
                {
                    InitializeLogs();
                    displacements = new double[totalsteps];
                    //DateTime start = DateTime.Now;
                    this.loadFactor = 0.0;
                    this.initializeInternalVectors();
                }
                this.incrementalLoadFactor = this.loadIncrementsList[currentincrement];
                this.loadFactor += this.incrementalLoadFactor;
                this.clearIncrementalSolutionVector();
                this.UpdateRightHandSide();

                for (int step = 0; step < MAX_ITERATIOMS; step++)
                {

                    this.solver.Solve();
                    double errorNorm = this.CalculateInternalRightHandSideNorm() / (this.rightHandSideNorm); //correction of the convergence criterion

                    if ((Double.IsNaN(errorNorm)) || (errorNorm >= MAX_RESIDUAL_NORM_ALLOWED))
                    {
                        fail = 1;
                        failinc = currentincrement;
                        break;
                    }
                    if (Math.Abs(errorNorm) < this.Tolerance)
                    {
                        break;
                    }
                    //    Console.WriteLine("NR {0}, iteration: {1}, error: {2}", currentincrement, step, errorNorm);

                    this.SplitResidualForcesToSubdomains();
                    if (((step + 1) % this.StepForMatrixRebuild) == 0)
                    {
                        provider.Reset();
                        this.BuildMatrices();
                        this.solver.Initialize();
                    }

                }
            }
            if (fail == 0)
            {
                this.SaveMaterialStateAndUpdateSolution();
                Console.WriteLine(currentincrement);
                Console.WriteLine(this.displacementMap[1][dofid]);
                displacements[currentincrement] = this.displacementMap[1][dofid];


                this.copySolutionToSubdomains(); //REMEMBER THE IF CLAUSE (CURRENTINCREMENT==TOTALSTEPS-1) NO NEED THIS IF CLAUSE

                //DateTime end = DateTime.Now;
                //StoreLogResults(start, end);
            }
        }
        private double CalculateInternalRightHandSideNorm()
        {
            this.globalRightHandSide.Clear();
            foreach (var subdomain in this.subdomains.Values)
            {
                int id = subdomain.ID;
                Vector<double> iterativeDisplacements = this.iterativeDisplacementMap[id];
                iterativeDisplacements.Clear();
                iterativeDisplacements.Add(subdomain.Solution as Vector<double>);

                Vector<double> incrementalDisplacements = this.incrementalDisplacementMap[id];
                incrementalDisplacements.Add(iterativeDisplacements);

                Vector<double> totalDisplacements = new Vector<double>(incrementalDisplacements.Length);
                Vector<double> displacements = this.displacementMap[id];
                totalDisplacements.Add(displacements);
                totalDisplacements.Add(incrementalDisplacements);

                Vector<double> internalRightHandSide = subdomain.GetRHSFromSolution(displacements, incrementalDisplacements) as Vector<double>;
                this.provider.ProcessInternalRHS(subdomain, internalRightHandSide.Data, totalDisplacements.Data);

                if (this.parentAnalyzer != null)
                {
                    Vector<double> additionalRightHandside =new Vector<double>(this.parentAnalyzer.GetOtherRHSComponents(subdomain, totalDisplacements));
                    internalRightHandSide.Add(additionalRightHandside);
                }

                Vector<double> subdomainRightHandSide = subdomain.RHS as Vector<double>;
                Vector<double> rightHandSide = this.rightHandSides[id];
                subdomainRightHandSide.Clear();
                subdomainRightHandSide.Add(rightHandSide);
                subdomainRightHandSide.Multiply(this.loadFactor);
                subdomainRightHandSide.Subtract(internalRightHandSide);
                subdomain.SubdomainToGlobalVector(subdomainRightHandSide.Data, this.globalRightHandSide.Data);
            }
            return this.provider.RHSNorm(this.globalRightHandSide.Data);
        }


        private void clearIncrementalSolutionVector()
        {
            foreach (var subdomain in this.subdomains.Values)
            {
                this.incrementalDisplacementMap[subdomain.ID].Clear();
            }
        }

        private void ClearMaterialStresses()
        {
            foreach (var subdomain in this.subdomains.Values)
            {
                subdomain.ClearMaterialStresses();
            }
        }


        private void copySolutionToSubdomains()
        {
            foreach (var subdomain in this.subdomains.Values)
            {
                int id = subdomain.ID;
                double[] displacementData = new double[this.displacementMap[id].Length];
                this.displacementMap[id].CopyTo(displacementData, 0);
                subdomain.Solution = new Vector<double>(displacementData);
            }
        }


        private void initializeInternalVectors()
        {
            this.globalRightHandSide.Clear();
            this.rightHandSides.Clear();
            this.incrementalDisplacementMap.Clear();
            this.iterativeDisplacementMap.Clear();

            foreach (var subdomain in this.subdomains.Values)
            {
                int vectorLength = subdomain.RHS.Length;
                int id = subdomain.ID;
                double[] vecRData = new double[vectorLength];
                subdomain.RHS.CopyTo(vecRData, 0);
                Vector<double> vecR = new Vector<double>(vecRData);
                this.rightHandSides[id] = vecR;
                this.incrementalDisplacementMap[id] = new Vector<double>(vectorLength);
                this.iterativeDisplacementMap[id] = new Vector<double>(vectorLength);
                subdomain.SubdomainToGlobalVector((subdomain.RHS as Vector<double>).Data, this.globalRightHandSide.Data);
            }

            this.rightHandSideNorm = this.provider.RHSNorm(this.globalRightHandSide.Data);
        }


        private void SaveMaterialStateAndUpdateSolution()
        {
            foreach (var subdomain in this.subdomains.Values)
            {
                subdomain.SaveMaterialState();
                int id = subdomain.ID;
                Vector<double> displacements = this.displacementMap[id];
                Vector<double> incrementalDisplacements = this.incrementalDisplacementMap[id];
                displacements.Add(incrementalDisplacements);
            }
        }


        private void SplitResidualForcesToSubdomains()
        {
            foreach (var subdomain in this.subdomains.Values)
            {
                Vector<double> subdomainRightHandSide = subdomain.RHS as Vector<double>;
                subdomainRightHandSide.Clear();
                subdomain.SplitGlobalVectorToSubdomain(this.globalRightHandSide.Data, subdomainRightHandSide.Data);
            }
        }


        private void UpdateRightHandSide()
        {
            this.globalRightHandSide.Clear();
            foreach (var subdomain in this.subdomains.Values)
            {
                int id = subdomain.ID;
                IVector<double> displacements = this.displacementMap[id];
                IVector<double> incrementalDisplacements = this.incrementalDisplacementMap[id];
                Vector<double> internalRightHandSide = subdomain.GetRHSFromSolution(displacements, incrementalDisplacements) as Vector<double>;

                if (this.parentAnalyzer != null)
                {
                    Vector<double> additionalRightHandside = new Vector<double>(
                            this.parentAnalyzer.GetOtherRHSComponents(subdomain, displacements));
                    internalRightHandSide.Add(additionalRightHandside);
                }

                Vector<double> subdomainRightHandSide = subdomain.RHS as Vector<double>;
                Vector<double> rightHandSide = this.rightHandSides[id];
                subdomainRightHandSide.Clear();
                subdomainRightHandSide.Add(rightHandSide);
                subdomainRightHandSide.Multiply(this.loadFactor);
                subdomainRightHandSide.Subtract(internalRightHandSide);
                subdomain.SubdomainToGlobalVector(subdomainRightHandSide.Data, this.globalRightHandSide.Data);
            }
        }

        public void BuildMatrices()
        {
            if (this.parentAnalyzer == null)
                throw new InvalidOperationException("This Newton-Raphson non-linear analyzer has no parent.");

            this.parentAnalyzer.BuildMatrices();
            this.solver.Initialize();
        }

        public void Initialize()
        {
            this.solver.Initialize();
        }

        private void InitializeLogs()
        {
            logs.Clear();
            foreach (int id in logFactories.Keys) logs.Add(id, logFactories[id].CreateLogs());
        }

        private void StoreLogResults(DateTime start, DateTime end)
        {
            foreach (int id in logs.Keys)
                foreach (var l in logs[id])
                    l.StoreResults(start, end, subdomains[id].Solution);
        }
    }
}
