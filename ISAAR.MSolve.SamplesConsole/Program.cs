using System;
using ISAAR.MSolve.Analyzers;
using ISAAR.MSolve.Analyzers.Dynamic;
using ISAAR.MSolve.Discretization.FreedomDegrees;
using ISAAR.MSolve.FEM.Entities;
using ISAAR.MSolve.Logging;
using ISAAR.MSolve.Problems;
using ISAAR.MSolve.Solvers;
using ISAAR.MSolve.Solvers.Direct;
using ISAAR.MSolve.LinearAlgebra;
using MGroup.Stochastic;
using MGroup.Stochastic.Structural;
using MGroup.Stochastic.Structural.Example;
using ISAAR.MSolve.Analyzers.NonLinear;
using ISAAR.MSolve.PreProcessor;
using ISAAR.MSolve.PreProcessor.Materials;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.IO;
using System.Threading;
using System.Threading.Tasks;
using ISAAR.MSolve.Analyzers.Interfaces;
using System.Diagnostics;

namespace ISAAR.MSolve.SamplesConsole
{
    class Program
    {
        private const int subdomainID = 0;
        #region Mgroupstuff
        static void MgroupMain(string[] args)
        {
            //SolveBuildingInNoSoilSmall();
            //TrussExample.Run();
            //FEM.Cantilever2D.Run();
            //FEM.Cantilever2DPreprocessor.Run();
            //FEM.WallWithOpenings.Run();
            //SeparateCodeCheckingClass.Check06();
            //SolveBuildingInNoSoilSmall();
            //SolveBuildingInNoSoilSmallDynamic();
            //SolveStochasticMaterialBeam2DWithBruteForceMonteCarlo();
            //CNTExamples.CNT_4_4_DisplacementControl();
            //CNTExamples.CNT_4_4_NewtonRaphson();
            //Tests.FEM.Shell8andCohesiveNonLinear.RunTest();
            //AppliedDisplacementExample.Run();

            //Logging.PrintForceDisplacementCurve.CantileverBeam2DCorotationalLoadControl();

            //SuiteSparseBenchmarks.MemoryConsumptionDebugging();
            //SolverBenchmarks.SuiteSparseMemoryConsumptionDebugging();
            //NRNLAnalyzerDevelopTest.SolveDisplLoadsExample();
            //SeparateCodeCheckingClass4.Check05bStressIntegrationObje_Integration();
            //SeparateCodeCheckingClass4.Check_Graphene_rve_Obje_Integration();
            //IntegrationElasticCantileverBenchmark.RunExample();
            //OneRveExample.Check_Graphene_rve_serial();
            //BondSlipTest.CheckStressStrainBonSlipMaterial();
            //OneRveExample.Check_Graphene_rve_parallel();
            //LinearRves.CheckShellScaleTransitionsAndMicrostructure();
            //SolveCantileverWithStochasticMaterial();

            //MeshPartitioningExamples.PartitionMeshes();

        }

        private static void SolveBuildingInNoSoilSmall()
        {
            var model = new Model();
            model.SubdomainsDictionary.Add(subdomainID, new Subdomain(subdomainID));
            BeamBuildingBuilder.MakeBeamBuilding(model, 20, 20, 20, 5, 4, model.NodesDictionary.Count + 1,
                model.ElementsDictionary.Count + 1, subdomainID, 4, false, false);
            model.Loads.Add(new Load() { Amount = -100, Node = model.Nodes[21], DOF = StructuralDof.TranslationX });

            // Solver
            var solverBuilder = new SkylineSolver.Builder();
            ISolver solver = solverBuilder.BuildSolver(model);

            // Structural problem provider
            var provider = new ProblemStructural(model, solver);

            // Linear static analysis
            var childAnalyzer = new LinearAnalyzer(model, solver, provider);
            var parentAnalyzer = new StaticAnalyzer(model, solver, provider, childAnalyzer);

            // Request output
            int monitorDof = 420;
            childAnalyzer.LogFactories[subdomainID] = new LinearAnalyzerLogFactory(new int[] { monitorDof });

            // Run the analysis
            parentAnalyzer.Initialize();
            parentAnalyzer.Solve();

            // Write output
            DOFSLog log = (DOFSLog)childAnalyzer.Logs[subdomainID][0]; //There is a list of logs for each subdomain and we want the first one
            Console.WriteLine($"dof = {monitorDof}, u = {log.DOFValues[monitorDof]}");
        }

        private static void SolveBuildingInNoSoilSmallDynamic()
        {
            var model = new Model();
            model.SubdomainsDictionary.Add(subdomainID, new Subdomain(subdomainID));
            BeamBuildingBuilder.MakeBeamBuilding(model, 20, 20, 20, 5, 4, model.NodesDictionary.Count + 1,
                model.ElementsDictionary.Count + 1, subdomainID, 4, false, false);

            // Solver
            var solverBuilder = new SkylineSolver.Builder();
            ISolver solver = solverBuilder.BuildSolver(model);

            // Structural problem provider
            var provider = new ProblemStructural(model, solver);

            // Linear static analysis
            var childAnalyzer = new LinearAnalyzer(model, solver, provider);
            var parentAnalyzerBuilder = new NewmarkDynamicAnalyzer.Builder(model, solver, provider, childAnalyzer, 0.01, 0.1);
            parentAnalyzerBuilder.SetNewmarkParametersForConstantAcceleration(); // Not necessary. This is the default
            NewmarkDynamicAnalyzer parentAnalyzer = parentAnalyzerBuilder.Build();

            // Request output
            int monitorDof = 420;
            childAnalyzer.LogFactories[subdomainID] = new LinearAnalyzerLogFactory(new int[] { monitorDof });

            // Run the analysis
            parentAnalyzer.Initialize();
            parentAnalyzer.Solve();

            // Write output
            DOFSLog log = (DOFSLog)childAnalyzer.Logs[subdomainID][0]; //There is a list of logs for each subdomain and we want the first one
            Console.WriteLine($"dof = {monitorDof}, u = {log.DOFValues[monitorDof]}");

            //TODO: No loads have been defined so the result is bound to be 0.
        }

        private static void SolveCantileverWithStochasticMaterial()
        {
            const int iterations = 1000;
            const double youngModulus = 2.1e8;

            var domainMapper = new CantileverStochasticDomainMapper(new[] { 0d, 0d, 0d });
            var evaluator = new StructuralStochasticEvaluator(youngModulus, domainMapper);
            var m = new MonteCarlo(iterations, evaluator, evaluator);
            m.Evaluate();
        }
        #endregion
        #region HexaProgramCode
        public static double[,] Stoch1;
        public static double[,] Stoch2;
        public static double[,] Stoch3;
        public static int montecarlosim;
        public static int indexbegin;
        public static bool hasfailed = false;
        public static int stepoffail = 0;
        public static double dispfail = 0.0;
        public static double[] increments;
        public static double[] dispstoch;
        public static double[] dispstoch150;
        public static double[] stresstoch;
        public static int numofsample = 0;
        #region readwritemethods
        public static void readData(string DataFileName, out double[] array)
        {
            string dataLine;
            string[] dataFields;
            string[] numSeparators1 = { ":" };
            string[] numSeparators2 = { " " };
            StreamReader rStream;
            rStream = File.OpenText(DataFileName);
            int dim = 1;
            dataLine = rStream.ReadLine();
            dataFields = dataLine.Split(numSeparators1, StringSplitOptions.RemoveEmptyEntries);
            dim = int.Parse(dataFields[0]);
            array = new double[dim];
            for (int i = 0; i < dim; i++)
            {
                dataLine = rStream.ReadLine();
                dataFields = dataLine.Split(numSeparators1, StringSplitOptions.RemoveEmptyEntries);
                array[i] = double.Parse(dataFields[0]);
            }
            rStream.Close();

        }
        public static void writeData(double[] array, int identifier)
        {
            string filename, dataLine;
            // The identifier is for telling if you want to write the whole array (1) or the last element (0) (for example the whole displacement curve or the last increment)
            // To insert spaces, use the simple space character " ", not tabs (i.e. "\t"). 
            // the editors do not 'interpret' the tabs in the same way, 
            // so if you open the file with different editors can be a mess.
            //string spaces1 = "        ";
            //string spaces2 = "              ";

            // format specifier to write the real numbers
            string fmtSpecifier = "{0: 0.0000E+00;-0.0000E+00}";

            StreamWriter wStream;
            filename = "displacements.txt";
            wStream = File.CreateText(filename);
            if (identifier == 1)
            {
                for (int i = 0; i < array.GetLength(0); i++)
                {
                    dataLine = String.Format(fmtSpecifier, array[i]);
                    wStream.WriteLine(dataLine);
                }
                wStream.Close();
            }
            else
            {
                dataLine = String.Format(fmtSpecifier, array[array.GetLength(0) - 1]);
                wStream.WriteLine(dataLine);
                wStream.Close();
            }
        }
        public static void writeData(double[] array, int identifier, string filename)
        {
            string dataLine;
            // The identifier is for telling if you want to write the whole array (1) or the last element (0) (for example the whole displacement curve or the last increment)
            // To insert spaces, use the simple space character " ", not tabs (i.e. "\t"). 
            // the editors do not 'interpret' the tabs in the same way, 
            // so if you open the file with different editors can be a mess.
            //string spaces1 = "        ";
            //string spaces2 = "              ";

            // format specifier to write the real numbers
            string fmtSpecifier = "{0: 0.0000E+00;-0.0000E+00}";

            StreamWriter wStream;

            wStream = File.CreateText(filename);
            if (identifier == 1)
            {
                for (int i = 0; i < array.GetLength(0); i++)
                {
                    dataLine = String.Format(fmtSpecifier, array[i]);
                    wStream.WriteLine(dataLine);
                }
                wStream.Close();
            }
            else
            {
                dataLine = String.Format(fmtSpecifier, array[array.GetLength(0) - 1]);
                wStream.WriteLine(dataLine);
                wStream.Close();
            }
        }
        public static void writeTime(DateTime begintime, DateTime endtime)
        {
            string filename;

            StreamWriter wStream;
            filename = "time.txt";
            wStream = File.CreateText(filename);
            wStream.WriteLine(begintime);
            wStream.WriteLine(endtime);
            wStream.Close();
        }
        public static void readMatrixData(string DataFileName, out double[,] array)
        {
            string dataLine;
            string[] dataFields;
            string[] numSeparators1 = { ":" };
            string[] numSeparators2 = { " " };
            StreamReader rStream;
            rStream = File.OpenText(DataFileName);
            int dim = 1;
            int dim1 = 1;
            dataLine = rStream.ReadLine();
            dataFields = dataLine.Split(numSeparators1, StringSplitOptions.RemoveEmptyEntries);
            dim = int.Parse(dataFields[0]);
            dataLine = rStream.ReadLine();
            dataFields = dataLine.Split(numSeparators1, StringSplitOptions.RemoveEmptyEntries);
            dim1 = int.Parse(dataFields[0]);
            double[,] array1 = new double[dim, dim1];
            for (int i = 0; i < dim; i++)
            {
                dataLine = rStream.ReadLine();
                dataFields = dataLine.Split(numSeparators1, StringSplitOptions.RemoveEmptyEntries);
                for (int j = 0; j < dim1; j++)
                {
                    array1[i, j] = double.Parse(dataFields[j]);
                }
            }
            rStream.Close();
            array = array1;
        }
        public static double[] readMatrixDataPartially(double[,] Matrix, int rowbegin, int rowend, int colbegin, int colend)
        {
            var array = new double[(rowend + 1 - rowbegin) * (colend + 1 - colbegin)];
            var k = 0;
            for (int i = rowbegin; i < rowend + 1; i++)
            {
                for (int j = colbegin; j < colend + 1; j++)
                {
                    array[k] = Matrix[i, j];
                    k++;
                }
            }
            return array;
        }
        #endregion
        #region solvermethods
        private static void CollectMonteCarloFailDetails(double lambda)
        {
            Console.WriteLine("Displacement of failure: ");
            Console.WriteLine(dispfail);
            Console.WriteLine("Fail load in Kpa");
            Console.WriteLine(lambda * 150.0);
            Console.WriteLine("End of {0} Monte Carlo Sample.",numofsample+1);
            numofsample++;
        }
        private static void SolveHexaSoil()
        {
        }
        private static void SolveStochasticHexaSoil(int samplenumber, double[] Stoch1, double[] Stoch2, double[] Stoch3, double[] omega,double lambda)
        {
            Model model = new Model();
            model.SubdomainsDictionary.Add(1, new Subdomain(1));

            HexaSoil2.MakeHexaSoil(model, Stoch1, Stoch2, Stoch3, omega,lambda);

            model.ConnectDataStructures();

            var solverBuilder = new SuiteSparseSolver.Builder();
            //var solverBuilder = new SkylineSolver.Builder();
            ISolver solver = solverBuilder.BuildSolver(model);

            var provider = new ProblemPorous(model, solver);
            LibrarySettings.LinearAlgebraProviders = LinearAlgebraProviderChoice.MKL;

            int increments = 1;
            var childAnalyzerBuilder = new LoadControlAnalyzer.Builder(model, solver, provider, increments);
            childAnalyzerBuilder.ResidualTolerance = 1E-5;
            childAnalyzerBuilder.MaxIterationsPerIncrement = 100;
            childAnalyzerBuilder.NumIterationsForMatrixRebuild = 1;
            LoadControlAnalyzer childAnalyzer = childAnalyzerBuilder.Build();
            var parentAnalyzerBuilder = new NewmarkDynamicAnalyzer.Builder(model, solver, provider, childAnalyzer, 0.001, 1);
            NewmarkDynamicAnalyzer parentAnalyzer = parentAnalyzerBuilder.Build();
           
            parentAnalyzer.Initialize();
            parentAnalyzer.Solve();
            int monitorDof = HexaSoil2.ProvideIdMonitor(model);
            Node nn = model.NodesDictionary[203];
            var hhhh = model.GlobalDofOrdering.GlobalFreeDofs[nn, StructuralDof.TranslationZ];
            if (model.Subdomains[0].hasfailed==true)
            {
                dispfail = childAnalyzer.dispfail;
                stepoffail = parentAnalyzer.failstep;
                hasfailed = true;
                Console.WriteLine("XXXXXXXXXXXXXXXXXXXX"); //In order to erase it as a previous iteration.
            }
            var d = 0;
            //dispstoch[samplenumber] = analyzer.displacements[analyzer.failinc-1];
            //dispstoch150[samplenumber] = analyzer.displacements[149];
            //stresstoch[samplenumber] = analyzer.failinc;
            //writeData(analyzer.displacements, 1);
            //Hexa8 h1 = (Hexa8)model.Elements[92].ElementType;
            //IFiniteElementMaterial3D[] matGP = new IFiniteElementMaterial3D[6];
            //matGP = h1.materialsAtGaussPoints;
            //var StressesCheck = matGP[7].Stresses[2];
            //stresstoch[samplenumber] = StressesCheck;
        }
        #endregion
        //static void Main(string[] args)
        //{
        //    SolveHexaSoil();
        //}
        static void Main(string[] args)
        {
            DateTime begin = DateTime.Now;
            readMatrixData("input1.txt", out Stoch1);
            readMatrixData("input2.txt", out Stoch2);
            readMatrixData("input3.txt", out Stoch3);
            Console.WriteLine("Provide the initial index. Dont forget we have zero indexing.");
            indexbegin = Int32.Parse(Console.ReadLine());
            Console.WriteLine("Provide the final index.");
            montecarlosim = Int32.Parse(Console.ReadLine());
            StreamWriter writer;
            TextWriter oldOut = Console.Out;
            writer = File.CreateText("Output.txt");
            Console.SetOut(writer);
            dispstoch = new double[montecarlosim - indexbegin];
            dispstoch150 = new double[montecarlosim - indexbegin];
            stresstoch = new double[montecarlosim - indexbegin];
            for (int index = indexbegin; index < montecarlosim; index++)
            {
               double lambda = 1000.0 / 150.0;
                double lambdaprev = 1.1*lambda;
                double maxlambdaofnofailure = 0.0;
                double thislambdac = lambda;
                double previouslambdac = 0.0;
                bool isfirstiter = true;
                hasfailed = false;
                while (!(hasfailed == true && Math.Abs((thislambdac - previouslambdac) / thislambdac)<0.01))
                {
                    if (hasfailed == true)
                    {
                        lambdaprev = lambda;
                        if (maxlambdaofnofailure == 0)
                        {
                            lambda = 0.5 * (lambda * stepoffail / 1000.0 + lambdaprev);
                        }
                        else
                        {
                            lambda = 0.5 * (maxlambdaofnofailure + lambdaprev);
                        }
                    }
                    else
                    {
                        if (isfirstiter)
                        {
                            isfirstiter = false;
                        }
                        else
                        {
                            if (maxlambdaofnofailure != 0 && maxlambdaofnofailure < lambda)
                            {
                                maxlambdaofnofailure = lambda;
                            }
                            if (maxlambdaofnofailure == 0)                             
                            {
                                maxlambdaofnofailure = lambda;
                                lambda = 0.5 * (lambda + lambdaprev);
                            }
                            else
                            {
                                lambda = 0.5 * (maxlambdaofnofailure + lambdaprev);
                            }
                        }
                    }
                    Debug.WriteLine("Previous Lambda {0} Current Lambda {1}", lambdaprev, lambda);
                    hasfailed = false;
                    SolveStochasticHexaSoil(index, readMatrixDataPartially(Stoch1, index, index, 0, 7), readMatrixDataPartially(Stoch2, index, index, 0, 7), readMatrixDataPartially(Stoch3, index, index, 0, 7), readMatrixDataPartially(Stoch3, 0, 7, 8, 8), lambda);
                    if (hasfailed == true)
                    {
                        previouslambdac = thislambdac;
                        thislambdac = lambdaprev;
                    }
                }
                CollectMonteCarloFailDetails(lambda);
            }
            //for (int i = 0; i < 1; i++)
            //{
            //    SolveStochasticHexaSoil(1, Stoch1[1], Stoch2[1]);
            //}
            //Parallel.For(0, montecarlosim-1,
            //      index =>
            //      {
            //          SolveStochasticHexaSoil(index, Stoch1[index], Stoch2[index], readMatrixDataPartially(Stoch3, index, index, 0, 7), readMatrixDataPartially(Stoch3, 0, 7, 8, 8));
            //          Console.WriteLine(index);
            //      });
            DateTime end = DateTime.Now;
            writeTime(begin, end);
            //writeData(dispstoch, 1, "displacements.txt");
            //writeData(dispstoch150, 1, "displacements150.txt");
            //writeData(stresstoch, 1, "stresses.txt");
            Console.SetOut(oldOut);
            writer.Close();
        }
        #endregion
    }
}
