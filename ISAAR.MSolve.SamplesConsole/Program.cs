using ISAAR.MSolve.Analyzers;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.FEM.Elements;
using ISAAR.MSolve.FEM.Entities;
using ISAAR.MSolve.FEM.Materials;
using ISAAR.MSolve.Logging;
using ISAAR.MSolve.Materials.Interfaces;
using ISAAR.MSolve.Numerical.LinearAlgebra;
using ISAAR.MSolve.Problems;
using ISAAR.MSolve.Solvers.Interfaces;
using ISAAR.MSolve.Solvers.Skyline;
using MGroup.Stochastic;
using MGroup.Stochastic.Structural;
using MGroup.Stochastic.Structural.Example;
using System.Collections.Generic;


using ISAAR.MSolve.FEM;
using ISAAR.MSolve.FEM.Elements.SupportiveClasses;
using ISAAR.MSolve.FEM.Interfaces;
using Moq;
using System;
using System.Text;
using Xunit;
using ISAAR.MSolve.FEM.Entities.TemporalFunctions;
using System.IO;

namespace ISAAR.MSolve.SamplesConsole
{
    class Program
    {
        public static double[] Stoch1;
        public static double[] Stoch2;
        public static double[] Stoch3;
        public static double[] increments;
        public static int montecarlosim = 1000;
        public static double[] dispstoch = new double[montecarlosim];
        public static double[] stresstoch = new double[montecarlosim];
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
            string[] numSeparators2 = { ":" };
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
        private static void SolveBuildingInNoSoilSmall()
        {
            VectorExtensions.AssignTotalAffinityCount();
            Model model = new Model();
            model.SubdomainsDictionary.Add(1, new Subdomain() { ID = 1 });
            BeamBuildingBuilder.MakeBeamBuilding(model, 20, 20, 20, 5, 4, model.NodesDictionary.Count + 1,
                model.ElementsDictionary.Count + 1, 1, 4, false, false);
            model.Loads.Add(new Load() { Amount = -100, Node = model.Nodes[21], DOF = DOFType.X });
            model.ConnectDataStructures();

            var linearSystems = new Dictionary<int, ILinearSystem>(); //I think this should be done automatically
            linearSystems[1] = new SkylineLinearSystem(1, model.Subdomains[0].Forces);
            SolverSkyline solver = new SolverSkyline(linearSystems[1]);
            ProblemStructural provider = new ProblemStructural(model, linearSystems);
            LinearAnalyzer analyzer = new LinearAnalyzer(solver, linearSystems);
            StaticAnalyzer parentAnalyzer = new StaticAnalyzer(provider, analyzer, linearSystems);

            analyzer.LogFactories[1] = new LinearAnalyzerLogFactory(new int[] { 420 });

            parentAnalyzer.BuildMatrices();
            parentAnalyzer.Initialize();
            parentAnalyzer.Solve();
        }

        private static void SolveBuildingInNoSoilSmallDynamic()
        {
            VectorExtensions.AssignTotalAffinityCount();
            Model model = new Model();
            model.SubdomainsDictionary.Add(1, new Subdomain() { ID = 1 });
            BeamBuildingBuilder.MakeBeamBuilding(model, 20, 20, 20, 5, 4, model.NodesDictionary.Count + 1,
                model.ElementsDictionary.Count + 1, 1, 4, false, false);
            model.ConnectDataStructures();

            var linearSystems = new Dictionary<int, ILinearSystem>(); //I think this should be done automatically
            linearSystems[1] = new SkylineLinearSystem(1, model.Subdomains[0].Forces);
            SolverSkyline solver = new SolverSkyline(linearSystems[1]);
            ProblemStructural provider = new ProblemStructural(model, linearSystems);
            LinearAnalyzer analyzer = new LinearAnalyzer(solver, linearSystems);
            NewmarkDynamicAnalyzer parentAnalyzer = new NewmarkDynamicAnalyzer(provider, analyzer, linearSystems, 0.25, 0.5, 0.01, 0.1);

            analyzer.LogFactories[1] = new LinearAnalyzerLogFactory(new int[] { 420 });

            parentAnalyzer.BuildMatrices();
            parentAnalyzer.Initialize();
            parentAnalyzer.Solve();
        }

        private static void SolveStochasticMaterialBeam2DWithBruteForceMonteCarlo()
        {
            #region Beam2D Geometry Data
            VectorExtensions.AssignTotalAffinityCount();
            double youngModulus = 2.0e08;
            double poissonRatio = 0.3;
            double nodalLoad = 10.0;

            IStochasticMaterialCoefficientsProvider coefficientProvider = new PowerSpectrumTargetEvaluatorCoefficientsProvider(10, 0.1, .05, 20, 200, DOFType.X, 0.1, 200, 1e-10);
            StochasticElasticMaterial material = new StochasticElasticMaterial(coefficientProvider)
            {
                YoungModulus = youngModulus,
                PoissonRatio = poissonRatio,
            };

            // Model creation
            Model model = new Model();

            // Add a single subdomain to the model
            model.SubdomainsDictionary.Add(0, new Subdomain() { ID = 0 });

            // Add nodes to the nodes dictonary of the model
            for (int i = 0; i < 11; i++)
            {
                model.NodesDictionary.Add(i, new Node
                {
                    ID = i,
                    X = i * 1,
                    Y = 0,
                    Z = 0
                });
            }

            // Fix cantilever left end node of the model
            model.NodesDictionary[0].Constraints.Add(DOFType.X);
            model.NodesDictionary[0].Constraints.Add(DOFType.Y);
            model.NodesDictionary[0].Constraints.Add(DOFType.RotZ);

            for (int i = 0; i < model.NodesDictionary.Count - 1; i++)
            {
                var element = new Element()
                {
                    ID = i,
                    ElementType = new Beam2DWithStochasticMaterial(material)
                    {
                        SectionArea = 1,
                        MomentOfInertia = 0.1
                    }
                };
                element.AddNode(model.NodesDictionary[i]);
                element.AddNode(model.NodesDictionary[i + 1]);
                model.ElementsDictionary.Add(i, element);
                model.SubdomainsDictionary[0].ElementsDictionary.Add(i, element);
            }


            // Add nodal load values at the right end of the cantilever
            model.Loads.Add(new Load() { Amount = -nodalLoad, Node = model.NodesDictionary[model.NodesDictionary.Count - 1], DOF = DOFType.Y });

            // Needed in order to make all the required data structures
            model.ConnectDataStructures();
            #endregion

            var linearSystems = new Dictionary<int, ILinearSystem>();
            linearSystems[0] = new SkylineLinearSystem(0, model.Subdomains[0].Forces);
            SolverSkyline solver = new SolverSkyline(linearSystems[0]);
            ProblemStructural provider = new ProblemStructural(model, linearSystems);
            LinearAnalyzer childAnalyzer = new LinearAnalyzer(solver, linearSystems);
            StaticAnalyzer parentAnalyzer = new StaticAnalyzer(provider, childAnalyzer, linearSystems);
            MonteCarloAnalyzerWithStochasticMaterial stohasticAnalyzer =
            new MonteCarloAnalyzerWithStochasticMaterial(model, provider, parentAnalyzer, linearSystems,
                    coefficientProvider, 1, 10000);
            stohasticAnalyzer.Initialize();
            stohasticAnalyzer.Solve();

            //Assert.Equal(-2.08333333333333333e-5, stohasticAnalyzer.MonteCarloMeanValue, 8);
        }

        //private static void Solve()
        //{
        //    const int iterations = 1000;
        //    const double youngModulus = 2.1e8;

        //    var domainMapper = new CantileverStochasticDomainMapper(new[] { 0d, 0d, 0d });
        //    var realizer = new GiannisStructuralStochasticRealizer(youngModulus, domainMapper);
        //    var evaluator = new StructuralStochasticEvaluator(youngModulus, domainMapper);
        //    var m = new MonteCarlo(iterations, realizer, evaluator);
        //    m.Evaluate();
        //}

        private static void Solve()
        {
            VectorExtensions.AssignTotalAffinityCount();
            const int iterations = 1000;
            const double youngModulus = 2.1e8;

            var domainMapper = new CantileverStochasticDomainMapper(new[] { 0d, 0d, 0d });
            //var evaluator = new StructuralStochasticEvaluator(youngModulus, domainMapper);
            var evaluator = new StructuralStochasticEvaluator(youngModulus, domainMapper);
            var m = new MonteCarlo(iterations, evaluator, evaluator);
            m.Evaluate();
        }

        public static void TestBeam3DElasticGeometryNonlinearNewmarkDynamicAnalysisExample()
        {
            VectorExtensions.AssignTotalAffinityCount();
            double youngModulus = 210000000.0;
            double poissonRatio = 0.3;
            double nodalLoad = 100.0;
            double area = 0.25 * 0.25;
            double inertiaY = Math.Pow(0.25, 4) / 12;
            double inertiaZ = Math.Pow(0.25, 4) / 12;
            double torsionalInertia = Math.Pow(0.25, 4) / 6;
            double effectiveAreaY = area;
            double effectiveAreaZ = area;
            double density = 7.85;
            int nNodes = 3;
            int nElems = 2;
            int monitorNode = 3;

            // Create new 3D material
            ElasticMaterial3D material = new ElasticMaterial3D
            {
                YoungModulus = youngModulus,
                PoissonRatio = poissonRatio,
            };

            // Node creation
            IList<Node> nodes = new List<Node>();
            Node node1 = new Node { ID = 1, X = 0.0, Y = 0.0, Z = 0.0 };
            Node node2 = new Node { ID = 2, X = 5.0, Y = 0.0, Z = 0.0 };
            Node node3 = new Node { ID = 3, X = 10.0, Y = 0.0, Z = 0.0 };

            nodes.Add(node1);
            nodes.Add(node2);
            nodes.Add(node3);

            // Model creation
            Model model = new Model();

            // Add a single subdomain to the model
            model.SubdomainsDictionary.Add(1, new Subdomain() { ID = 1 });

            // Add nodes to the nodes dictonary of the model
            for (int i = 0; i < nodes.Count; ++i)
            {
                model.NodesDictionary.Add(i + 1, nodes[i]);
            }

            // Constrain bottom nodes of the model
            model.NodesDictionary[1].Constraints.Add(DOFType.X);
            model.NodesDictionary[1].Constraints.Add(DOFType.Y);
            model.NodesDictionary[1].Constraints.Add(DOFType.Z);
            model.NodesDictionary[1].Constraints.Add(DOFType.RotX);
            model.NodesDictionary[1].Constraints.Add(DOFType.RotY);
            model.NodesDictionary[1].Constraints.Add(DOFType.RotZ);
            // Generate elements of the structure
            int iNode = 1;
            for (int iElem = 0; iElem < nElems; iElem++)
            {
                // element nodes
                IList<Node> elementNodes = new List<Node>();
                elementNodes.Add(model.NodesDictionary[iNode]);
                elementNodes.Add(model.NodesDictionary[iNode + 1]);

                // Create new Beam3D section and element
                var beamSection = new BeamSection3D(area, inertiaY, inertiaZ, torsionalInertia, effectiveAreaY, effectiveAreaZ);
                var beam = new EulerBeam3D(youngModulus, poissonRatio, null, null);
                beam.MomentOfInertiaY = inertiaY;
                beam.MomentOfInertiaZ = inertiaZ;
                beam.MomentOfInertiaPolar = torsionalInertia;
                beam.SectionArea = area;
                beam.Density = density;
                // Create elements
                var element = new Element()
                {
                    ID = iElem + 1,
                    ElementType = beam
                };

                // Add nodes to the created element
                element.AddNode(model.NodesDictionary[iNode]);
                element.AddNode(model.NodesDictionary[iNode + 1]);

                var a = beam.StiffnessMatrix(element);
                var b = beam.MassMatrix(element);

                // Add beam element to the element and subdomains dictionary of the model
                model.ElementsDictionary.Add(element.ID, element);
                model.SubdomainsDictionary[1].ElementsDictionary.Add(element.ID, element);
                iNode++;
            }

            // Add nodal load values at the top nodes of the model
            double totalDuration = 0.27;
            double timeStepDuration = 0.01;
            //var timeFunction = new RampTemporalFunction(nodalLoad, totalDuration, timeStepDuration, 1.0 * totalDuration);
            //model.TimeDependentNodalLoads.Add(new GeneralDynamicNodalLoad(model.NodesDictionary[monitorNode], DOFType.Y, timeFunction));
            model.Loads.Add(new Load() { Amount = nodalLoad, Node = model.NodesDictionary[monitorNode], DOF = DOFType .Y });

            // Needed in order to make all the required data structures
            model.ConnectDataStructures();

            // Choose linear equation system solver
            var linearSystems = new Dictionary<int, ILinearSystem>();
            linearSystems[1] = new SkylineLinearSystem(1, model.Subdomains[0].Forces);
            SolverSkyline solver = new SolverSkyline(linearSystems[1]);

            // Choose the provider of the problem -> here a structural problem
            ProblemStructural provider = new ProblemStructural(model, linearSystems);

            // Choose child analyzer -> Child: NewtonRaphsonNonLinearAnalyzer
            var linearSystemsArray = new[] { linearSystems[1] };
            var subdomainUpdaters = new[] { new NonLinearSubdomainUpdater(model.Subdomains[0]) };
            var subdomainMappers = new[] { new SubdomainGlobalMapping(model.Subdomains[0]) };
            int increments = 1;
            int totalDOFs = model.TotalDOFs;
            int maximumIteration = 120;
            int iterationStepsForMatrixRebuild = 500;
            NewtonRaphsonNonLinearAnalyzer childAnalyzer = new NewtonRaphsonNonLinearAnalyzer(solver, linearSystemsArray, subdomainUpdaters, subdomainMappers,
                provider, increments, totalDOFs);

            // Choose parent analyzer -> Parent: Static or Dynamic
            //StaticAnalyzer parentAnalyzer = new StaticAnalyzer(provider, childAnalyzer, linearSystems);
            NewmarkDynamicAnalyzer parentAnalyzer = new NewmarkDynamicAnalyzer(provider, childAnalyzer, linearSystems, 0.25, 0.5, timeStepDuration, totalDuration); //GRAB THE MAX.

            parentAnalyzer.BuildMatrices();
            parentAnalyzer.Initialize();
            parentAnalyzer.Solve();

            if ((0.975238095 - linearSystems[1].Solution[7]) / 0.975238095 < 0.01) Console.WriteLine("The test problem works.");
        }
        public static void TestBeam3DElasticGeometryNonlinearNewmarkDynamicAnalysisRAMPExample()
        {
            VectorExtensions.AssignTotalAffinityCount();
            double youngModulus = 210000000.0;
            double poissonRatio = 0.3;
            double nodalLoad = 100.0;
            double area = 0.25 * 0.25;
            double inertiaY = Math.Pow(0.25, 4) / 12;
            double inertiaZ = Math.Pow(0.25, 4) / 12;
            double torsionalInertia = Math.Pow(0.25, 4) / 6;
            double effectiveAreaY = area;
            double effectiveAreaZ = area;
            double density = 7.85;
            int nNodes = 3;
            int nElems = 2;
            int monitorNode = 3;

            // Create new 3D material
            ElasticMaterial3D material = new ElasticMaterial3D
            {
                YoungModulus = youngModulus,
                PoissonRatio = poissonRatio,
            };

            // Node creation
            IList<Node> nodes = new List<Node>();
            Node node1 = new Node { ID = 1, X = 0.0, Y = 0.0, Z = 0.0 };
            Node node2 = new Node { ID = 2, X = 5.0, Y = 0.0, Z = 0.0 };
            Node node3 = new Node { ID = 3, X = 10.0, Y = 0.0, Z = 0.0 };

            nodes.Add(node1);
            nodes.Add(node2);
            nodes.Add(node3);

            // Model creation
            Model model = new Model();

            // Add a single subdomain to the model
            model.SubdomainsDictionary.Add(1, new Subdomain() { ID = 1 });

            // Add nodes to the nodes dictonary of the model
            for (int i = 0; i < nodes.Count; ++i)
            {
                model.NodesDictionary.Add(i + 1, nodes[i]);
            }

            // Constrain bottom nodes of the model
            model.NodesDictionary[1].Constraints.Add(DOFType.X);
            model.NodesDictionary[1].Constraints.Add(DOFType.Y);
            model.NodesDictionary[1].Constraints.Add(DOFType.Z);
            model.NodesDictionary[1].Constraints.Add(DOFType.RotX);
            model.NodesDictionary[1].Constraints.Add(DOFType.RotY);
            model.NodesDictionary[1].Constraints.Add(DOFType.RotZ);
            // Generate elements of the structure
            int iNode = 1;
            for (int iElem = 0; iElem < nElems; iElem++)
            {
                // element nodes
                IList<Node> elementNodes = new List<Node>();
                elementNodes.Add(model.NodesDictionary[iNode]);
                elementNodes.Add(model.NodesDictionary[iNode + 1]);

                // Create new Beam3D section and element
                var beamSection = new BeamSection3D(area, inertiaY, inertiaZ, torsionalInertia, effectiveAreaY, effectiveAreaZ);
                var beam = new EulerBeam3D(youngModulus, poissonRatio, null, null);
                beam.MomentOfInertiaY = inertiaY;
                beam.MomentOfInertiaZ = inertiaZ;
                beam.MomentOfInertiaPolar = torsionalInertia;
                beam.SectionArea = area;
                beam.Density = density;
                // Create elements
                var element = new Element()
                {
                    ID = iElem + 1,
                    ElementType = beam
                };

                // Add nodes to the created element
                element.AddNode(model.NodesDictionary[iNode]);
                element.AddNode(model.NodesDictionary[iNode + 1]);

                var a = beam.StiffnessMatrix(element);
                var b = beam.MassMatrix(element);

                // Add beam element to the element and subdomains dictionary of the model
                model.ElementsDictionary.Add(element.ID, element);
                model.SubdomainsDictionary[1].ElementsDictionary.Add(element.ID, element);
                iNode++;
            }

            // Add nodal load values at the top nodes of the model
            double totalDuration = 10;
            double timeStepDuration = 0.005;
            var timeFunction = new RampTemporalFunction(nodalLoad, totalDuration, timeStepDuration, 1.0 * totalDuration);
            model.TimeDependentNodalLoads.Add(new GeneralDynamicNodalLoad(model.NodesDictionary[monitorNode], DOFType.Y, timeFunction));
            //model.Loads.Add(new Load() { Amount = nodalLoad, Node = model.NodesDictionary[monitorNode], DOF = DOFType .Y });

            // Needed in order to make all the required data structures
            model.ConnectDataStructures();

            // Choose linear equation system solver
            var linearSystems = new Dictionary<int, ILinearSystem>();
            linearSystems[1] = new SkylineLinearSystem(1, model.Subdomains[0].Forces);
            SolverSkyline solver = new SolverSkyline(linearSystems[1]);

            // Choose the provider of the problem -> here a structural problem
            ProblemStructural provider = new ProblemStructural(model, linearSystems);

            // Choose child analyzer -> Child: NewtonRaphsonNonLinearAnalyzer
            var linearSystemsArray = new[] { linearSystems[1] };
            var subdomainUpdaters = new[] { new NonLinearSubdomainUpdater(model.Subdomains[0]) };
            var subdomainMappers = new[] { new SubdomainGlobalMapping(model.Subdomains[0]) };
            int increments = 1;
            int totalDOFs = model.TotalDOFs;
            int maximumIteration = 120;
            int iterationStepsForMatrixRebuild = 500;
            NewtonRaphsonNonLinearAnalyzer childAnalyzer = new NewtonRaphsonNonLinearAnalyzer(solver, linearSystemsArray, subdomainUpdaters, subdomainMappers,
                provider, increments, totalDOFs);

            // Choose parent analyzer -> Parent: Static or Dynamic
            //StaticAnalyzer parentAnalyzer = new StaticAnalyzer(provider, childAnalyzer, linearSystems);
            NewmarkDynamicAnalyzer parentAnalyzer = new NewmarkDynamicAnalyzer(provider, childAnalyzer, linearSystems, 0.25, 0.5, timeStepDuration, totalDuration); //GRAB THE MAX.

            parentAnalyzer.BuildMatrices();
            parentAnalyzer.Initialize();
            parentAnalyzer.Solve();

            if ((0.975238095/2 - linearSystems[1].Solution[7]) / (0.975238095/2) < 0.01) Console.WriteLine("The test problem works.");
        }
        private static void SolveHexaSoil()
        {
            VectorExtensions.AssignTotalAffinityCount();
            Model model = new Model();
            model.SubdomainsDictionary.Add(1, new Subdomain() { ID = 1 });

            HexaSoil2.MakeHexaSoil(model);
            var nodeid = HexaSoil2.ProvideIdMonitor(model);
            model.ConnectDataStructures();
            var linearSystems = new Dictionary<int, ILinearSystem>();
            linearSystems[1] = new SkylineLinearSystem(1, model.Subdomains[0].Forces);
            SolverSkyline solver = new SolverSkyline(linearSystems[1]);
            ProblemPorous provider = new ProblemPorous(model, linearSystems);
            var linearSystemsArray = new[] { linearSystems[1] };
            var subdomainUpdaters = new[] { new NonLinearSubdomainUpdater(model.Subdomains[0]) };
            var subdomainMappers = new[] { new SubdomainGlobalMapping(model.Subdomains[0]) };
            int increments = 1;
            int totalDOFs = model.TotalDOFs;
            int maximumIteration = 120;
            int iterationStepsForMatrixRebuild = 500;
            NewtonRaphsonNonLinearAnalyzer childAnalyzer = new NewtonRaphsonNonLinearAnalyzer(solver, linearSystemsArray, subdomainUpdaters, subdomainMappers,
            provider, increments, totalDOFs);
            NewmarkDynamicAnalyzer parentAnalyzer = new NewmarkDynamicAnalyzer(provider, childAnalyzer, linearSystems, 0.25, 0.5, 0.1, 2.0);

            parentAnalyzer.BuildMatrices();
            parentAnalyzer.Initialize();
            parentAnalyzer.Solve(); 
        }
        private static void SolveStochasticHexaSoil(int samplenumber, double Stoch1, double Stoch2, double Stoch3)
        {
            VectorExtensions.AssignTotalAffinityCount();
            Model model = new Model();
            model.SubdomainsDictionary.Add(1, new Subdomain() { ID = 1 });

            HexaSoil2.MakeHexaSoil(model, Stoch1, Stoch2, Stoch3);

            model.ConnectDataStructures();
            var linearSystems = new Dictionary<int, ILinearSystem>();
            linearSystems[1] = new SkylineLinearSystem(1, model.Subdomains[0].Forces);
            SolverSkyline solver = new SolverSkyline(linearSystems[1]);
            ProblemPorous provider = new ProblemPorous(model, linearSystems);
            var linearSystemsArray = new[] { linearSystems[1] };
            var subdomainUpdaters = new[] { new NonLinearSubdomainUpdater(model.Subdomains[0]) };
            var subdomainMappers = new[] { new SubdomainGlobalMapping(model.Subdomains[0]) };
            int increments = 1;
            int totalDOFs = model.TotalDOFs;
            int maximumIteration = 120;
            int iterationStepsForMatrixRebuild = 500;
            NewtonRaphsonNonLinearAnalyzer childAnalyzer = new NewtonRaphsonNonLinearAnalyzer(solver, linearSystemsArray, subdomainUpdaters, subdomainMappers,
            provider, increments, totalDOFs);
            childAnalyzer.totnuminc = 4;
            NewmarkDynamicAnalyzer parentAnalyzer = new NewmarkDynamicAnalyzer(provider, childAnalyzer, linearSystems, 0.25, 0.5, 0.1, 0.50);
            var dofid = HexaSoil2.ProvideIdMonitor(model);
            parentAnalyzer.BuildMatrices();
            parentAnalyzer.Initialize();
            parentAnalyzer.Solve();
            dispstoch[samplenumber] = linearSystems[1].Solution[dofid];
            //Hexa8 h1 = (Hexa8)model.Elements[0].ElementType;
            //IFiniteElementMaterial3D[] matGP = new IFiniteElementMaterial3D[6];
            //matGP = h1.materialsAtGaussPoints;
            //var StressesCheck = matGP[7].ConstitutiveMatrix[5,5];
            //stresstoch[samplenumber] = StressesCheck;
        }
        public static void Main(string[] args)
        {
            //SolveHexaSoil();
            //TestBeam3DElasticGeometryNonlinearNewmarkDynamicAnalysisRAMPExample();
            //Solve(); 
            //SolveBuildingInNoSoilSmall();
            //TrussExample.Run();
            //FEM.Cantilever2D.Run();
            //FEM.Cantilever2DPreprocessor.Run();
            //FEM.WallWithOpenings.Run();
            //SeparateCodeCheckingClass.Check06();
            DateTime begin = DateTime.Now;
            readData("input1.txt", out Stoch1);
            readData("input2.txt", out Stoch2);
            readData("input3.txt", out Stoch3);
            for (int index = 0; index < montecarlosim; index++)
            {
                //SolveStochasticHexaSoil(index, Stoch1[index], Stoch2[index], readMatrixDataPartially(Stoch3, index, index, 0, 7), readMatrixDataPartially(Stoch3, 0, 7, 8, 8));
                SolveStochasticHexaSoil(index, Stoch1[index], Stoch2[index], Stoch3[index]);
                Console.WriteLine(index);
            }
            DateTime end = DateTime.Now;
            writeTime(begin, end);
            writeData(dispstoch, 1);
        }
    }
}
