﻿using ISAAR.MSolve.Analyzers;
using ISAAR.MSolve.Logging;
using ISAAR.MSolve.Matrices;
using ISAAR.MSolve.PreProcessor;
using ISAAR.MSolve.Problems;
using ISAAR.MSolve.Solvers.Skyline;
using ISAAR.MSolve.PreProcessor.Elements;
using ISAAR.MSolve.PreProcessor.Materials;
using ISAAR.MSolve.PreProcessor.Interfaces;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.IO;
using System.Threading;
using System.Threading.Tasks;
namespace ISAAR.MSolve.SamplesConsole
{
    class Program2
    {
        #region HexaProgramCode
        public static double[] Stoch1;
        public static double[] Stoch2;
        public static double[,] Stoch3;
        public static int montecarlosim = 1;
        public static double[] increments;
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
        private static void SolveHexaSoil()
        {
            VectorExtensions.AssignTotalAffinityCount();
            Model model = new Model();
            model.SubdomainsDictionary.Add(1, new Subdomain() { ID = 1 });

            HexaSoil2.MakeHexaSoil(model);

            model.ConnectDataStructures();

            SolverSkyline solver = new SolverSkyline(model);

            ProblemStructural provider = new ProblemStructural(model, solver.SubdomainsDictionary);
            NonLinearAnalyzerNewtonRaphsonNew analyzer = NonLinearAnalyzerNewtonRaphsonNew.NonLinearAnalyzerWithFixedLoadIncrements(solver, solver.SubdomainsDictionary, provider, 100, model.TotalDOFs);
            StaticAnalyzer parentAnalyzer = new StaticAnalyzer(provider, analyzer, solver.SubdomainsDictionary);
            analyzer.dofid = 14;
            parentAnalyzer.BuildMatrices();
            parentAnalyzer.Initialize();
            parentAnalyzer.Solve();
            writeData(analyzer.displacements, 0);
        }
        private static void SolveStochasticHexaSoil(int samplenumber, double Stoch1, double Stoch2, double[] Stoch3, double[] omega)
        {
            VectorExtensions.AssignTotalAffinityCount();
            Model model = new Model();
            model.SubdomainsDictionary.Add(1, new Subdomain() { ID = 1 });

            HexaSoil2.MakeHexaSoil(model, Stoch1, Stoch2, Stoch3, omega);

            model.ConnectDataStructures();


            SolverSkyline solver = new SolverSkyline(model);
            //ProblemStructural provider = new ProblemStructural(model, solver.SubdomainsDictionary);
            ProblemPorous provider = new ProblemPorous(model, solver.SubdomainsDictionary);
            NonLinearAnalyzerNewtonRaphsonNew analyzer = NonLinearAnalyzerNewtonRaphsonNew.nonLinearAnalyzerWithPrescribedIncrements(solver, solver.SubdomainsDictionary, provider, model.TotalDOFs, increments);
            //StaticAnalyzer parentAnalyzer = new StaticAnalyzer(provider, analyzer, solver.SubdomainsDictionary);

            NewmarkDynamicAnalyzer parentAnalyzer = new NewmarkDynamicAnalyzer(provider, analyzer, solver.SubdomainsDictionary, 0.25, 0.5, 0.1, 0.51);
            analyzer.dofid = HexaSoil2.ProvideIdMonitor(model);
            parentAnalyzer.BuildMatrices();
            parentAnalyzer.Initialize();
            parentAnalyzer.Solve();
            dispstoch[samplenumber] = analyzer.displacements[4];
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
            readData("input1.txt", out Stoch1);
            readData("input2.txt", out Stoch2);
            readMatrixData("input3.txt", out Stoch3);
            readData("timefun.txt", out increments);
            for (int index = 0; index < montecarlosim; index++)
            {
                SolveStochasticHexaSoil(index, Stoch1[index], Stoch2[index], readMatrixDataPartially(Stoch3, index, index, 0, 7), readMatrixDataPartially(Stoch3, 0, 7, 8, 8));
            }
            //for (int i = 0; i < 1; i++)
            //{
            //    SolveStochasticHexaSoil(1, Stoch1[1], Stoch2[1]);
            //}
            //Parallel.For (0, montecarlosim,
            //      index =>
            //      {
            //          SolveStochasticHexaSoil(index, Stoch1[index+250], Stoch2[index+250], readMatrixDataPartially(Stoch3,index+250,index+250,0,7),readMatrixDataPartially(Stoch3,0,7,8,8));
            //      }) ;
            DateTime end = DateTime.Now;
            writeTime(begin, end);
            writeData(dispstoch, 1);
            //writeData(stresstoch, 1);
        }
        #endregion
        #region FiberBeamCode
        //private static void SolveFibers()
        //{
        //    VectorExtensions.AssignTotalAffinityCount();
        //    Model model = new Model();
        //    model.SubdomainsDictionary.Add(1, new Subdomain() { ID = 1 });
        //    FiberBeam.MakeFiberBeamModel(model);
        //    model.ConnectDataStructures();

        //    SolverSkyline solver = new SolverSkyline(model);

        //    ProblemStructural provider = new ProblemStructural(model, solver.SubdomainsDictionary);
        //    NonLinearAnalyzerNewtonRaphsonNew analyzer = NonLinearAnalyzerNewtonRaphsonNew.NonLinearAnalyzerWithFixedLoadIncrements(solver, solver.SubdomainsDictionary, provider, 20, model.TotalDOFs);
        //    StaticAnalyzer parentAnalyzer = new StaticAnalyzer(provider, analyzer, solver.SubdomainsDictionary);
        //    analyzer.dofid = 4;
        //    parentAnalyzer.BuildMatrices();
        //    parentAnalyzer.Initialize();
        //    parentAnalyzer.Solve();
        //    writeData(analyzer.displacements, 1);
        //}
        //static void Main(string[] args)
        //{
        //    SolveFibers();
        //}
        #endregion
    }
}
