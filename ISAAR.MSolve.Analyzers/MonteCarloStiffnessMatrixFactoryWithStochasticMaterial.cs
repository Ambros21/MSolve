﻿using System;
using System.Collections.Generic;
using ISAAR.MSolve.Logging.Interfaces;
using ISAAR.MSolve.Analyzers.Interfaces;
using ISAAR.MSolve.Solvers.Interfaces;
using ISAAR.MSolve.PreProcessor.Interfaces;
using ISAAR.MSolve.PreProcessor;
using ISAAR.MSolve.PreProcessor.Stochastic;
using ISAAR.MSolve.Matrices.Interfaces;
using Troschuetz.Random.Distributions.Continuous;

namespace ISAAR.MSolve.Analyzers
{
    public class MonteCarloAnalyzerStiffnessMatrixFactoryWithStochasticMaterial : IAnalyzer
    {
        private int currentSimulation = -1;
        private readonly int expansionOrder;
        private readonly int simulations;
        private readonly IDictionary<int, ISolverSubdomain> subdomains;
        //private readonly IDictionary<int, IMatrix2D<double>> matrices;
        private readonly IDictionary<int, IMatrix2D<double>>[] matrices;
        private readonly Model model;
        private readonly Dictionary<int, IAnalyzerLog[]> logs = new Dictionary<int, IAnalyzerLog[]>();
        private readonly IAnalyzerProvider provider;
        private readonly double[][] randomNumbers;
        //private readonly double[] stochasticDomain;
        private IAnalyzer childAnalyzer;
        private IAnalyzer parentAnalyzer = null;
        private readonly IStochasticMaterialCoefficientsProvider coefficientsProvider;

        public MonteCarloAnalyzerStiffnessMatrixFactoryWithStochasticMaterial(Model model, IAnalyzerProvider provider, IAnalyzer embeddedAnalyzer, IDictionary<int, ISolverSubdomain> subdomains, GaussianFileStochasticCoefficientsProvider coefficientsProvider, int expansionOrder, int simulations)
        {
            this.childAnalyzer = embeddedAnalyzer;
            this.provider = provider;
            this.model = model;
            this.subdomains = subdomains;
            this.expansionOrder = expansionOrder;
            this.simulations = simulations;
            this.childAnalyzer.ParentAnalyzer = this;
            //this.matrices = new Dictionary<int, IMatrix2D<double>>(subdomains.Count);
            this.matrices = new Dictionary<int, IMatrix2D<double>>[expansionOrder + 1];
            this.randomNumbers = new double[simulations][];
            this.coefficientsProvider = coefficientsProvider;
            //this.stochasticDomain = stochasticDomain;

            NormalDistribution n = new NormalDistribution();
            n.Mu = 0;
            n.Sigma = 1;
            string[] randoms = new string[simulations];
            for (int i = 0; i < simulations; i++)
            {
                randomNumbers[i] = new double[expansionOrder];
                for (int j = 0; j < expansionOrder; j++)
                    randomNumbers[i][j] = n.NextDouble();
                randoms[i] = randomNumbers[i][0].ToString();
            }
            //File.WriteAllLines(String.Format(@"randoms.txt", expansionOrder), randoms);
        }
        public void Solve(int h, int j)
        {

        }
        #region IAnalyzer Members

        public Dictionary<int, IAnalyzerLog[]> Logs { get { return logs; } }
        public IAnalyzer ParentAnalyzer
        {
            get { return parentAnalyzer; }
            set { parentAnalyzer = value; }
        }

        public IAnalyzer ChildAnalyzer
        {
            get { return childAnalyzer; }
            set { childAnalyzer = value; }
        }

        public void BuildMatrices()
        {
            if (childAnalyzer == null) throw new InvalidOperationException("Monte Carlo analyzer must contain an embedded analyzer.");
            if (currentSimulation < 0) return;

            provider.Reset();
            childAnalyzer.BuildMatrices();
        }

        public void Initialize()
        {
            if (childAnalyzer == null) throw new InvalidOperationException("Monte Carlo analyzer must contain an embedded analyzer.");

            for (int i = 0; i < simulations; i++)
            {
                currentSimulation = i;
                coefficientsProvider.RandomVariables = randomNumbers[currentSimulation];
                BuildMatrices();
            }
        }

        public void Solve()
        {
            throw new InvalidOperationException("Monte Carlo stifness factory cannot solve.");
        }

        #endregion
    }
}
