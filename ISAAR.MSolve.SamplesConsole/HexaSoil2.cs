﻿#region USING
using ISAAR.MSolve.Analyzers;
using ISAAR.MSolve.Discretization;
using ISAAR.MSolve.FEM;
using ISAAR.MSolve.Geometry;
using ISAAR.MSolve.IGA;
using ISAAR.MSolve.Logging;
using ISAAR.MSolve.Materials;
using ISAAR.MSolve.Numerical;
using ISAAR.MSolve.Numerical.LinearAlgebra;
using ISAAR.MSolve.PreProcessor;
using ISAAR.MSolve.PreProcessor.Materials;
using ISAAR.MSolve.Problems;
using ISAAR.MSolve.Solvers.Skyline;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.FEM.Elements;
using ISAAR.MSolve.FEM.Entities;
using ISAAR.MSolve.FEM.Materials;
using ISAAR.MSolve.Materials.Interfaces;
using ISAAR.MSolve.Solvers.Interfaces;
using MGroup.Stochastic;
using MGroup.Stochastic.Structural;
using MGroup.Stochastic.Structural.Example;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.IO;
using System.Threading;
using System.Threading.Tasks;
using ISAAR.MSolve.FEM.Entities.TemporalFunctions;
#endregion

namespace ISAAR.MSolve.SamplesConsole
{
    public class HexaSoil2
    {
        public static void MakeHexaSoil(Model model)
        {
            var startX = 0.0;
            var startY = 0.0;
            var startZ = 0.0;
            var LengthX = 10.0;
            var LengthY = 10.0;
            var LengthZ = 10.0;
            int nodeID = 1;
            var hx = 80.0;
            var hy = 80.0;
            var hz = 20.0;
            var imax = (int)Math.Truncate(hx / LengthX) + 1;
            var jmax = (int)Math.Truncate(hy / LengthY) + 1;
            var kmax = (int)Math.Truncate(hz / LengthZ) + 1;
            for (int l = 0; l < kmax; l++)
            {
                for (int k = 0; k < jmax; k++)
                {
                    for (int j = 0; j < imax; j++)
                    {
                        model.NodesDictionary.Add(nodeID, new Node() { ID = nodeID, X = startX + j * LengthX, Y = startY + k * LengthY, Z = startZ + l * LengthZ });

                        nodeID++;
                    }
                }
            }
            nodeID = 1;
            for (int j = 0; j < jmax; j++)
            {
                for (int k = 0; k < imax; k++)
                {
                    model.NodesDictionary[nodeID].Constraints.Add(DOFType.X);
                    model.NodesDictionary[nodeID].Constraints.Add(DOFType.Y);
                    model.NodesDictionary[nodeID].Constraints.Add(DOFType.Z);
                    nodeID++;
                }
            }
            ElasticMaterial3D material1 = new ElasticMaterial3D()
            {
                YoungModulus = 2.1e5,
                PoissonRatio = 0.35,
            };
            Element e1;
            int IDhelp = 1;
            for (int ii = 1; ii < model.NodesDictionary.Count + 1; ii++)
            {
                var nodecheck = model.NodesDictionary[ii];
                if (nodecheck.X != hx && nodecheck.Y != hy && nodecheck.Z != hz)
                {
                    const int gpNo = 8;
                    var initialStresses = new StressStrainVectorContinuum3D();
                    var element1Nodes = new Node[8];
                    element1Nodes[0] = model.NodesDictionary[ii];
                    element1Nodes[1] = model.NodesDictionary[ii + 1];
                    element1Nodes[2] = model.NodesDictionary[ii + 1 + imax];
                    element1Nodes[3] = model.NodesDictionary[ii + imax];
                    element1Nodes[4] = model.NodesDictionary[ii + jmax * imax];
                    element1Nodes[5] = model.NodesDictionary[ii + jmax * imax + 1];
                    element1Nodes[6] = model.NodesDictionary[ii + jmax * imax + 1 + imax];
                    element1Nodes[7] = model.NodesDictionary[ii + jmax * imax + imax];
                    var nodeCoordinates = new double[8, 3];
                    for (int i = 0; i < 8; i++)
                    {
                        nodeCoordinates[i, 0] = element1Nodes[i].X;
                        nodeCoordinates[i, 1] = element1Nodes[i].Y;
                        nodeCoordinates[i, 2] = element1Nodes[i].Z;
                    }
                    var gaussPointMaterials = new KavvadasClays[8];
                    var young = 2.1e5;
                    var poisson = 0.3;
                    var alpha = 1.0;
                    var ksi = 0.02;
                    var gamma = 10; //effective stress
                    var Htot = hz;
                    for (int i = 0; i < gpNo; i++)
                        gaussPointMaterials[i] = new KavvadasClays(young, poisson, alpha, ksi);
                    var elementType1 = new Hexa8(gaussPointMaterials);
                    var gaussPoints = elementType1.CalculateGaussMatrices(nodeCoordinates);
                    var elementType2 = new Hexa8u8p(gaussPointMaterials); //this because hexa8u8p has not all the fortran that hexa8 has.
                    for (int i = 0; i < gpNo; i++)
                    {
                        var ActualZeta = 0.0;
                        var help = elementType1.CalcH8Shape(gaussPoints[i].Xi, gaussPoints[i].Eta, gaussPoints[i].Zeta);
                        for (int j = 0; j < gpNo; j++)
                        {
                            ActualZeta += help[j] * nodeCoordinates[j, 2];
                        }
                        gaussPointMaterials[i].Zeta = ActualZeta;
                        initialStresses[2] = -gamma * (Htot - gaussPointMaterials[i].Zeta);
                        initialStresses[0] = 0.85 * initialStresses[2];
                        initialStresses[1] = 0.85 * initialStresses[2];
                        initialStresses[3] = 0;
                        initialStresses[4] = 0;
                        initialStresses[5] = 0;
                        gaussPointMaterials[i] = new KavvadasClays(young, poisson, 2, ksi, initialStresses); //NO STOCHASTIC USE.
                    }
                    elementType2 = new Hexa8u8p(gaussPointMaterials);
                    for (int i = 0; i < gpNo; i++)
                    {
                        elementType2.Permeability[i] = Math.Pow(10, -8);
                    }
                    e1 = new Element()
                    {
                        ID = IDhelp
                    };
                    IDhelp++;
                    e1.ElementType = elementType2; //yliko meta diorthosi to material1
                    e1.NodesDictionary.Add(1, model.NodesDictionary[ii]);
                    e1.NodesDictionary.Add(2, model.NodesDictionary[ii + 1]);
                    e1.NodesDictionary.Add(4, model.NodesDictionary[ii + 1 + imax]);
                    e1.NodesDictionary.Add(3, model.NodesDictionary[ii + imax]);
                    e1.NodesDictionary.Add(5, model.NodesDictionary[ii + jmax * imax]);
                    e1.NodesDictionary.Add(6, model.NodesDictionary[ii + jmax * imax + 1]);
                    e1.NodesDictionary.Add(8, model.NodesDictionary[ii + jmax * imax + 1 + imax]);
                    e1.NodesDictionary.Add(7, model.NodesDictionary[ii + jmax * imax + imax]);
                    int subdomainID = 1;
                    //e1.initialForces = e1.ElementType.CalculateForces(e1, initialStresses, initialStresses); no initial forces at this problem we use
                    for (int i = 0; i < gpNo; i++)
                        for (int j = 0; j < 6; j++)
                        {
                            gaussPointMaterials[i].Stresses[j] = gaussPointMaterials[i].Stresses[j] - gaussPointMaterials[i].initialStresses[j];
                        }
                    model.ElementsDictionary.Add(e1.ID, e1);
                    model.SubdomainsDictionary[subdomainID].ElementsDictionary.Add(e1.ID, e1);
                }
            };
            //DOFType doftype1;
            //doftype1 = new DOFType();
            #region initalloads
            //Load loadinitialx = new Load();
            //Load loadinitialy = new Load();
            //Load loadinitialz = new Load();
            //foreach (Node nodecheck in model.NodesDictionary.Values)
            //{
            //    loadinitialx = new Load()
            //    {
            //        Node = nodecheck,
            //        //DOF = doftype1,
            //        DOF = DOFType.X
            //    };
            //    foreach (Element elementcheck in model.ElementsDictionary.Values)
            //    {
            //        var bool1 = elementcheck.NodesDictionary.ContainsValue(nodecheck);
            //        if (bool1 == true)
            //        {
            //        var help1 = 0;
            //        var help = 0;
            //        foreach(Node nodeelement in elementcheck.NodesDictionary.Values)
            //        {
            //            if (nodeelement==nodecheck)
            //            {
            //                    help = help1;
            //            }
            //                help1 += 1;
            //        }
            //            loadinitialx.Amount += -elementcheck.initialForces[3 * help];
            //        }
            //    }
            //    model.Loads.Add(loadinitialx);
            //}
            //foreach (Node nodecheck in model.NodesDictionary.Values)
            //{
            //    loadinitialy = new Load()
            //    {
            //        Node = nodecheck,
            //        //DOF = doftype1,
            //        DOF = DOFType.Y
            //    };
            //    foreach (Element elementcheck in model.ElementsDictionary.Values)
            //    {
            //        var bool1 = elementcheck.NodesDictionary.ContainsValue(nodecheck);
            //        if (bool1 == true)
            //        {
            //            var help1 = 0;
            //            var help = 0;
            //            foreach (Node nodeelement in elementcheck.NodesDictionary.Values)
            //            {
            //                if (nodeelement == nodecheck)
            //                {
            //                    help = help1;
            //                }
            //                help1 += 1;
            //            }
            //            loadinitialy.Amount += -elementcheck.initialForces[3*help+1];
            //        }
            //    }
            //    model.Loads.Add(loadinitialy);
            //}
            //foreach (Node nodecheck in model.NodesDictionary.Values)
            //{
            //    loadinitialz = new Load()
            //    {
            //        Node = nodecheck,
            //        //DOF = doftype1,
            //        DOF = DOFType.Z
            //    };
            //    foreach (Element elementcheck in model.ElementsDictionary.Values)
            //    {
            //        var Pa = -5.0;
            //        var P2a = -10.0/2;
            //        var P4a = -20.0/4;
            //        var bool1 = elementcheck.NodesDictionary.ContainsValue(nodecheck);
            //        var bool2 = nodecheck.Z == hz;
            //        var bool3 = (nodecheck.X == 0 || nodecheck.X == hx) && (nodecheck.Y == 0 || nodecheck.Y == hy);
            //        var bool4 = nodecheck.X != 0 && nodecheck.X != hx && nodecheck.Y != 0 && nodecheck.Y != hy;
            //        if (bool1 == true)
            //        {
            //            var help1 = 0;
            //            var help = 0;
            //            foreach (Node nodeelement in elementcheck.NodesDictionary.Values)
            //            {
            //                if (nodeelement == nodecheck)
            //                {
            //                    help = help1;
            //                }
            //                help1 += 1;
            //            }
            //            if (bool2 == true)
            //            {
            //                if (bool3 == true)
            //                {
            //                    loadinitialz.Amount += -Pa - elementcheck.initialForces[3 * help+2];
            //                }
            //                else if (bool4 == true)
            //                {
            //                    loadinitialz.Amount += -P4a - elementcheck.initialForces[3 * help+2];
            //                }
            //                else
            //                {
            //                    loadinitialz.Amount += -P2a - elementcheck.initialForces[3 * help+2];
            //                }
            //            }
            //            else
            //            {
            //                loadinitialz.Amount += -elementcheck.initialForces[3 * help+2];
            //            }
            //        }
            //    }
            //    model.Loads.Add(loadinitialz);
            //}
            #endregion
            double nodalLoad = 0.0;
            double totalDuration = 2.0;
            double timeStepDuration = 0.1;
            double constantsegmentdurationratio = 0.25;
            GeneralDynamicNodalLoad loadinitialz;
            foreach (Node nodecheck in model.NodesDictionary.Values)
            {
                nodalLoad = 0.0;
                foreach (Element elementcheck in model.ElementsDictionary.Values)
                {
                    var Pa = -3750.0;
                    var P2a = -7500.0 / 2;
                    var P4a = -15000.0 / 4;
                    var bool1 = elementcheck.NodesDictionary.ContainsValue(nodecheck);
                    var bool2 = nodecheck.Z == hz;
                    var bool3 = (nodecheck.X == 0 || nodecheck.X == hx) && (nodecheck.Y == 0 || nodecheck.Y == hy);
                    var bool4 = nodecheck.X != 0 && nodecheck.X != hx && nodecheck.Y != 0 && nodecheck.Y != hy;
                    if (bool1 == true)
                    {
                        var help1 = 0;
                        var help = 0;
                        foreach (Node nodeelement in elementcheck.NodesDictionary.Values)
                        {
                            if (nodeelement == nodecheck)
                            {
                                help = help1;
                            }
                            help1 += 1;
                        }
                        if (bool2 == true)
                        {
                            if (bool3 == true)
                            {
                                nodalLoad += Pa;
                            }
                            else if (bool4 == true)
                            {
                                nodalLoad += P4a;
                            }
                            else
                            {
                                nodalLoad += P2a;
                            }
                        }
                        else
                        {
                            nodalLoad += 0;
                        }
                    }
                }
                if (nodalLoad != 0)
                {
                    var timeFunction = new RampTemporalFunction(nodalLoad, totalDuration, timeStepDuration, constantsegmentdurationratio * totalDuration);
                    loadinitialz = new GeneralDynamicNodalLoad(nodecheck, DOFType.Z, timeFunction);
                    model.TimeDependentNodalLoads.Add(loadinitialz);
                }
            }
        }
        public static void MakeHexaSoil(Model model,double Stoch1,double Stoch2,double[] Stoch3, double[] omega)
        {
            // xreiazetai na rythmizei kaneis ta megethi me auto to configuration
            var startX = 0.0;
            var startY = 0.0;
            var startZ = 0.0;
            var LengthX = 10.0;
            var LengthY = 10.0;
            var LengthZ = 10.0;
            int nodeID = 1;
            var hx = 80.0;
            var hy = 80.0;
            var hz = 20.0;
            var imax = (int)Math.Truncate(hx / LengthX) + 1;
            var jmax = (int)Math.Truncate(hy / LengthY) + 1;
            var kmax = (int)Math.Truncate(hz / LengthZ) + 1;
            for (int l = 0; l < kmax; l++)
            {
                for (int k = 0; k < jmax; k++)
                {
                    for (int j = 0; j < imax; j++)
                    {
                        model.NodesDictionary.Add(nodeID, new Node() { ID = nodeID, X = startX + j * LengthX, Y = startY + k * LengthY, Z = startZ + l * LengthZ });

                        nodeID++;
                    }
                }
            }
            nodeID = 1;
            for (int j = 0; j < jmax; j++)
            {
                for (int k = 0; k < imax; k++)
                {
                    model.NodesDictionary[nodeID].Constraints.Add(DOFType.X);
                    model.NodesDictionary[nodeID].Constraints.Add(DOFType.Y);
                    model.NodesDictionary[nodeID].Constraints.Add(DOFType.Z);
                    nodeID++;
                }
            }
            ElasticMaterial3D material1 = new ElasticMaterial3D()
            {
                YoungModulus = 2.1e5,
                PoissonRatio = 0.35,
            };
            Element e1;
            int IDhelp = 1;
            for (int ii = 1; ii < model.NodesDictionary.Count + 1; ii++)
            {
                var nodecheck = model.NodesDictionary[ii];
                if (nodecheck.X != hx && nodecheck.Y != hy && nodecheck.Z != hz)
                {
                    const int gpNo = 8;
                    var initialStresses = new StressStrainVectorContinuum3D();
                    var element1Nodes = new Node[8];
                    element1Nodes[0] = model.NodesDictionary[ii];
                    element1Nodes[1] = model.NodesDictionary[ii + 1];
                    element1Nodes[2] = model.NodesDictionary[ii + 1 + imax];
                    element1Nodes[3] = model.NodesDictionary[ii + imax];
                    element1Nodes[4] = model.NodesDictionary[ii + jmax * imax];
                    element1Nodes[5] = model.NodesDictionary[ii + jmax * imax + 1];
                    element1Nodes[6] = model.NodesDictionary[ii + jmax * imax + 1 + imax];
                    element1Nodes[7] = model.NodesDictionary[ii + jmax * imax + imax];
                    var nodeCoordinates = new double[8, 3];
                    for (int i = 0; i < 8; i++)
                    {
                        nodeCoordinates[i, 0] = element1Nodes[i].X;
                        nodeCoordinates[i, 1] = element1Nodes[i].Y;
                        nodeCoordinates[i, 2] = element1Nodes[i].Z;
                    }
                    var gaussPointMaterials = new KavvadasClays[8];
                    var young = 2.1e5;
                    var poisson = 0.3;
                    var alpha = 1.0;
                    var ksi = 0.02;
                    var gamma = 10; //effective stress
                    var Htot = hz;
                    for (int i = 0; i < gpNo; i++)
                        gaussPointMaterials[i] = new KavvadasClays(young, poisson, alpha, ksi);
                    var elementType1 = new Hexa8(gaussPointMaterials);
                    var gaussPoints = elementType1.CalculateGaussMatrices(nodeCoordinates);
                    var elementType2 = new Hexa8u8p(gaussPointMaterials); //this because hexa8u8p has not all the fortran that hexa8 has.
                    for (int i = 0; i < gpNo; i++)
                    {
                        var ActualZeta = 0.0;
                        var help = elementType1.CalcH8Shape(gaussPoints[i].Xi, gaussPoints[i].Eta, gaussPoints[i].Zeta);
                        for (int j = 0; j < gpNo; j++)
                        {
                            ActualZeta += help[j] * nodeCoordinates[j, 2];
                        }
                        gaussPointMaterials[i].Zeta = ActualZeta;
                        initialStresses[2] = -gamma * (Htot - gaussPointMaterials[i].Zeta);
                        initialStresses[0] = 0.85*initialStresses[2];
                        initialStresses[1] = 0.85*initialStresses[2];
                        initialStresses[3] = 0;
                        initialStresses[4] = 0;
                        initialStresses[5] = 0;
                        //double comp = 0.0;
                        //double csl = 0.0;
                        //for (int j = 0; j < 8; j = j + 2)
                        //{
                        //   comp += Stoch1[j] * Math.Cos(omega[j] * ActualZeta);
                        //   csl  += Stoch2[j] * Math.Cos(omega[j] * ActualZeta);
                        //}
                        //for (int j = 1; j < 8; j = j + 2)
                        //{
                        //    comp += Stoch1[j] * Math.Sin(omega[j] * ActualZeta);
                        //    csl  += Stoch2[j] * Math.Sin(omega[j] * ActualZeta);
                        //}
                        //comp= Math.Abs((comp) * 0.25 * 0.008686 + 0.008686) / 1;
                        //csl = Math.Abs((csl) * 0.25 * 0.733609251 + 0.733609251) / 1;
                        gaussPointMaterials[i] = new KavvadasClays(Stoch1,Stoch2, 1, ksi, initialStresses);
                    }
                    elementType2 = new Hexa8u8p(gaussPointMaterials);
                    for (int i = 0; i < gpNo; i++)
                    {
                        var ActualZeta = 0.0;
                        var help = elementType1.CalcH8Shape(gaussPoints[i].Xi, gaussPoints[i].Eta, gaussPoints[i].Zeta);
                        for (int j = 0; j < gpNo; j++)
                        {
                            ActualZeta += help[j] * nodeCoordinates[j, 2];
                        }
                        for (int j = 0; j < 8; j = j + 2)
                        {
                            elementType2.Permeability[i] += Stoch3[j] * Math.Cos(omega[j] * ActualZeta);
                        }
                        for (int j = 1; j < 8; j = j + 2)
                        {
                            elementType2.Permeability[i] += Stoch3[j] * Math.Sin(omega[j] * ActualZeta);
                        }
                        elementType2.Permeability[i] = Math.Abs((elementType2.Permeability[i]) * 0.25 * Math.Pow(10, -8) + Math.Pow(10, -8)) / 1;
                    }
                    e1 = new Element()
                    {
                        ID = IDhelp
                    };
                    IDhelp++;
                    e1.ElementType = elementType2; //yliko meta diorthosi to material1
                    e1.NodesDictionary.Add(1, model.NodesDictionary[ii]);
                    e1.NodesDictionary.Add(2, model.NodesDictionary[ii + 1]);
                    e1.NodesDictionary.Add(4, model.NodesDictionary[ii + 1 + imax]);
                    e1.NodesDictionary.Add(3, model.NodesDictionary[ii + imax]);
                    e1.NodesDictionary.Add(5, model.NodesDictionary[ii + jmax * imax]);
                    e1.NodesDictionary.Add(6, model.NodesDictionary[ii + jmax * imax + 1]);
                    e1.NodesDictionary.Add(8, model.NodesDictionary[ii + jmax * imax + 1 + imax]);
                    e1.NodesDictionary.Add(7, model.NodesDictionary[ii + jmax * imax + imax]);
                    int subdomainID = 1;
                    //e1.initialForces = e1.ElementType.CalculateForces(e1, initialStresses, initialStresses); no initial forces at this problem we use
                    for (int i = 0; i < gpNo; i++)
                        for (int j = 0; j < 6; j++)
                        {
                            gaussPointMaterials[i].Stresses[j] = gaussPointMaterials[i].Stresses[j] - gaussPointMaterials[i].initialStresses[j];
                        }
                    model.ElementsDictionary.Add(e1.ID, e1);
                    model.SubdomainsDictionary[subdomainID].ElementsDictionary.Add(e1.ID, e1);
                }
            };
            //DOFType doftype1;
            //doftype1 = new DOFType();
            #region initalloads
            //Load loadinitialx = new Load();
            //Load loadinitialy = new Load();
            //Load loadinitialz = new Load();
            //foreach (Node nodecheck in model.NodesDictionary.Values)
            //{
            //    loadinitialx = new Load()
            //    {
            //        Node = nodecheck,
            //        //DOF = doftype1,
            //        DOF = DOFType.X
            //    };
            //    foreach (Element elementcheck in model.ElementsDictionary.Values)
            //    {
            //        var bool1 = elementcheck.NodesDictionary.ContainsValue(nodecheck);
            //        if (bool1 == true)
            //        {
            //        var help1 = 0;
            //        var help = 0;
            //        foreach(Node nodeelement in elementcheck.NodesDictionary.Values)
            //        {
            //            if (nodeelement==nodecheck)
            //            {
            //                    help = help1;
            //            }
            //                help1 += 1;
            //        }
            //            loadinitialx.Amount += -elementcheck.initialForces[3 * help];
            //        }
            //    }
            //    model.Loads.Add(loadinitialx);
            //}
            //foreach (Node nodecheck in model.NodesDictionary.Values)
            //{
            //    loadinitialy = new Load()
            //    {
            //        Node = nodecheck,
            //        //DOF = doftype1,
            //        DOF = DOFType.Y
            //    };
            //    foreach (Element elementcheck in model.ElementsDictionary.Values)
            //    {
            //        var bool1 = elementcheck.NodesDictionary.ContainsValue(nodecheck);
            //        if (bool1 == true)
            //        {
            //            var help1 = 0;
            //            var help = 0;
            //            foreach (Node nodeelement in elementcheck.NodesDictionary.Values)
            //            {
            //                if (nodeelement == nodecheck)
            //                {
            //                    help = help1;
            //                }
            //                help1 += 1;
            //            }
            //            loadinitialy.Amount += -elementcheck.initialForces[3*help+1];
            //        }
            //    }
            //    model.Loads.Add(loadinitialy);
            //}
            //foreach (Node nodecheck in model.NodesDictionary.Values)
            //{
            //    loadinitialz = new Load()
            //    {
            //        Node = nodecheck,
            //        //DOF = doftype1,
            //        DOF = DOFType.Z
            //    };
            //    foreach (Element elementcheck in model.ElementsDictionary.Values)
            //    {
            //        var Pa = -5.0;
            //        var P2a = -10.0/2;
            //        var P4a = -20.0/4;
            //        var bool1 = elementcheck.NodesDictionary.ContainsValue(nodecheck);
            //        var bool2 = nodecheck.Z == hz;
            //        var bool3 = (nodecheck.X == 0 || nodecheck.X == hx) && (nodecheck.Y == 0 || nodecheck.Y == hy);
            //        var bool4 = nodecheck.X != 0 && nodecheck.X != hx && nodecheck.Y != 0 && nodecheck.Y != hy;
            //        if (bool1 == true)
            //        {
            //            var help1 = 0;
            //            var help = 0;
            //            foreach (Node nodeelement in elementcheck.NodesDictionary.Values)
            //            {
            //                if (nodeelement == nodecheck)
            //                {
            //                    help = help1;
            //                }
            //                help1 += 1;
            //            }
            //            if (bool2 == true)
            //            {
            //                if (bool3 == true)
            //                {
            //                    loadinitialz.Amount += -Pa - elementcheck.initialForces[3 * help+2];
            //                }
            //                else if (bool4 == true)
            //                {
            //                    loadinitialz.Amount += -P4a - elementcheck.initialForces[3 * help+2];
            //                }
            //                else
            //                {
            //                    loadinitialz.Amount += -P2a - elementcheck.initialForces[3 * help+2];
            //                }
            //            }
            //            else
            //            {
            //                loadinitialz.Amount += -elementcheck.initialForces[3 * help+2];
            //            }
            //        }
            //    }
            //    model.Loads.Add(loadinitialz);
            //}
            #endregion
            double nodalLoad = 0.0;
            double totalDuration = 0.50;
            double timeStepDuration = 0.1;
            double constantsegmentdurationratio = 0.8;
            GeneralDynamicNodalLoad loadinitialz;
            foreach (Node nodecheck in model.NodesDictionary.Values)
            {
                nodalLoad = 0.0;
                foreach (Element elementcheck in model.ElementsDictionary.Values)
                {
                    var Pa = -375.0;
                    var P2a = -750.0 / 2;
                    var P4a = -1500.0 / 4;
                    var bool1 = elementcheck.NodesDictionary.ContainsValue(nodecheck);
                    var bool2 = nodecheck.Z == hz;
                    var bool3 = (nodecheck.X == 0 || nodecheck.X == hx) && (nodecheck.Y == 0 || nodecheck.Y == hy);
                    var bool4 = nodecheck.X != 0 && nodecheck.X != hx && nodecheck.Y != 0 && nodecheck.Y != hy;
                    if (bool1 == true)
                    {
                        var help1 = 0;
                        var help = 0;
                        foreach (Node nodeelement in elementcheck.NodesDictionary.Values)
                        {
                            if (nodeelement == nodecheck)
                            {
                                help = help1;
                            }
                            help1 += 1;
                        }
                        if (bool2 == true)
                        {
                            if (bool3 == true)
                            {
                                nodalLoad += Pa;
                            }
                            else if (bool4 == true)
                            {
                                nodalLoad += P4a;
                            }
                            else
                            {
                                nodalLoad += P2a;
                            }
                        }
                        else
                        {
                            nodalLoad += 0;
                        }
                    }
                }
                if (nodalLoad != 0)
                {
                    var timeFunction = new RampTemporalFunction(nodalLoad, totalDuration, timeStepDuration, constantsegmentdurationratio * totalDuration);
                    loadinitialz = new GeneralDynamicNodalLoad(nodecheck, DOFType.Z, timeFunction);
                    model.TimeDependentNodalLoads.Add(loadinitialz);
                }
            }

        }
        public static int ProvideIdMonitor(Model model)
        {
            //works only in the case when all the points are tottaly constrained. Note that the enumeration is beginning for the free dofs and
            //counts only for them. Later we may make it more general this method.
            var startX = 0.0;
            var startY = 0.0;
            var startZ = 0.0;
            var LengthX = 10.0;
            var LengthY = 10.0;
            var LengthZ = 10.0;
            var hx = 80.0;
            var hy = 80.0;
            var hz = 20.0;
            var imax = (int)Math.Truncate(hx / LengthX) + 1;
            var jmax = (int)Math.Truncate(hy / LengthY) + 1;
            var kmax = (int)Math.Truncate(hz / LengthZ) + 1;
            var nodeid = 1;
            var mx = hx / 2;
            var my = hy / 2;
            var mz = hz;
            var Monitor = 0;
            for (int l = 0; l < kmax; l++)
            {
                for (int k = 0; k < jmax; k++)
                {
                    for (int j = 0; j < imax; j++)
                    {
                        var X = startX + j * LengthX;
                        var Y = startY + k * LengthY;
                        var Z = startZ + l * LengthZ;
                        var boolx = (mx == X);
                        var booly = (my == Y);
                        var boolz = (mz == Z);
                        var conx = model.NodesDictionary[nodeid].Constraints.Contains(DOFType.X);
                        var cony = model.NodesDictionary[nodeid].Constraints.Contains(DOFType.Y);
                        var conz = model.NodesDictionary[nodeid].Constraints.Contains(DOFType.Z);
                        if (conx == false || cony == false || conz == false)
                        {
                            if (boolx && booly && boolz == true)
                            {
                                Monitor = 4 * nodeid - 2+81;
                            }
                        }
                        if (Z != 0)
                        {
                            nodeid++;
                        }
                    }
                }
            }
            return Monitor;
        }
    }
}

