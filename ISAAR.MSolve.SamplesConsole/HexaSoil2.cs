using ISAAR.MSolve.Analyzers;
using ISAAR.MSolve.Analyzers.NonLinear;
using ISAAR.MSolve.PreProcessor;
using ISAAR.MSolve.Discretization;
using ISAAR.MSolve.FEM;
using ISAAR.MSolve.Discretization.Commons;
using ISAAR.MSolve.Discretization.FreedomDegrees;
using ISAAR.MSolve.Discretization.Integration.Quadratures;
using ISAAR.MSolve.FEM.Elements;
using ISAAR.MSolve.FEM.Entities;
using ISAAR.MSolve.Logging;
using ISAAR.MSolve.Materials;
using ISAAR.MSolve.Problems;
using ISAAR.MSolve.Solvers.Direct;
using ISAAR.MSolve.PreProcessor.Materials;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.IO;
using ISAAR.MSolve.FEM.Entities.TemporalFunctions;

namespace ISAAR.MSolve.SamplesConsole
{
    public class HexaSoil2
    {
        public static double startX = 0.0;
        public static double startY = 0.0;
        public static double startZ = 0.0;
        public static double LengthX = 10.0;
        public static double LengthY = 10.0;
        public static double LengthZ = 10.0;
        int nodeID = 1;
        public static double hx = 80.0;
        public static double hy = 80.0;
        public static double hz = 20.0;
        public static int imax = (int)Math.Truncate(hx / LengthX) + 1;
        public static int jmax = (int)Math.Truncate(hy / LengthY) + 1;
        public static int kmax = (int)Math.Truncate(hz / LengthZ) + 1;
        public static double qfail=0.0;
        public static double ufail = 0.0;
        public static void MakeHexaSoil(Model model, double Stoch1, double Stoch2, double[] Stoch3,double[] omega,double lambda)
        {
            // xreiazetai na rythmizei kaneis ta megethi me auto to configuration
            int nodeID = 1;
            for (int l = 0; l < kmax; l++)
            {
                for (int k = 0; k < jmax; k++)
                {
                    for (int j = 0; j < imax; j++)
                    {
                        model.NodesDictionary.Add(nodeID,new Node(nodeID,startX + j * LengthX, startY + k * LengthY, startZ + l * LengthZ));
                        nodeID++;
                    }
                }
            }
            nodeID = 1;
            for (int j = 0; j < jmax; j++)
            {
                for (int k = 0; k < imax; k++)
                {
                    model.NodesDictionary[nodeID].Constraints.Add(new Constraint { DOF = StructuralDof.TranslationX });
                    model.NodesDictionary[nodeID].Constraints.Add(new Constraint { DOF = StructuralDof.TranslationY });
                    model.NodesDictionary[nodeID].Constraints.Add(new Constraint { DOF = StructuralDof.TranslationZ });
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
                    var initialStresses = new double[6];
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
                    var gamma = 20;  //effective stress
                    var Htot = hz;
                    for (int i = 0; i < gpNo; i++)
                        gaussPointMaterials[i] = new KavvadasClays(young, poisson, alpha, ksi);
                    var elementType1 = new Hexa8Fixed(gaussPointMaterials);
                    var gaussPoints = elementType1.CalculateGaussMatrices(nodeCoordinates);
                    var elementType2 = new Hexa8Fixed(gaussPointMaterials); //this because hexa8u8p has not all the fortran that hexa8 has.
                    for (int i = 0; i < gpNo; i++)
                    {
                        var ActualZeta = 0.0;
                        var ActualXi = 0.0;
                        var ActualPsi = 0.0;
                        double[] Coord = new double[3];
                        var help = elementType1.CalcH8Shape(gaussPoints[i].Xi, gaussPoints[i].Eta, gaussPoints[i].Zeta);
                        for (int j = 0; j < gpNo; j++)
                        {
                            ActualXi += help[j] * nodeCoordinates[j, 0];
                            ActualPsi += help[j] * nodeCoordinates[j, 1];
                            ActualZeta += help[j] * nodeCoordinates[j, 2];
                        }
                        gaussPointMaterials[i].Zeta = ActualZeta;
                        Coord[0] = ActualXi;
                        Coord[1] = ActualPsi;
                        Coord[2] = ActualZeta;
                        initialStresses[2] = -gamma * (Htot - gaussPointMaterials[i].Zeta);
                        initialStresses[0] = 0.85 * initialStresses[2];
                        initialStresses[1] = 0.85 * initialStresses[2];
                        initialStresses[3] = 0;
                        initialStresses[4] = 0;
                        initialStresses[5] = 0;
                        gaussPointMaterials[i] = new KavvadasClays(Stoch1, Stoch2, 1, ksi, initialStresses,Htot,Coord);
                    }
                    //elementType1 = new Hexa8(gaussPointMaterials);
                    elementType2 = new Hexa8Fixed(gaussPointMaterials);
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
                    model.SubdomainsDictionary[subdomainID].Elements.Add(e1);
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
            foreach (Node nodecheck in model.NodesDictionary.Values)
            {
                Load loadinitialz = new Load();
                nodalLoad = 0.0;
                foreach (Element elementcheck in model.ElementsDictionary.Values)
                {
                    var Pa = lambda*-3750.0;
                    var P2a = lambda*-7500.0 / 2;
                    var P4a = lambda *-15000.0 / 4;
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
                    //var timeFunction = new RampTemporalFunction(nodalLoad, totalDuration, timeStepDuration, constantsegmentdurationratio * totalDuration);
                    //loadinitialz = new GeneralDynamicNodalLoad(nodecheck, StructuralDof.TranslationZ, timeFunction);
                    //model.TimeDependentNodalLoads.Add(loadinitialz);
                    loadinitialz.Amount = nodalLoad;
                    model.Loads.Add(loadinitialz);
                }
            }
        }
        public static int ProvideIdMonitor(Model model)
        {
            //works only in the case when all the points are tottaly constrained. Note that the enumeration is beginning for the free dofs and
            //counts only for them. Later we may make it more general this method.
            var nodeid = 1;
            var mx = hx/2;
            var my = hy/2;
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
                        bool conx = false;
                        bool cony = false;
                        bool conz = false;
                        //ATTENTION. THIS WORKS IF AND ONLY IF THE CONSTRAINTS ARE INPUTED IN THE MODEL LIKE THIS CONSTRAINTS[0]==STRUCTURAL.TRANSLATIONX,CONSTRAINTS[1]==STRUCTURAL.TRANSLATIONY,CONSTRAINTS[2]==STRUCTURAL.TRANSLATIONZ 
                        if (model.NodesDictionary[nodeid].Constraints.Count == 3)
                        {
                            conx = model.NodesDictionary[nodeid].Constraints[0].DOF.Equals(StructuralDof.TranslationX);
                            cony = model.NodesDictionary[nodeid].Constraints[1].DOF.Equals(StructuralDof.TranslationY);
                            conz = model.NodesDictionary[nodeid].Constraints[2].DOF.Equals(StructuralDof.TranslationZ);
                            //END OF ATTENTION
                        }
                        if (conx == false || cony == false || conz == false)
                        {
                            if (boolx && booly && boolz == true)
                            {
                                Monitor = 4 * nodeid - 2 + ((int)(hx / LengthX) + 1) * ((int)(hy / LengthY) + 1);
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

