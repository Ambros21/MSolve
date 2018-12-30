using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using ISAAR.MSolve.PreProcessor.Interfaces;
using ISAAR.MSolve.Matrices.Interfaces;
using ISAAR.MSolve.Matrices;
using ISAAR.MSolve.PreProcessor.Materials;
using ISAAR.MSolve.PreProcessor.Elements.SupportiveClasses;
using static ISAAR.MSolve.PreProcessor.Elements.SupportiveClasses.GaussLobatto;


namespace ISAAR.MSolve.PreProcessor.Elements
{
    public class FiberBeam3D : IFiberFiniteElement
    {
        //private SymmetricMatrix2D<double>[] sectionStiffnesses; // instead of storing them in the Gauss Point
        //modified by Ambrosios Savvides. The calculate rot transformation stuff is not checked thus the connectivity between rotations and nearby
        // displacements needs to be checked and added.It is the Euler Bernoulli Beam with Torsion Torque and Shear Forces to be calculated via equillibrium.
        // Also here the "fibers" are considered points in the rectangles widths and in the center.
        private static readonly DOFType[] noDOFTypes = new DOFType[0];
        private static readonly DOFType[] nodalDOFTypes = new DOFType[6] { DOFType.X, DOFType.Y, DOFType.Z, DOFType.RotX, DOFType.RotY, DOFType.RotZ };
        //private static readonly DOFType[][] dofs = new DOFType[][] { nodalDOFTypes, nodalDOFTypes, noDOFTypes };
        private static readonly DOFType[][] dofs = new DOFType[][] { nodalDOFTypes, nodalDOFTypes };
        private GaussLobattoPoint1D[] points = new GaussLobattoPoint1D[4];
        private int noOfDOFs = 12;
        private Matrix2D<double> elementStiffnessMatrix = new Matrix2D<double>(12, 12);
        private Matrix2D<double> currentStiffnessMatrix = new Matrix2D<double>(6, 6);
        private Matrix2D<double> currentFlexibilityMatrix = new Matrix2D<double>(6, 6);
        private bool transformationMatricesInitialized = false;
        private bool fibersInitialized = false;
        private IFiberFiniteElementMaterial material;
        private readonly List<List<Fiber>> sections;
        private List<Fiber> fibers;
        private readonly double b, h;
        private double length, sectionArea, iyy, izz, iyz, it;
        private Matrix2D<double> rotTransformation;
        private readonly Matrix2D<double> t04 = new Matrix2D<double>(12, 12);
        private readonly Matrix2D<double> an = new Matrix2D<double>(6, 12);
        private readonly Matrix2D<double> b04 = new Matrix2D<double>(3, 6);
        private readonly Node[][] rotNodes;
        private IFiniteElementDOFEnumerator dofEnumerator = new GenericDOFEnumerator();
        private bool firstiteration = false;
        private Vector<double> currentInternalForces = new Vector<double>(12);
        private Vector<double> currentInternalForcesBalanced = new Vector<double>(12);
        private double[] dcurrentInternalForces = new double[12];
        private double[] dcurrentInternalForcesBalanced = new double[12];
        private double TorsionS;
        public FiberBeam3D(IFiberFiniteElementMaterial material, int noOfFibers, double b, double h)
        {
            this.b = b;
            this.h = h;
            this.material = material;
            points= GaussLobatto.GetGaussLobattoPoints();
            sections = new List<List<Fiber>>(1); // The row of the sections is the row given in Gauss Lobatto class.
            InitializeFibers(noOfFibers);
        }

        public FiberBeam3D(IFiberFiniteElementMaterial material1, IFiniteElementDOFEnumerator dofEnumerator, int noOfFibers, double b, double h,
            Node[] rot1Nodes, Node[] rot2Nodes) : this(material1, noOfFibers, b, h)
        {
            this.dofEnumerator = dofEnumerator;
        }

        public IFiniteElementDOFEnumerator DOFEnumerator
        {
            get { return dofEnumerator; }
            set { dofEnumerator = value; }
        }

        public double Density { get; set; }

        public double SectionArea { get { return sectionArea; } }
        public double Iyy { get { return iyy; } }
        public double Izz { get { return izz; } }
        public double Iyz { get { return iyz; } }
        public double It { get { return it; } }

        private void InitializeFibers(int noOfFibers)
        {
                int rows = 51; //user defined. Always choose odd number in order for the Simpson Integration to be defined.
                int cols = 11; // user defined. Always choose odd number in order for the Simpson Integration to be defined.
                noOfFibers = 4 * rows * cols;
                fibers = new List<Fiber>(noOfFibers);
                double fiberB = b / (cols - 1);
                double fiberH = h / (rows - 1);
                int pos = 0;
                for (int hhh = 0; hhh < 4; hhh++)
                {
                    for (int i = 0; i < rows; i++)
                        for (int j = 0; j < cols; j++)
                        {
                            double x = -b * 0.5 + j * fiberB;
                            double y = -h * 0.5 + i * fiberH;
                            Fiber fiber = new Fiber(fiberB, fiberH, x, y) { Material = material.FiberMaterials[pos] };
                            fibers.Add(fiber);
                            pos++;
                        }
                }
                sections.Add(fibers);
                sectionArea = b * h;
                iyy = b * Math.Pow(h, 3) / 12;
                izz = h * Math.Pow(b, 3) / 12;
                iyz = 0; //As a rectangular shape of section. 
                it = iyy + izz;
            foreach (GaussLobattoPoint1D point in this.points)
            {
                point.sectionforces = new Vector<double>(3);
                point.sectionforcesbalanced = new Vector<double>(3);
                point.residualdef = new double[3];
                point.deformations = new Vector<double>(3);
                point.deformationsbalanced = new Vector<double>(3);
                point.sectionflexibilityMatrix[0, 0] = 1.0 / (Material.YoungModulus * SectionArea);
                point.sectionflexibilityMatrix[1, 1] = 1.0 / (Material.YoungModulus * Iyy);
                point.sectionflexibilityMatrix[2, 2] = 1.0 / (Material.YoungModulus * Izz);
            }
            fibersInitialized = true;
        }
        private void CalculateTransformationMatrices(Element element)
        {
            // Points A,B,C,D are respectively the j and k nodal points of the element and two points that define the local y and z axes of the element. 
            // Theese points C and D shall be input from the user.
            double Ax = element.Nodes[0].X;
            double Ay = element.Nodes[0].Y;
            double Az = element.Nodes[0].Z;
            double Bx = element.Nodes[1].X;
            double By = element.Nodes[1].Y;
            double Bz = element.Nodes[1].Z;
            double Cx = element.Nodes[0].X;
            double Cy = element.Nodes[0].Y+b/2;
            double Cz = element.Nodes[0].Z;
            double Dx = element.Nodes[0].X;
            double Dy = element.Nodes[0].Y;
            double Dz = element.Nodes[0].Z+h/2;
            double Bx1 = Bx - Ax;
            double Bx2 = By - Ay;
            double Bx3 = Bz - Az;
            double Cx1 = Cx - Ax;
            double Cx2 = Cy - Ay;
            double Cx3 = Cz - Az;
            double Dx1 = Dx - Ax;
            double Dx2 = Dy - Ay;
            double Dx3 = Dz - Az;
            double L1 = Math.Sqrt(Math.Pow(Bx1, 2) + Math.Pow(Bx2, 2) + Math.Pow(Bx3, 2));
            double L2= Math.Sqrt(Math.Pow(Cx1, 2) + Math.Pow(Cx2, 2) + Math.Pow(Cx3, 2));
            double L3 = Math.Sqrt(Math.Pow(Dx1, 2) + Math.Pow(Dx2, 2) + Math.Pow(Dx3, 2));
            length = L1;
            t04.Clear();
            t04[0, 0] = Bx1 / L1;
            t04[0, 1] = Bx2 / L1;
            t04[0, 2] = Bx3 / L1;
            t04[1, 0] = Cx1 / L2;
            t04[1, 1] = Cx2 / L2;
            t04[1, 2] = Cx3 / L2;
            t04[2, 0] = Dx1 / L3;
            t04[2, 1] = Dx2 / L3;
            t04[2, 2] = Dx3 / L3;
            for (int i = 0; i < 3; i++)
                for (int j = 0; j < 3; j++)
                {
                    t04[i + 3, j + 3] = t04[i, j];
                    t04[i + 6, j + 6] = t04[i, j];
                    t04[i + 9, j + 9] = t04[i, j];
                }
            an.Clear();
            an[0, 0] = -1;
            an[0, 6] = 1;
            an[1, 5] = -1;
            an[1, 11] = 1;
            an[2, 1] = -2 / length;
            an[2, 5] = -1;
            an[2, 7] = 2 / length;
            an[2, 11] = -1;
            an[3, 4] = 1;
            an[3, 10] = -1;
            an[4, 2] = 2 / length;
            an[4, 4] = 1;
            an[4, 8] = -2 / length;
            an[4, 10] = 1;
            an[5, 3] = -1;
            an[5, 9] = 1; 
            transformationMatricesInitialized = true;
        }

        #region IElementType Members

        public int ID
        {
            get { throw new NotImplementedException(); }
        }

        public ElementDimensions ElementDimensions
        {
            get { return ElementDimensions.ThreeD; }
        }

        public IList<IList<DOFType>> GetElementDOFTypes(Element element)
        {
            return dofs;
        }

        public IList<Node> GetNodesForMatrixAssembly(Element element)
        {
            return element.Nodes;
        }

        public IMatrix2D<double> StiffnessMatrix(Element element)
        {
            if (!transformationMatricesInitialized) // in the first iteration
            {
                CalculateTransformationMatrices(element);
                elementStiffnessMatrix =BuildInitialStiffnessMatrix();
                //return dofEnumerator.GetTransformedMatrix(t04.Transpose() * an.Transpose() * BuildInitialStiffnessMatrix() * an * t04);
                return elementStiffnessMatrix; //ASK SER
            }
            else // in the next iterations CalculateStresses() has been called by the Analyzer (the Subdomain actually) and currentStiffnessMatrix has been set.
            {
                // return dofEnumerator.GetTransformedMatrix(t04.Transpose() * an.Transpose() * currentStiffnessMatrix * an * t04);
                return elementStiffnessMatrix; //ASK SER
            }
        }

        //public SymmetricMatrix2D<double> InverseMatrix(SymmetricMatrix2D<double> input)
        //{


        //    int size = input.Rows;
        //    var inputtt = input.ToMatrix2D();
        //    var helllppp = new double[size, size];
        //    for (int i1=0;i1<size;i1++)
        //        for(int j1=0;j1<size;j1++)
        //        {
        //            helllppp[i1, j1] = inputtt[i1, j1];
        //        }
        //    var helpppp = new SkylineMatrix2D<double>(helllppp);
        //    //helpppp.Factorize()
        //    SymmetricMatrix2D<double> result = new SymmetricMatrix2D<double>(size);
        //    for (int i=0; i<size; i++)
        //    {
        //        Vector <double> v= new Vector<double>(size);
        //        v[i] = 1;
        //        double[] d = new double[size];
        //        helpppp.Solve(v, d); //NEED TO FACTORIZE IT. BUT DONT KNOW HOW....
        //        for (int j = 0; j < size; j++)
        //        {
        //            result[j, i] = d[j];
        //        }
        //    }
        //    return result;
        //}

        //public Matrix2D<double> InverseMatrix(Matrix2D<double> input)
        //{
        //    int size = input.Rows;
        //    Matrix2D<double> result = new Matrix2D<double>(size,size);
        //    for (int i = 0; i < size; i++)
        //    {
        //        Vector<double> v = new Vector<double>(size);
        //        v[i] = 1;
        //        double[] d = new double[size];
        //        input.Solve(v, d);
        //        for (int j = 0; j < size; j++)
        //        {
        //            result[j, i] = d[j];
        //        }
        //    }
        //    return result;
        //}

        public IMatrix2D<double> MassMatrix(Element element)
        {
            if (!transformationMatricesInitialized) CalculateTransformationMatrices(element);
            double massDivided = Density * sectionArea * length / 4;

            Matrix2D<double> snb = new Matrix2D<double>(6, 6);
            snb[0, 0] = massDivided;
            snb[1, 1] = massDivided;
            snb[3, 3] = massDivided;
            snb[4, 4] = massDivided;

            return  t04.Transpose() * an.Transpose() * snb * an * t04;
        }

        public IMatrix2D<double> DampingMatrix(Element element)
        {
            //need to check it
            return   (new Matrix2D<double>(12, 12)) ;
        }

        public Tuple<double[], double[]> CalculateStresses(Element element, double[] localDisplacements, double[] rhoEl)
        {
            ElementStateDetermination(element, localDisplacements, rhoEl);
            return null;
        }

        public double[] CalculateForcesForLogging(Element element, double[] localDisplacements)
        {
            return CalculateForces(element, localDisplacements, new double[localDisplacements.Length]);
        }

        public double[] CalculateForces(Element element, double[] localDisplacements, double[] localdDisplacements)
        {
            return dcurrentInternalForces;
        }

        public void ElementStateDetermination(Element element, double[] localDisplacements, double[] rhoEl)
        {
            if (!fibersInitialized) throw new InvalidOperationException("Fibers have not been initialized");
            if (!transformationMatricesInitialized) throw new InvalidOperationException("Transformation matrices have not been initialized");
            Vector<double> residual = new Vector<double>(6);
            Vector<double> rhoEn = an * t04 * (new Vector<double>(rhoEl));
            double[] dPn = new double[6];
            firstiteration = true;
            for (int i = 0; i < 6; i++)
            {
                residual[i] = rhoEn[i];
            }
            while (residual.Norm > Math.Pow(10, -16))
            {
                int pointid = -1;
                Matrix2D<double> stiffnessMatrix = new Matrix2D<double>(6, 6);
                Vector<double> help1 = new Vector<double>(6);
                if (firstiteration)
                {
                    for (int i = 0; i < 6; i++)
                    {
                        help1[i] = rhoEn[i];
                    }
                    currentStiffnessMatrix.Multiply(help1, dPn);
                }
                else
                {
                    for (int i = 0; i < 6; i++)
                    {
                        help1[i] = - residual[i];
                    }
                    currentStiffnessMatrix.Multiply(help1, dPn);
                }
                residual = new Vector<double>(6);
                Vector<double> ddPn = new Vector<double>(6);
                for (int i = 0; i < 6; i++)
                {
                    ddPn[i] = dPn[i];
                }
                Vector<double> help2 = an.Transpose() * ddPn;
                if (firstiteration)
                {
                    currentInternalForces.Clear();
                    currentInternalForces.Add(currentInternalForcesBalanced);
                    currentInternalForces.Add(help2);
                    currentInternalForces.CopyTo(dcurrentInternalForces,0);
                }
                else
                {
                    currentInternalForces.Add(help2);
                    currentInternalForces.CopyTo(dcurrentInternalForces, 0);
                }
                foreach (GaussLobattoPoint1D point in points)
                {
                    pointid++;
                    b04[0, 0] = 1;
                    b04[1, 1] = -1;
                    b04[1, 2] = point.Coordinate;
                    b04[2, 3] = -1;
                    b04[2, 4] = point.Coordinate;
                    double[] dPsec = new double[3];
                    b04.Multiply(ddPn, dPsec);
                    double[] dispsec = new double[3];
                    Vector<double> ddPsec = new Vector<double>(3);
                    for (int i = 0; i < 3; i++)
                    {
                        ddPsec[i] = dPsec[i];
                    }
                    if (firstiteration)
                    {
                        point.sectionforces.Clear();
                        point.sectionforces.Add(point.sectionforcesbalanced);
                        point.sectionforces.Add(ddPsec);
                    }
                    else
                    {
                        point.sectionforces.Add(ddPsec);
                    }
                    Vector<double> df1 = new Vector<double>(3); 
                    for (int i = 0; i < 3; i++)
                    {
                        df1[i] = dPsec[i];
                    }
                    var help3 = new double[3];
                    point.sectionflexibilityMatrix.Multiply(df1, help3);
                    Vector<double> h3 = new Vector<double>(3);
                    for (int i = 0; i < 3; i++)
                    {
                        h3[i] = help3[i];
                    }
                    if (firstiteration)
                    {
                        point.deformations.Clear();
                        point.deformations.Add(point.deformationsbalanced);
                        point.deformations.Add(h3);
                    }
                    else
                    {
                        point.deformations.Add(h3);
                    }
                    var EDA = new double[fibers.Count()/4, 3];
                    var EYDA = new double[fibers.Count()/4, 3];
                    var EZDA = new double[fibers.Count()/4, 3];
                    var EY2DA = new double[fibers.Count()/4, 3];
                    var EZ2DA = new double[fibers.Count()/4, 3];
                    var EYZDA = new double[fibers.Count()/4, 3];
                    var GR2DA = new double[fibers.Count()/4, 3];
                    var SDA = new double[fibers.Count()/4, 3];
                    var SYDA = new double[fibers.Count()/4, 3];
                    var SZDA = new double[fibers.Count()/4, 3];
                    var hc = 0;
                    double dx = 0;
                    double dy = 0;
                    var size = fibers.Count() / 4;
                    Fiber[] Fibcheck = new Fiber[fibers.Count() / 4];
                    sections[0].CopyTo(pointid * fibers.Count() / 4, Fibcheck, 0, fibers.Count() / 4); //watch out. 
                    foreach (Fiber fiber in Fibcheck)
                    {
                        double x = fiber.X;
                        double y = fiber.Y;
                        double bb = fiber.B;
                        double hh = fiber.H;
                        dx = bb;
                        dy = hh;
                        double dStrain = (rhoEn[0] + y * rhoEn[1] + x * rhoEn[3] - 3 * y * point.Coordinate * rhoEn[2] -
                            3 * x * point.Coordinate * rhoEn[4]) / length;
                        fiber.Material.UpdateMaterial(dStrain); 
                        double e = fiber.Material.YoungModulus;
                        double n = fiber.Material.PoissonRatio;
                        double g = e / (2 * (1 + n));
                        double stresss = fiber.Material.Stress;
                        EDA[hc, 0] = fiber.X;
                        EDA[hc, 1] = fiber.Y;
                        EDA[hc, 2] = e;
                        EYDA[hc, 0] = fiber.X;
                        EYDA[hc, 1] = fiber.Y;
                        EYDA[hc, 2] = e*fiber.Y;
                        EZDA[hc, 0] = fiber.X;
                        EZDA[hc, 1] = fiber.Y;
                        EZDA[hc, 2] = e*fiber.X;
                        EY2DA[hc, 0] = fiber.X;
                        EY2DA[hc, 1] = fiber.Y;
                        EY2DA[hc, 2] = e*Math.Pow(fiber.Y,2);
                        EZ2DA[hc, 0] = fiber.X;
                        EZ2DA[hc, 1] = fiber.Y;
                        EZ2DA[hc, 2] = e * Math.Pow(fiber.X, 2);
                        EYZDA[hc, 0] = fiber.X;
                        EYZDA[hc, 1] = fiber.Y;
                        EYZDA[hc, 2] = e*fiber.X*fiber.Y;
                        GR2DA[hc, 0] = fiber.X;
                        GR2DA[hc, 1] = fiber.Y;
                        GR2DA[hc, 2] = g*(Math.Pow(fiber.X, 2)+ Math.Pow(fiber.Y, 2));
                        SDA[hc, 0] = fiber.X;
                        SDA[hc, 1] = fiber.Y;
                        SDA[hc, 2] = stresss;
                        SYDA[hc, 0] = fiber.X;
                        SYDA[hc, 1] = fiber.Y;
                        SYDA[hc, 2] = stresss*fiber.Y;
                        SZDA[hc, 0] = fiber.X;
                        SZDA[hc, 1] = fiber.Y;
                        SZDA[hc, 2] = stresss*fiber.X;
                        hc++;
                    }
                    point.sectionstiffnessMatrix[0, 0] = SimpsonIntegral(EDA,dx,dy,EDA[0,0],EDA[0,1],EDA[fibers.Count()/4-1,0],EDA[fibers.Count()/4-1,1]);
                    point.sectionstiffnessMatrix[0, 1] = -SimpsonIntegral(EYDA, dx, dy, EYDA[0, 0], EDA[0, 1], EYDA[fibers.Count()/4-1, 0], EYDA[fibers.Count()/4-1, 1]);
                    point.sectionstiffnessMatrix[0, 2] = SimpsonIntegral(EZDA, dx, dy, EZDA[0, 0], EZDA[0, 1], EZDA[fibers.Count()/4-1, 0], EZDA[fibers.Count()/4-1, 1]); 
                    point.sectionstiffnessMatrix[1, 0] = point.sectionstiffnessMatrix[0, 1];
                    point.sectionstiffnessMatrix[1, 1] = SimpsonIntegral(EY2DA, dx, dy, EY2DA[0, 0], EY2DA[0, 1], EY2DA[fibers.Count()/4-1, 0], EY2DA[fibers.Count()/4-1, 1]); 
                    point.sectionstiffnessMatrix[1, 2] = -SimpsonIntegral(EYZDA, dx, dy, EYZDA[0, 0], EYZDA[0, 1], EYZDA[fibers.Count()/4-1, 0], EYZDA[fibers.Count()/4-1, 1]);
                    point.sectionstiffnessMatrix[2, 0] = point.sectionstiffnessMatrix[0, 2];
                    point.sectionstiffnessMatrix[2, 1] = point.sectionstiffnessMatrix[1, 2];
                    point.sectionstiffnessMatrix[2, 2] = SimpsonIntegral(EZ2DA, dx, dy, EZ2DA[0, 0], EZ2DA[0, 1], EZ2DA[fibers.Count()/4-1, 0], EZ2DA[fibers.Count()/4-1, 1]);
                    //point.sectionflexibilityMatrix = InverseMatrix(point.sectionstiffnessMatrix);
                    point.sectionflexibilityMatrix = point.sectionstiffnessMatrix.Invert();
                    var sectionresisting = new Vector<double>(3);
                    sectionresisting[0]= SimpsonIntegral(SDA,dx,dy, SDA[0, 0], SDA[0, 1], SDA[fibers.Count()/4-1, 0], SDA[fibers.Count()/4-1, 1]);
                    sectionresisting[1] = SimpsonIntegral(SYDA,dx,dy, SYDA[0, 0], SYDA[0, 1], SYDA[fibers.Count()/4-1, 0], SYDA[fibers.Count()/4-1, 1]); 
                    sectionresisting[2] = SimpsonIntegral(SZDA,dx,dy, SZDA[0, 0], SZDA[0, 1], SZDA[fibers.Count()/4-1, 0], SZDA[fibers.Count()/4-1, 1]);
                    TorsionS= SimpsonIntegral(GR2DA,dx,dy, GR2DA[0, 0], GR2DA[0, 1], GR2DA[fibers.Count()/4-1, 0], GR2DA[fibers.Count()/4-1, 1]);
                    sectionresisting.Add(point.sectionforces);
                    point.sectionflexibilityMatrix.Multiply(sectionresisting, point.residualdef);
                    b04[0, 0] = 1;
                    b04[1, 1] = -1;
                    b04[1, 2] = point.Coordinate;
                    b04[2, 3] = -1;
                    b04[2, 4] = point.Coordinate;
                    Matrix2D<double> fl = point.sectionflexibilityMatrix.ToMatrix2D();
                    var h = b04.Transpose() * fl * b04;
                    h.Scale(length*point.WeightFactor/2);
                    stiffnessMatrix = Matrix2D<double>.Add(h,stiffnessMatrix);
                    Matrix2D<double> rf = new Matrix2D<double>(3,1);
                    for (int i = 0; i < 3; i++)
                    {
                        rf[i,0] = point.residualdef[i];
                    }
                    Matrix2D<double> h1 = new Matrix2D<double>(6, 1) ;
                    h1 = b04.Transpose()*rf;
                    h1.Scale(length * point.WeightFactor / 2);
                    for (int i = 0; i < 6; i++)
                    {
                        residual[i] += h1[i,0];
                    }
                }
                this.currentFlexibilityMatrix = stiffnessMatrix; //in natural modes terms.
                Matrix2D<double> matrixinv1 = new Matrix2D<double>(5, 5);
                Matrix2D<double> matrixinv2 = new Matrix2D<double>(5, 5);
                for (int ii = 0; ii < 5; ii++)
                {
                    for (int jj = 0; jj < 5; jj++)
                    {
                        matrixinv1[ii, jj] = this.currentFlexibilityMatrix[ii, jj];
                    }
                }
                //this.currentStiffnessMatrix = InverseMatrix(stiffnessMatrix);
                matrixinv2=matrixinv1.Invert();
                for (int ii = 0; ii < 5; ii++)
                {
                    for (int jj = 0; jj < 5; jj++)
                    {
                        this.currentStiffnessMatrix[ii, jj] = matrixinv2[ii, jj];
                    }
                }
                this.currentStiffnessMatrix[5, 5] = TorsionS / length;
                firstiteration = false;
                var help12 = an * t04;
                this.elementStiffnessMatrix =help12.Transpose() * currentStiffnessMatrix * help12;
            }
        }
        public static double SimpsonIntegral(double[,] fxy, double dx, double dy, double x0, double y0, double xf, double yf)
        {
            double SimInt = 0;
            var h1 = fxy.GetLength(0);
            var h2 = fxy.GetLength(1);
            double af = (xf - x0) / dx;
            int aaf = Convert.ToInt32(af);
            double bf = (yf - y0) / dy;
            int bbf = Convert.ToInt32(bf);
            for (int i = 0; i < h1; i++)
            {
                bool hasbeenadded = false;
                double a = (fxy[i, 0] - x0) / dx;
                int aa = Convert.ToInt32(a);
                double b = (fxy[i, 1] - y0) / dy;
                int bb = Convert.ToInt32(b);
                if ((aa == 0 || aa == aaf) && (bb == 0 || bb == bbf) && !hasbeenadded)
                {
                    SimInt += fxy[i, 2];
                    hasbeenadded = true;
                }
                if (aa == 0 && bb % 2 == 0 && !hasbeenadded)
                {
                    SimInt += fxy[i, 2] * 2;
                    hasbeenadded = true;
                }
                if (aa == 0 && bb % 2 == 1 && !hasbeenadded)
                {
                    SimInt += fxy[i, 2] * 4;
                    hasbeenadded = true;
                }
                if (aa == aaf && bb % 2 == 0 && !hasbeenadded)
                {
                    SimInt += fxy[i, 2] * 2;
                    hasbeenadded = true;
                }
                if (aa == aaf && bb % 2 == 1 && !hasbeenadded)
                {
                    SimInt += fxy[i, 2] * 4;
                    hasbeenadded = true;
                }
                if (bb == 0 && aa % 2 == 0 && !hasbeenadded)
                {
                    SimInt += fxy[i, 2] * 2;
                    hasbeenadded = true;
                }
                if (bb == 0 && aa % 2 == 1 && !hasbeenadded)
                {
                    SimInt += fxy[i, 2] * 4;
                    hasbeenadded = true;
                }
                if (bb == bbf && aa % 2 == 0 && !hasbeenadded)
                {
                    SimInt += fxy[i, 2] * 2;
                    hasbeenadded = true;
                }
                if (bb == bbf && aa % 2 == 1 && !hasbeenadded)
                {
                    SimInt += fxy[i, 2] * 4;
                    hasbeenadded = true;
                }
                if (aa % 2 == 0 && bb % 2 == 0 && !hasbeenadded)
                {
                    SimInt += fxy[i, 2] * 4;
                    hasbeenadded = true;
                }
                if (aa % 2 == 1 && bb % 2 == 0 && !hasbeenadded)
                {
                    SimInt += fxy[i, 2] * 8;
                    hasbeenadded = true;
                }
                if (aa % 2 == 0 && bb % 2 == 1 && !hasbeenadded)
                {
                    SimInt += fxy[i, 2] * 8;
                    hasbeenadded = true;
                }
                if (aa % 2 == 1 && bb % 2 == 1 && !hasbeenadded)
                {
                    SimInt += fxy[i, 2] * 16;
                    hasbeenadded = true;
                }
            }
            SimInt = SimInt * dx * dy / 9;
            return SimInt;
        }
        public double TrapezoidIntegral(double[,] fxy, double dx, double dy, double x0, double y0, double xf, double yf)
        {
            double TrapInt = 0;
            var h1 = fxy.GetLength(0);
            var h2 = fxy.GetLength(1);
            double af = (xf - x0) / dx;
            int aaf = Convert.ToInt32(af);
            double bf = (yf - y0) / dy;
            int bbf = Convert.ToInt32(bf);
            for (int i = 0; i < h1; i++)
            {
                bool hasbeenadded = false;
                double a = (fxy[i, 0] - x0) / dx;
                int aa = Convert.ToInt32(a);
                double b = (fxy[i, 1] - y0) / dy;
                int bb = Convert.ToInt32(b);
                if ((aa == 0 || aa == aaf) && (bb == 0 || bb == bbf) && !hasbeenadded)
                {
                    TrapInt += fxy[i, 2];
                    hasbeenadded = true;
                }
                if (aa == 0 && !hasbeenadded)
                {
                    TrapInt += fxy[i, 2] * 2;
                    hasbeenadded = true;
                }
                if (aa == aaf && !hasbeenadded)
                {
                    TrapInt += fxy[i, 2] * 2;
                    hasbeenadded = true;
                }
                if (bb == 0 && !hasbeenadded)
                {
                    TrapInt += fxy[i, 2] * 2;
                    hasbeenadded = true;
                }
                if (bb == bbf && !hasbeenadded)
                {
                    TrapInt += fxy[i, 2] * 2;
                    hasbeenadded = true;
                }
                if (aa!=0 && bb!=0 && aa!=aaf && bb!=bbf && !hasbeenadded)
                {
                    TrapInt += fxy[i, 2] * 4;
                    hasbeenadded = true;
                }
            }
            TrapInt = TrapInt * dx * dy / 4;
            return TrapInt;
        }
        public double[] CalculateAccelerationForces(Element element, IList<MassAccelerationLoad> loads)
        {
            Vector<double> accelerations = new Vector<double>(noOfDOFs);
            IMatrix2D<double> massMatrix = MassMatrix(element);

            foreach (MassAccelerationLoad load in loads)
            {
                int index = 0;
                foreach (DOFType[] nodalDOFTypes in GetElementDOFTypes(element))
                    foreach (DOFType dofType in nodalDOFTypes)
                    {
                        if (dofType == load.DOF) accelerations[index] += load.Amount;
                        index++;
                    }
            }
            double[] forces = new double[noOfDOFs];
            massMatrix.Multiply(accelerations, forces);
            return forces;
        }

        public bool MaterialModified
        {
          get { return material.Modified; }
        }

        public void ResetMaterialModified()
        {
            material.ResetModified();
        }

        public void ClearMaterialState()
        {
            //foreach (IFiberMaterial m in material.FiberMaterials) m.ClearState();
        }

        public void SaveMaterialState()
        {
            currentInternalForces.CopyTo(dcurrentInternalForcesBalanced, 0);
            currentInternalForcesBalanced.Clear();
            currentInternalForcesBalanced.Add(currentInternalForces);
            foreach (GaussLobattoPoint1D point in this.points)
            {
                point.deformationsbalanced.Clear();
                point.deformationsbalanced.Add(point.deformations);
                point.sectionforcesbalanced.Clear();
                point.sectionforcesbalanced.Add(point.sectionforces);
            }
            foreach (Fiber fiber in sections[0])
                {
                    fiber.Material.SaveState();
                }
            
        }
        #endregion

        #region IFiberFiniteElement Members

        public IList<IFiber> Fibers
        {
            get { return (IList<IFiber>)fibers; }
        }

        public IFiberFiniteElementMaterial Material
        {
            get { return material; }
        }

        #endregion

        //#region IFiniteElement Members


        //IFiniteElementMaterial IFiniteElement.Material
        //{
        //    get { return material; }
        //}

        //#endregion

        #region IFiniteElement Members


        public void ClearMaterialStresses()
        {
            foreach (Fiber fiber in fibers) fiber.Material.ClearStresses();
        }

        #endregion

        private Matrix2D<double> BuildInitialStiffnessMatrix()
        {
            var help=new Matrix2D<double>(6, 6);
            help[0, 0] = Material.YoungModulus * SectionArea / length;
            help[1, 1] = Material.YoungModulus * Iyy / length;
            help[2, 2] = 3 * help[1, 1];
            help[3, 3] = Material.YoungModulus * Izz / length;
            help[4, 4] = 3 * help[3, 3];
            help[5, 5] = Material.YoungModulus * It / ((2*(1+Material.PoissonRatio))*length);
            this.currentStiffnessMatrix = help;
            elementStiffnessMatrix= t04.Transpose() * an.Transpose() * help * an * t04;
            return elementStiffnessMatrix;
        }
    }
}
