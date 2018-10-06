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
        //modified by Ambrosios Savvides. The calculate rot transformation stuff is not checked thus the connectivity between rotations and nearby
        // displacements needs to be checked and added.
        private static readonly DOFType[] noDOFTypes = new DOFType[0];
        private static readonly DOFType[] nodalDOFTypes = new DOFType[6] { DOFType.X, DOFType.Y, DOFType.Z, DOFType.RotX, DOFType.RotY, DOFType.RotZ };
        //private static readonly DOFType[][] dofs = new DOFType[][] { nodalDOFTypes, nodalDOFTypes, noDOFTypes };
        private static readonly DOFType[][] dofs = new DOFType[][] { nodalDOFTypes, nodalDOFTypes };
        GaussLobattoPoint1D[] points = new GaussLobattoPoint1D[4];
        private int noOfDOFs = 12;
        private  Matrix2D<double> tempstiffnessMatrix = new Matrix2D<double>(6, 6);
        private Matrix2D<double> tempflexibilityMatrix = new Matrix2D<double>(6, 6);
        private bool transformationMatricesInitialized = false;
        private bool fibersInitialized = false;
        private readonly IFiberFiniteElementMaterial material;
        private readonly List<List<Fiber>> sections;
        private readonly List<Fiber> fibers;
        private readonly double b, h;
        private double length, sectionArea, iyy, izz, iyz;
        private Matrix2D<double> rotTransformation;
        private readonly Matrix2D<double> t04 = new Matrix2D<double>(12, 12);
        private readonly Matrix2D<double> an = new Matrix2D<double>(6, 12);
        private readonly Matrix2D<double> b04 = new Matrix2D<double>(6, 6);
        private readonly Node[][] rotNodes;
        private IFiniteElementDOFEnumerator dofEnumerator = new GenericDOFEnumerator();

        public FiberBeam3D(IFiberFiniteElementMaterial material, int noOfFibers, double b, double h)
        {
            this.b = b;
            this.h = h;
            this.material = material;
            points= GaussLobatto.GetGaussLobattoPoints();
            sections = new List<List<Fiber>>(4);
            fibers = new List<Fiber>(noOfFibers); // The row of the sections is the row given in Gauss Lobatto class.
            InitializeFibers();
        }

        public FiberBeam3D(IFiberFiniteElementMaterial material, IFiniteElementDOFEnumerator dofEnumerator, int noOfFibers, double b, double h,
            Node[] rot1Nodes, Node[] rot2Nodes) : this(material, noOfFibers, b, h)
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

        private void InitializeFibers()
        {
            int noOfFibers = fibers.Capacity;
		    int rows = Math.Max(1, (int)(Math.Sqrt(noOfFibers * h / b)));
		    int cols = Math.Max(1, (int)(noOfFibers / rows));
		    noOfFibers = rows * cols;
            double fiberB = b / cols;
            double fiberH = h / rows;

            int pos = 0;
            for (int i = 0; i < rows; i++)
                for (int j = 0; j < cols; j++)
                {
			        double x = -b * 0.5 + j * fiberB;
			        double y = -h * 0.5 + i * fiberH;
                    Fiber fiber = new Fiber(fiberB, fiberH, x, y) { Material = material.FiberMaterials[pos] };
                    fibers.Add(fiber);
                    pos++;
                }

            sectionArea = 0;
            iyy = 0;
            izz = 0;
            foreach (Fiber fiber in fibers)
            {
                sectionArea += fiber.B * fiber.H;
        		iyy += fiber.B * fiber.H * Math.Pow(fiber.Y - fiber.H / 2, 2);
        		izz += fiber.B * fiber.H * Math.Pow(fiber.X - fiber.B / 2, 2);
            }
            iyz = Math.Pow(b, 3) * h / 3;
            sections[0] = fibers;
            sections[1] = fibers;
            sections[2] = fibers;
            sections[3] = fibers;
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
            double Cx = element.Nodes[0].X-h/2;
            double Cy = element.Nodes[0].Y;
            double Cz = element.Nodes[0].Z;
            double Dx = element.Nodes[0].X;
            double Dy = element.Nodes[0].Y-b/2;
            double Dz = element.Nodes[0].Z;
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
            if (!fibersInitialized) throw new InvalidOperationException("Fibers have not been initialized");
            //if (!transformationMatricesInitialized) throw new InvalidOperationException("Transformation matrices have not been initialized");
            if (!transformationMatricesInitialized) CalculateTransformationMatrices(element);
            Matrix2D<double> stiffnessMatrix = new Matrix2D<double>(12,12);
            Matrix2D<double> flexibilityMatrix = new Matrix2D<double>(6,6);
            int pointid = -1;
            foreach (GaussLobattoPoint1D point in points)
            {
                pointid++;
                double k11 = 0;
                double k12 = 0;
                double k13 = 0;
                double k22 = 0;
                double k23 = 0;
                double k33 = 0;
                double kg = 0;
                double kt = 0;
                foreach (Fiber fiber in sections[pointid])
                {
                    double x = fiber.X;
                    double y = fiber.Y;
                    double bb = fiber.B;
                    double hh = fiber.H;
                    double ym = y - hh / 2;
                    double xm = x - bb / 2;
                    double e = fiber.Material.YoungModulus;
                    double n = fiber.Material.PoissonRatio;
                    double g = e / (2 * (1 + n));
                    double A = bb * hh;
                    k11 += e * A;
                    k12 += -ym * e * A;
                    k13 += xm * e * A;
                    k22 += Math.Pow(ym, 2) * e * A;
                    k23 += -ym * xm * e * A;
                    k33+= Math.Pow(xm, 2) * e * A;
                    kg += g*A;
                    kt += g * A * Math.Sqrt(Math.Pow(xm, 2) + Math.Pow(ym, 2));                    
                }
                point.sectionstiffnessMatrix[0, 0] = k11;
                point.sectionstiffnessMatrix[0, 1] = k12;
                point.sectionstiffnessMatrix[0, 2] = k13;
                point.sectionstiffnessMatrix[1, 0] = point.sectionstiffnessMatrix[0, 1];
                point.sectionstiffnessMatrix[1, 1] = k22;
                point.sectionstiffnessMatrix[1, 2] = k23;
                point.sectionstiffnessMatrix[2, 0] = point.sectionstiffnessMatrix[0, 2];
                point.sectionstiffnessMatrix[2, 1] = point.sectionstiffnessMatrix[1, 2];
                point.sectionstiffnessMatrix[2, 2] = k33;
                point.sectionstiffnessMatrix[3, 3] = kg;
                point.sectionstiffnessMatrix[4, 4] = kg;
                point.sectionstiffnessMatrix[5, 5] = kt;
                point.sectionflexibilityMatrix = InverseMatrix(point.sectionstiffnessMatrix);
                b04[0, 0] = 1;
                b04[1, 1] = -1;
                b04[1, 2] = point.Coordinate;
                b04[2, 3] = -1;
                b04[2, 4] = point.Coordinate;
                b04[3, 2] = 2 / length;
                b04[4, 4] = 2 / length;
                b04[5, 5] = 1;
                Matrix2D<double> fl = point.sectionflexibilityMatrix.ToMatrix2D(); 
                var h = b04.Transpose() * fl * b04;
                h.Scale(length / 2);
                Matrix2D<double>.Add(h,stiffnessMatrix);
            }
            this.tempstiffnessMatrix = InverseMatrix(stiffnessMatrix);
            return dofEnumerator.GetTransformedMatrix(t04.Transpose() * an.Transpose() * tempstiffnessMatrix * an * t04);
        }

        public SymmetricMatrix2D<double> InverseMatrix(SymmetricMatrix2D<double> input)
        {
            int size = input.Rows;
            SymmetricMatrix2D<double> result = new SymmetricMatrix2D<double>(size);
            for (int i=0; i<size; i++)
            {
                Vector <double> v= new Vector<double>(size);
                v[size] = 1;
                double[] d = new double[size];
                input.Solve(v, d);
                for (int j = 0; j < size; j++)
                {
                    result[j, i] = d[j];
                }
            }
            return result;
        }

        public Matrix2D<double> InverseMatrix(Matrix2D<double> input)
        {
            int size = input.Rows;
            Matrix2D<double> result = new Matrix2D<double>(size,size);
            for (int i = 0; i < size; i++)
            {
                Vector<double> v = new Vector<double>(size);
                v[size] = 1;
                double[] d = new double[size];
                input.Solve(v, d);
                for (int j = 0; j < size; j++)
                {
                    result[j, i] = d[j];
                }
            }
            return result;
        }

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
            if (!fibersInitialized) throw new InvalidOperationException("Fibers have not been initialized");
            if (!transformationMatricesInitialized) throw new InvalidOperationException("Transformation matrices have not been initialized");
            Vector<double> rhoEn = an * t04 * (new Vector<double>(rhoEl));
            double[] dPn = new double[6];
            tempstiffnessMatrix.Multiply(rhoEn, dPn);
            Vector<double> ddPn = new Vector<double>(6);
            for (int i = 0; i < 6; i++)
            {
                ddPn[i] = dPn[i];
            }
            int pointid = -1;
            foreach (GaussLobattoPoint1D point in points)
            {
                pointid++;
                b04[0, 0] = 1;
                b04[1, 1] = -1;
                b04[1, 2] = point.Coordinate;
                b04[2, 3] = -1;
                b04[2, 4] = point.Coordinate;
                b04[3, 2] = 2 / length;
                b04[4, 4] = 2 / length;
                b04[5, 5] = 1;
                double[] dPsec = new double[6];
                b04.Multiply(ddPn, dPsec);
                double[] dispsec = new double[6];
                Vector<double> ddPsec = new Vector<double>(6);
                for (int i = 0; i < 6; i++)
                {
                    ddPsec[i] = dPsec[i];
                }
                point.sectionflexibilityMatrix.Multiply(ddPsec, dispsec);
                foreach (Fiber fiber in sections[pointid])
                {
                    double x = fiber.X;
                    double y = fiber.Y;
                    double bb = fiber.B;
                    double hh = fiber.H;
                    double ym = y - hh / 2;
                    double xm = x - bb / 2;
                    double dStrain = (dispsec[0] + ym * dispsec[1] + xm * dispsec[3] - 3 * ym * point.Coordinate * dispsec[2] -
                        3 * xm * point.Coordinate * dispsec[4]) / length;
                    fiber.Material.UpdateMaterial(dStrain);
                }
            }
            return new Tuple<double[], double[]>(new double[6], new double[6]); //check with serafeim what to do with it.
        }

        public double[] CalculateForcesForLogging(Element element, double[] localDisplacements)
        {
            return CalculateForces(element, localDisplacements, new double[localDisplacements.Length]);
        }

        public double[] CalculateForces(Element element, double[] localDisplacements, double[] localdDisplacements)
        {
            double[] f = new double[6];
            int pointid = -1;
            foreach (GaussLobattoPoint1D point in points)
            {
                pointid++;
                point.sectionforces = new double[6];                
              foreach (Fiber fiber in sections[pointid])
                {
                    double x = fiber.X;
                    double y = fiber.Y;
                    double bb = fiber.B;
                    double hh = fiber.H;
                    double ym = y - hh / 2;
                    double xm = x - bb / 2;
                    point.sectionforces[0] += fiber.Material.Stress * bb * hh;
                    point.sectionforces[1] += fiber.Material.Stress * bb * hh *(-ym);
                    point.sectionforces[2] += fiber.Material.Stress * bb * hh * (xm);
                }
              //this needs to be reset with serafeim.
            }
            return new double[6];
        }

        public void ElementStateDetermination(Element element, double[] localDisplacements, double[] rhoEl)
        {
            if (!fibersInitialized) throw new InvalidOperationException("Fibers have not been initialized");
            if (!transformationMatricesInitialized) throw new InvalidOperationException("Transformation matrices have not been initialized");
            Vector<double> residual = new Vector<double>(6);
            for (int i = 0; i < 6; i++)
            {
                residual[i] = rhoEl[i];
            }
            Vector<double> rhoEn = an * t04 * (new Vector<double>(rhoEl));
            Vector<double> rhoEnnew = rhoEn;
            double[] dPn = new double[6];
            tempstiffnessMatrix.Multiply(rhoEn, dPn);
            Vector<double> ddPn = new Vector<double>(6);
            for (int i = 0; i < 6; i++)
            {
                ddPn[i] = dPn[i];
            }
            while (residual.Norm > Math.Pow(10, -3))
            {
                int pointid = -1;
                residual = new Vector<double>(6);
                Matrix2D<double> stiffnessMatrix = new Matrix2D<double>(12, 12);
                rhoEn = rhoEnnew;
                tempstiffnessMatrix.Multiply(rhoEn, dPn);
                for (int i = 0; i < 6; i++)
                {
                    ddPn[i] = dPn[i];
                }
                foreach (GaussLobattoPoint1D point in points)
                {
                    pointid++;
                    double k11 = 0;
                    double k12 = 0;
                    double k13 = 0;
                    double k22 = 0;
                    double k23 = 0;
                    double k33 = 0;
                    double kg = 0;
                    double kt = 0;
                    b04[0, 0] = 1;
                    b04[1, 1] = -1;
                    b04[1, 2] = point.Coordinate;
                    b04[2, 3] = -1;
                    b04[2, 4] = point.Coordinate;
                    b04[3, 2] = 2 / length;
                    b04[4, 4] = 2 / length;
                    b04[5, 5] = 1;
                    double[] dPsec = new double[6];
                    b04.Multiply(ddPn, dPsec);
                    double[] dispsec = new double[6];
                    Vector<double> ddPsec = new Vector<double>(6);
                    for (int i = 0; i < 6; i++)
                    {
                        ddPsec[i] = dPsec[i];
                    }
                    point.sectionflexibilityMatrix.Multiply(ddPsec, dispsec);
                    point.sectionforcesprev = point.sectionforces;
                    point.sectionforces = new double[6];
                    foreach (Fiber fiber in sections[pointid])
                    {
                        double x = fiber.X;
                        double y = fiber.Y;
                        double bb = fiber.B;
                        double hh = fiber.H;
                        double ym = y - hh / 2;
                        double xm = x - bb / 2;
                        double dStrain = (dispsec[0] + ym * dispsec[1] + xm * dispsec[3] - 3 * ym * point.Coordinate * dispsec[2] -
                            3 * xm * point.Coordinate * dispsec[4]) / length;
                        fiber.Material.UpdateMaterial(dStrain);
                        double e = fiber.Material.YoungModulus;
                        double n = fiber.Material.PoissonRatio;
                        double g = e / (2 * (1 + n));
                        double A = bb * hh;
                        k11 += e * A;
                        k12 += -ym * e * A;
                        k13 += xm * e * A;
                        k22 += Math.Pow(ym, 2) * e * A;
                        k23 += -ym * xm * e * A;
                        k33 += Math.Pow(xm, 2) * e * A;
                        kg += g * A;
                        kt += g * A * Math.Sqrt(Math.Pow(xm, 2) + Math.Pow(ym, 2));
                        point.sectionforces[0] += fiber.Material.Stress * bb * hh;
                        point.sectionforces[1] += fiber.Material.Stress * bb * hh * (-ym);
                        point.sectionforces[2] += fiber.Material.Stress * bb * hh * (xm);
                    }
                    point.sectionstiffnessMatrix[0, 0] = k11;
                    point.sectionstiffnessMatrix[0, 1] = k12;
                    point.sectionstiffnessMatrix[0, 2] = k13;
                    point.sectionstiffnessMatrix[1, 0] = point.sectionstiffnessMatrix[0, 1];
                    point.sectionstiffnessMatrix[1, 1] = k22;
                    point.sectionstiffnessMatrix[1, 2] = k23;
                    point.sectionstiffnessMatrix[2, 0] = point.sectionstiffnessMatrix[0, 2];
                    point.sectionstiffnessMatrix[2, 1] = point.sectionstiffnessMatrix[1, 2];
                    point.sectionstiffnessMatrix[2, 2] = k33;
                    point.sectionstiffnessMatrix[3, 3] = kg;
                    point.sectionstiffnessMatrix[4, 4] = kg;
                    point.sectionstiffnessMatrix[5, 5] = kt;
                    point.sectionflexibilityMatrix = InverseMatrix(point.sectionstiffnessMatrix);
                    Vector<double> df = new Vector<double>(6);
                    for (int i = 0; i < 6; i++)
                    {
                        df[i] = point.sectionforcesprev[i]-point.sectionforces[i];
                    }
                    point.sectionflexibilityMatrix.Multiply(df, point.residualdef);
                    Vector<double> df1 = new Vector<double>(6);
                    for (int i = 0; i < 6; i++)
                    {
                        df1[i] = point.sectionforces[i];
                    }
                    point.sectionflexibilityMatrix.Multiply(df1, point.deformations);
                    b04[0, 0] = 1;
                    b04[1, 1] = -1;
                    b04[1, 2] = point.Coordinate;
                    b04[2, 3] = -1;
                    b04[2, 4] = point.Coordinate;
                    b04[3, 2] = 2 / length;
                    b04[4, 4] = 2 / length;
                    b04[5, 5] = 1;
                    Matrix2D<double> fl = point.sectionflexibilityMatrix.ToMatrix2D();
                    var h = b04.Transpose() * fl * b04;
                    h.Scale(length / 2);
                    Matrix2D<double>.Add(h, stiffnessMatrix);
                    Vector<double> rf = new Vector<double>(6);
                    for (int i = 0; i < 6; i++)
                    {
                        rf[i] = point.residualdef[i];
                    }
                    double[] h1 = new double[6];
                    b04.Transpose().Multiply(rf,h1);
                    for (int i = 0; i < 6; i++)
                    {
                        residual[i]+=h1[i];
                    }
                    Vector<double> ddf = new Vector<double>(6);
                    for (int i = 0; i < 6; i++)
                    {
                        ddf[i] = point.deformations[i];
                    }
                    h1 = new double[6];
                    b04.Transpose().Multiply(ddf, h1);
                    for (int i = 0; i < 6; i++)
                    {
                        rhoEnnew[i] += h1[i];
                    }
                }
                residual.Scale(length / 2);
                rhoEnnew.Scale(length / 2);
                this.tempflexibilityMatrix = stiffnessMatrix; //in natural modes terms.
                this.tempstiffnessMatrix = InverseMatrix(stiffnessMatrix);
                dofEnumerator.GetTransformedMatrix(t04.Transpose() * an.Transpose() * tempstiffnessMatrix * an * t04);
            }
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
            foreach (IFiberMaterial m in material.FiberMaterials) m.SaveState();
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
    }
}
