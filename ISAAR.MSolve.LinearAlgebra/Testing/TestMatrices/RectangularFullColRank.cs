﻿using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ISAAR.MSolve.LinearAlgebra.Factorizations;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.LinearAlgebra.Output;
using ISAAR.MSolve.LinearAlgebra.Testing.Utilities;
using ISAAR.MSolve.LinearAlgebra.Vectors;

namespace ISAAR.MSolve.LinearAlgebra.Testing.TestMatrices
{
    class RectangularFullColRank
    {
        public const int numRows = 10;
        public const int numCols = 5;

        public static readonly double[,] matrix = new double[,] {
            { 1.687338, 1.910625, 0.078085, 2.332444, 3.131424 },
            { 3.215230, 1.767679, 1.883659, 0.809185, 2.026953 },
            { 0.981285, 1.325970, 2.639193, 0.923645, 3.475809 },
            { 0.751857, 2.296387, 2.101308, 1.871156, 3.148936 },
            { 3.446088, 3.932793, 0.828995, 0.343996, 1.349877 },
            { 0.779645, 1.764859, 1.694237, 2.854230, 1.260338 },
            { 0.968975, 3.746153, 2.681011, 2.329230, 2.481479 },
            { 3.581092, 1.647008, 3.430712, 2.409152, 3.780291 },
            { 1.536053, 0.897209, 3.872541, 1.485733, 0.601071 },
            { 3.186167, 0.066731, 1.536828, 2.245892, 2.542407 }};

        public static readonly double[] lhs5 = { 3.4484, 1.9563, 2.7385, 4.2828, 5.3064 };
        public static readonly double[] rhs5 = MatrixOperations.MatrixTimesVector(matrix, lhs5);

        public static readonly double[] lhs10 = { 3.4484, 1.9563, 2.7385, 4.2828, 5.3064, 4.3251, 0.1117, 4.0487, 2.6311, 2.6269 };
        public static readonly double[] rhs10 = MatrixOperations.MatrixTimesVector(MatrixOperations.Transpose(matrix), lhs10);

        public static readonly double[,] qrFactorQ = {
            { -0.23075142099349710,  0.14858764284193901,  0.31792868485316900,  0.52033887561493763,  0.25119640190154863, -0.31609052075075511, -0.14778642143763676, -0.20713121727488037,  0.38894313590494700, -0.41690588846886784 },
            { -0.43969785029491520, -0.10772975041878161,  0.07417564455815109, -0.25365561595287567,  0.00742163793160727,  0.50045187775446243,  0.46428134845573432, -0.34723518223729494,  0.05832008220679226, -0.36766705946723366 },
            { -0.13419534684195095,  0.13124334022410436, -0.36103200437303257, -0.19632113866951575,  0.62338833887492229, -0.15238342744727723, -0.24643597223954369, -0.33041361660005214, -0.45550611799719420, -0.08159225990019033 },
            { -0.10281998694624772,  0.36775989343103105, -0.18811764198214728,  0.11534039321239833,  0.32455249539684250,  0.52520662918835392, -0.15271034272843115,  0.07115730533531130,  0.45778828012345169,  0.43169295990059880 },
            { -0.47126876942772489,  0.30986638208644840,  0.50635053078230396, -0.41979393109143270, -0.20998081083827991, -0.09764449684237536, -0.37794948782490234,  0.03273127635904442, -0.10261133630329747,  0.20297325141250672 },
            { -0.10662012686283072,  0.25272112351608017, -0.14061937078687231,  0.51532130006775501, -0.40692322124504376,  0.38027807781600009, -0.26111733434638934, -0.05636543553210413, -0.47734907759472461, -0.17324820800576807 },
            { -0.13251189634630042,  0.63811348815517166, -0.17161054696409667,  0.07128836704429709, -0.09596999203053819, -0.34460294522294443,  0.62009957715162045,  0.06108149541389936, -0.09053444009338751,  0.12769970133355366 },
            { -0.48973120246710766, -0.18714506899997330, -0.23356674775948780,  0.03581993650560919,  0.14709353964540772, -0.03588951204723363, -0.01271828595741990,  0.77029706318046287, -0.03789988799293117, -0.22676704045448226 },
            { -0.21006248450003742, -0.04046691405280455, -0.60620280595936349, -0.17752093128642521, -0.44985545375115621, -0.23659125205442419, -0.25089634999754107, -0.25248284816875666,  0.40429678229270199, -0.08211487942414421 },
            { -0.43572334812147162, -0.45837653401074380,  0.00607081461594355,  0.36523844289427687, -0.00167830105033348, -0.14286866951346483,  0.14102683949598455, -0.23917224941218027, -0.12956106417030591,  0.59355527534776520 }};

        public static readonly double[,] qrFactorR = {
            { -7.312362336644295, -5.194296836207316, -5.506302126480993, -4.455957642566284, -6.588771140414244 },
            {  0.000000000000000,  4.792011033788077,  1.820469593879060,  1.842539394813214,  2.284251661879478 },
            {  0.000000000000000,  0.000000000000000, -4.601665076558234, -1.960507008648889, -1.852780335086659 },
            {  0.000000000000000,  0.000000000000000,  0.000000000000000,  3.178213946676662,  2.013079671190600 },
            {  0.000000000000000,  0.000000000000000,  0.000000000000000,  0.000000000000000,  3.237357394936110 },
            {  0.000000000000000,  0.000000000000000,  0.000000000000000,  0.000000000000000,  0.000000000000000 },
            {  0.000000000000000,  0.000000000000000,  0.000000000000000,  0.000000000000000,  0.000000000000000 },
            {  0.000000000000000,  0.000000000000000,  0.000000000000000,  0.000000000000000,  0.000000000000000 },
            {  0.000000000000000,  0.000000000000000,  0.000000000000000,  0.000000000000000,  0.000000000000000 },
            {  0.000000000000000,  0.000000000000000,  0.000000000000000,  0.000000000000000,  0.000000000000000 }};


        /// <summary>
        /// Actually since I am taking LQ(transpose(A)), this is the transpose of <see cref="qrFactorR"/>
        /// </summary>
        public static readonly double[,] lqFactorL = {
            { -7.312362336644295, 0.000000000000000,  0.000000000000000, 0.000000000000000, 0.000000000000000, 0.000000000000000, 0.000000000000000, 0.000000000000000, 0.000000000000000, 0.000000000000000 },
            { -5.194296836207316, 4.792011033788077,  0.000000000000000, 0.000000000000000, 0.000000000000000, 0.000000000000000, 0.000000000000000, 0.000000000000000, 0.000000000000000, 0.000000000000000 },
            { -5.506302126480993, 1.820469593879060, -4.601665076558234, 0.000000000000000, 0.000000000000000, 0.000000000000000, 0.000000000000000, 0.000000000000000, 0.000000000000000, 0.000000000000000 },
            { -4.455957642566284, 1.842539394813214, -1.960507008648889, 3.178213946676662, 0.000000000000000, 0.000000000000000, 0.000000000000000, 0.000000000000000, 0.000000000000000, 0.000000000000000 },
            { -6.588771140414244, 2.284251661879478, -1.852780335086659, 2.013079671190600, 3.237357394936110, 0.000000000000000, 0.000000000000000, 0.000000000000000, 0.000000000000000, 0.000000000000000 }};

        /// <summary>
        /// Actually since I am taking LQ(transpose(A)), this is the transpose of <see cref="qrFactorQ"/>
        /// </summary>
        public static readonly double[,] lqFactorQ = {
            { -0.23075142099349710, -0.43969785029491520, -0.13419534684195095, -0.10281998694624772, -0.47126876942772489, -0.10662012686283072, -0.13251189634630042, -0.48973120246710766, -0.21006248450003742, -0.43572334812147162 },
            {  0.14858764284193901, -0.10772975041878161,  0.13124334022410436,  0.36775989343103105,  0.30986638208644840,  0.25272112351608017,  0.63811348815517166, -0.18714506899997330, -0.04046691405280455, -0.45837653401074380 },
            {  0.31792868485316900,  0.07417564455815109, -0.36103200437303257, -0.18811764198214728,  0.50635053078230396, -0.14061937078687231, -0.17161054696409667, -0.23356674775948780, -0.60620280595936349,  0.00607081461594355 },
            {  0.52033887561493763, -0.25365561595287567, -0.19632113866951575,  0.11534039321239833, -0.41979393109143270,  0.51532130006775501,  0.07128836704429709,  0.03581993650560919, -0.17752093128642521,  0.36523844289427687 },
            {  0.25119640190154863,  0.00742163793160727,  0.62338833887492229,  0.32455249539684250, -0.20998081083827991, -0.40692322124504376, -0.09596999203053819,  0.14709353964540772, -0.44985545375115621, -0.00167830105033348 },
            { -0.31609052075075511,  0.50045187775446243, -0.15238342744727723,  0.52520662918835392, -0.09764449684237536,  0.38027807781600009, -0.34460294522294443, -0.03588951204723363, -0.23659125205442419, -0.14286866951346483 },
            { -0.14778642143763676,  0.46428134845573432, -0.24643597223954369, -0.15271034272843115, -0.37794948782490234, -0.26111733434638934,  0.62009957715162045, -0.01271828595741990, -0.25089634999754107,  0.14102683949598455 },
            { -0.20713121727488037, -0.34723518223729494, -0.33041361660005214,  0.07115730533531130,  0.03273127635904442, -0.05636543553210413,  0.06108149541389936,  0.77029706318046287, -0.25248284816875666, -0.23917224941218027 },
            {  0.38894313590494700,  0.05832008220679226, -0.45550611799719420,  0.45778828012345169, -0.10261133630329747, -0.47734907759472461, -0.09053444009338751, -0.03789988799293117,  0.40429678229270199, -0.12956106417030591 },
            { -0.41690588846886784, -0.36766705946723366, -0.08159225990019033,  0.43169295990059880,  0.20297325141250672, -0.17324820800576807,  0.12769970133355366, -0.22676704045448226, -0.08211487942414421,  0.59355527534776520 }};

        /// <summary>
        /// This vector is not in the column space of the matrix. Thus the linear system can only be solved approximately, 
        /// via least squares.
        /// </summary>
        public static readonly double[] lsqRhs = { 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0 };

        /// <summary>
        /// Solution of the least squares system <see cref="matrix"/> * <see cref="lsqSolution"/> ~= <see cref="lsqRhs"/>.
        /// </summary>
        public static readonly double[] lsqSolution = { 1.073134012122104, -0.126789906548218, 1.420147230082810, 1.877832570455326, -1.117439931920276 };

        public static readonly double[] minNormRhs = { 1.0, 2.0, 3.0, 4.0, 5.0 };

        public static readonly double[] minNormSolution = { 0.4115012354314005, -0.1652228399115381, 0.3318849071280361, 0.4048451728664134, -0.4209082824823209, 0.3061133109183448, 0.2615025913792203, 0.1943785674066800, -0.0675278461506491, 0.1798101026385349 };



        public static void CheckFactorizationLQ()
        {
            var comparer = new Comparer(Comparer.PrintMode.Always);
            var A = Matrix.CreateFromArray(matrix).Transpose();
            LQFactorization LQ = A.FactorLQ();
            Matrix L = LQ.GetFactorL();
            Matrix Q = LQ.GetFactorQ();

            Console.WriteLine("Check LQ factorization: ");
            Console.WriteLine("L: ");
            comparer.CheckMatrixEquality(lqFactorL, L.CopyToArray2D());
            Console.WriteLine("Q: ");
            comparer.CheckMatrixEquality(lqFactorQ, Q.CopyToArray2D());
        }

        public static void CheckFactorizationQR()
        {
            var comparer = new Comparer(Comparer.PrintMode.Always);
            var A = Matrix.CreateFromArray(matrix);
            QRFactorization QR = A.FactorQR();
            Matrix Q = QR.GetFactorQ();
            Matrix R = QR.GetFactorR();
            
            Console.WriteLine("Check QR factorization: ");
            Console.WriteLine("Q: ");
            comparer.CheckMatrixEquality(qrFactorQ, Q.CopyToArray2D());
            Console.WriteLine("R: ");
            comparer.CheckMatrixEquality(qrFactorR, R.CopyToArray2D());
        }

        public static void CheckMatrixVectorMult()
        {
            var comparer = new Comparer(Comparer.PrintMode.Always);
            var A = Matrix.CreateFromArray(matrix);
            var x5 = Vector.CreateFromArray(lhs5);
            var x10 = Vector.CreateFromArray(lhs10);

            Vector b5 = A.MultiplyRight(x5, false);
            Vector b10 = A.MultiplyRight(x10, true);
            comparer.CheckMatrixVectorMult(matrix, lhs5, rhs5, b5.InternalData);
            comparer.CheckMatrixVectorMult(matrix, lhs10, rhs10, b10.InternalData);
        }

        public static void CheckSolutionLeastSquares()
        {
            var comparer = new Comparer(Comparer.PrintMode.Always);
            var A = Matrix.CreateFromArray(matrix);
            var QR = A.FactorQR();

            // RHS is in the column space
            var b1 = Vector.CreateFromArray(rhs5);
            var x1 = QR.SolveLeastSquares(b1);
            Console.WriteLine("Check least squares solution for a rhs vector in the column space: ");
            comparer.CheckVectorEquality(Vector.CreateFromArray(lhs5), x1);

            // RHS is not in the column space
            var b2 = Vector.CreateFromArray(lsqRhs);
            var x2 = QR.SolveLeastSquares(b2);
            Console.WriteLine("\nCheck least squares solution for a rhs vector outside the column space: ");
            comparer.CheckVectorEquality(Vector.CreateFromArray(lsqSolution), x2);
        }

        public static void CheckSolutionMinNorm()
        {
            var comparer = new Comparer(Comparer.PrintMode.Always);
            var A = Matrix.CreateFromArray(matrix).Transpose();
            var LQ = A.FactorLQ();

            var b = Vector.CreateFromArray(minNormRhs);
            var x = LQ.SolveMinNorm(b);
            Console.WriteLine("Check min norm solution: ");
            comparer.CheckVectorEquality(Vector.CreateFromArray(minNormSolution), x);
        }

        public static void Print()
        {
            var A = Matrix.CreateFromArray(matrix);
            Console.WriteLine("Rectangular matrix = ");
            var writer = new FullMatrixWriter(A);
            writer.WriteToConsole();
        }
    }
}
