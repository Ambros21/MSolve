﻿using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ISAAR.MSolve.Numerical.LinearAlgebra;
using ISAAR.MSolve.XFEM.CrackGeometry;
using ISAAR.MSolve.XFEM.Elements;
using ISAAR.MSolve.XFEM.Entities;
using ISAAR.MSolve.XFEM.Entities.FreedomDegrees;
using ISAAR.MSolve.XFEM.Geometry.CoordinateSystems;
using ISAAR.MSolve.XFEM.Geometry.Shapes;
using ISAAR.MSolve.XFEM.Geometry.Mesh;
using ISAAR.MSolve.XFEM.CrackPropagation.Direction;
using ISAAR.MSolve.XFEM.CrackPropagation.Jintegral;
using ISAAR.MSolve.XFEM.CrackPropagation.Length;
using ISAAR.MSolve.XFEM.Integration.Points;
using ISAAR.MSolve.XFEM.Integration.Strategies;
using ISAAR.MSolve.XFEM.Interpolation;
using ISAAR.MSolve.XFEM.LinearAlgebra;
using ISAAR.MSolve.XFEM.Tensors;
using ISAAR.MSolve.XFEM.Utilities;


namespace ISAAR.MSolve.XFEM.CrackPropagation
{
    class Propagator
    {
        private readonly IMesh2D<XNode2D, XContinuumElement2D> mesh;
        private readonly ICrackDescription crack;
        private readonly double magnificationOfJintegralRadius;
        private readonly IAuxiliaryStates auxiliaryStatesStrategy;
        private readonly ISIFCalculator sifCalculationStrategy;
        private readonly ICrackGrowthDirectionLaw2D growthDirectionLaw;
        private readonly ICrackGrowthLengthLaw2D growthLengthLaw;

        public PropagationLogger Logger { get; }

        /// <summary>
        /// 
        /// </summary>
        /// <param name="mesh"></param>
        /// <param name="crack"></param>
        /// <param name="magnificationOfJintegralRadius">The outer countour of the J-integral domain is defined as:
        ///     radius = magnification * sqrt(areaOfElementContainingTip). This parameter is the magnification. 
        ///     It should be at least 1.5 (see "Modeling quasi-static crack growth with the extended finite element 
        ///     method Part II: Numerical applications, Huang et al, 2003" page 7546). Usually values 2-3 are selected 
        ///     (see Ahmed thesis, 2009).</param>
        /// <param name="auxiliaryStatesStrategy"></param>
        /// <param name="sifCalculationStrategy"></param>
        public Propagator(IMesh2D<XNode2D, XContinuumElement2D> mesh, ICrackDescription crack,
            double magnificationOfJintegralRadius,
            IAuxiliaryStates auxiliaryStatesStrategy, ISIFCalculator sifCalculationStrategy,
            ICrackGrowthDirectionLaw2D growthDirectionLaw, ICrackGrowthLengthLaw2D growthLengthLaw)
        {
            this.mesh = mesh;
            this.crack = crack;
            this.magnificationOfJintegralRadius = magnificationOfJintegralRadius;
            this.auxiliaryStatesStrategy = auxiliaryStatesStrategy;
            this.sifCalculationStrategy = sifCalculationStrategy;
            this.growthDirectionLaw = growthDirectionLaw;
            this.growthLengthLaw = growthLengthLaw;
            this.Logger = new PropagationLogger();
        }

        public void Propagate(Model2D model, double[] totalFreeDisplacements, double[] totalConstrainedDisplacements,
            out double growthAngle, out double growthLength)
        {
            // TODO: Also check if the sifs do not violate the material toughness
            double sifMode1, sifMode2;
            ComputeSIFS(model, totalFreeDisplacements, totalConstrainedDisplacements, out sifMode1, out sifMode2);
            growthAngle = growthDirectionLaw.ComputeGrowthAngle(sifMode1, sifMode2);
            growthLength = growthLengthLaw.ComputeGrowthLength(sifMode1, sifMode2);
            Logger.GrowthAngles.Add(growthAngle);
            Logger.GrowthLengths.Add(growthLength);
        }

        private void ComputeSIFS(Model2D model, double[] totalFreeDisplacements, double[] totalConstrainedDisplacements,
             out double sifMode1, out double sifMode2)
        {
            double interactionIntegralMode1 = 0.0, interactionIntegralMode2 = 0.0;
            IReadOnlyDictionary<XContinuumElement2D, double[]> elementWeights = FindJintegralElementsAndNodalWeights();
            foreach (var pair in elementWeights)
            {
                XContinuumElement2D element = pair.Key;
                double[] nodalWeights = pair.Value;
                double[] standardElementDisplacements = model.DofEnumerator.ExtractDisplacementVectorOfElementFromGlobal(
                    element, totalFreeDisplacements, totalConstrainedDisplacements);
                double[] enrichedElementDisplacements = model.DofEnumerator.
                    ExtractEnrichedDisplacementsOfElementFromGlobal(element, totalFreeDisplacements);

                double partialIntegralMode1, partialIntegralMode2;
                ComputeInteractionIntegrals(element, standardElementDisplacements, enrichedElementDisplacements,
                    nodalWeights, out partialIntegralMode1, out partialIntegralMode2);

                interactionIntegralMode1 += partialIntegralMode1;
                interactionIntegralMode2 += partialIntegralMode2;
            }

            sifMode1 = sifCalculationStrategy.CalculateSIF(interactionIntegralMode1);
            sifMode2 = sifCalculationStrategy.CalculateSIF(interactionIntegralMode2);

            Logger.InteractionIntegralsMode1.Add(interactionIntegralMode1);
            Logger.InteractionIntegralsMode2.Add(interactionIntegralMode2);
            Logger.SIFsMode1.Add(sifMode1);
            Logger.SIFsMode2.Add(sifMode2);
        }

        private IReadOnlyDictionary<XContinuumElement2D, double[]> FindJintegralElementsAndNodalWeights()
        {
            Circle2D outerContour = new Circle2D(crack.CrackTip, ComputeRadiusOfJintegralOuterContour());
            IReadOnlyList<XContinuumElement2D> intersectedElements =
                mesh.FindElementsIntersectedByCircle(outerContour, crack.TipElements[0]);

            var elementsAndWeights = new Dictionary<XContinuumElement2D, double[]>();
            foreach (var element in intersectedElements)
            {
                // The relative position of the circle and the nodes was already calculated when checking the
                // circle-element intersection, but that method should be decoupled from assigning the nodal 
                // weights, even at the cost of some duplicate operations. What could be done more efficiently is 
                // caching the nodes and weights already processed by previous elements, but even then the cost of
                // processing each node will be increased by the lookup.
                double[] nodalWeights = new double[element.Nodes.Count];
                for (int nodeIdx = 0; nodeIdx < element.Nodes.Count; ++nodeIdx)
                {
                    if (outerContour.FindRelativePositionOfPoint(element.Nodes[nodeIdx]) == CirclePointPosition.Outside)
                    {
                        nodalWeights[nodeIdx] = 0.0;
                    }
                    else // Node lies inside or exactly on the circle
                    {
                        nodalWeights[nodeIdx] = 1.0;
                    }
                }
                elementsAndWeights.Add(element, nodalWeights);
            }
            return elementsAndWeights;
        }

        // This method should directly return the elements and take care of cases near the domain boundaries (see Ahmed)
        public double ComputeRadiusOfJintegralOuterContour()
        {
            double maxTipElementArea = -1.0;
            foreach (var element in crack.TipElements)
            {
                var outline = ConvexPolygon2D.CreateUnsafe(element.Nodes);
                double elementArea = outline.ComputeArea();
                if (elementArea > maxTipElementArea) maxTipElementArea = elementArea;
            }
            return magnificationOfJintegralRadius * Math.Sqrt(maxTipElementArea);
        }

        private void ComputeInteractionIntegrals(XContinuumElement2D element, double[] standardNodalDisplacements,
            double[] enrichedNodalDisplacements, double[] nodalWeights, 
            out double integralMode1, out double integralMode2)
        {
            integralMode1 = 0.0;
            integralMode2 = 0.0;

            foreach (GaussPoint2D naturalGP in element.JintegralStrategy.GenerateIntegrationPoints(element))
            {
                // Nomenclature: global = global cartesian system, natural = element natural system, 
                // local = tip local cartesian system  
                EvaluatedInterpolation2D evaluatedInterpolation = 
                    element.Interpolation.EvaluateAt(element.Nodes, naturalGP);
                ICartesianPoint2D globalGP = evaluatedInterpolation.TransformPointNaturalToGlobalCartesian(naturalGP);
                Matrix2D constitutive = 
                    element.Material.CalculateConstitutiveMatrixAt(naturalGP, evaluatedInterpolation);

                // State 1
                DenseMatrix globalDisplacementGradState1 = element.CalculateDisplacementFieldGradient(
                    naturalGP, evaluatedInterpolation, standardNodalDisplacements, enrichedNodalDisplacements);
                Tensor2D globalStressState1 = element.CalculateStressTensor(globalDisplacementGradState1, constitutive);
                DenseMatrix localDisplacementGradState1 = crack.TipSystem.
                    TransformVectorFieldDerivativesGlobalCartesianToLocalCartesian(globalDisplacementGradState1);
                Tensor2D localStressTensorState1 = crack.TipSystem.
                    TransformTensorGlobalCartesianToLocalCartesian(globalStressState1);

                // Weight Function
                // TODO: There should be a method InterpolateScalarGradient(double[] nodalValues) in EvaluatedInterpolation
                double[] globalWeightGradient = new double[2];
                for (int nodeIdx = 0; nodeIdx < element.Nodes.Count; ++nodeIdx)
                {
                    var interpolationGradient = 
                        evaluatedInterpolation.GetGlobalCartesianDerivativesOf(element.Nodes[nodeIdx]);
                    globalWeightGradient[0] += interpolationGradient.Item1 * nodalWeights[nodeIdx];
                    globalWeightGradient[1] += interpolationGradient.Item2 * nodalWeights[nodeIdx];
                }
                double[] localWeightGradient = crack.TipSystem.
                    TransformScalarFieldDerivativesGlobalCartesianToLocalCartesian(globalWeightGradient);

                // State 2
                // TODO: XContinuumElement shouldn't have to pass tipCoordinate system to auxiliaryStates. 
                // It would be better to have CrackTip handle this and the coordinate transformations. That would also 
                // obey LoD, but a lot of wrapper methods would be required.
                AuxiliaryStatesTensors auxiliary = auxiliaryStatesStrategy.ComputeTensorsAt(globalGP, crack.TipSystem);

                // Interaction integrals
                double integrandMode1 = ComputeJIntegrand(localWeightGradient, localDisplacementGradState1,
                    localStressTensorState1, auxiliary.DisplacementGradientMode1,
                    auxiliary.StrainTensorMode1, auxiliary.StressTensorMode1);
                double integrandMode2 = ComputeJIntegrand(localWeightGradient, localDisplacementGradState1,
                    localStressTensorState1, auxiliary.DisplacementGradientMode2,
                    auxiliary.StrainTensorMode2, auxiliary.StressTensorMode2);

                integralMode1 += integrandMode1 * evaluatedInterpolation.Jacobian.Determinant * naturalGP.Weight;
                integralMode2 += integrandMode2 * evaluatedInterpolation.Jacobian.Determinant * naturalGP.Weight;
            }
        }

        private static double ComputeJIntegrand(double[] weightGrad, DenseMatrix displGrad1, Tensor2D stress1,
            DenseMatrix displGrad2, Tensor2D strain2, Tensor2D stress2)
        {
            // Unrolled to greatly reduce mistakes. Alternatively Einstein notation products could be implementated
            // in Tensor2D (like the tensor-tensor multiplication is), but still some parts would have to be unrolled.
            // Perhaps vector (and scalar) gradients should also be accessed by component and derivative variable.

            double strainEnergy = stress1.MultiplyColon(strain2);
            double parenthesis0 = stress1.XX * displGrad2[0, 0] + stress1.XY * displGrad2[1, 0]
                - stress2.XX * displGrad1[0, 0] - stress2.XY * displGrad1[1, 0] - strainEnergy;
            double parenthesis1 = stress1.XY * displGrad2[0, 0] + stress1.YY * displGrad2[1, 0]
                - stress2.XY * displGrad1[0, 0] - stress2.YY * displGrad1[1, 0];
            return parenthesis0 * weightGrad[0] + parenthesis1 * weightGrad[1];

        }
    }
}
