using ISAAR.MSolve.LinearAlgebra.Matrices;

namespace ISAAR.MSolve.Materials.Interfaces
{
    public interface IContinuumMaterial3D : IFiniteElementMaterial
    {
        double[] Stresses { get; }
        IMatrixView ConstitutiveMatrix { get; }
        bool hasfailed { get; set; }

        void UpdateMaterial(double[] strains);
    }
}
