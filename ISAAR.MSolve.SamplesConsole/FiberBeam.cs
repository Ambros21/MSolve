using ISAAR.MSolve.PreProcessor;
using ISAAR.MSolve.PreProcessor.Elements;
using ISAAR.MSolve.PreProcessor.Materials;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.IO;

namespace ISAAR.MSolve.SamplesConsole
{
    class FiberBeam
    {
        public static void MakeFiberBeamModel(Model model)
        {
            int fibers = 4*561; //Here as fibers are considered the number of integration points in sections multiplied by 4 GLPoints. The row of the fibers sectionsare 
            // as follows. ksi=-1 ksi=-0.4 ksi=0.4 ksi=1
            SteelFiberElementMaterial material = new SteelFiberElementMaterial(fibers, 210000000, 0.3, 0.1, 275000, 275000, -275000);
            //Here it is an example of a cantilever beam with 2 elements of 5 m Each. The rest can be easily seen in the code.
            model.NodesDictionary.Add(1, new Node() { ID = 1, X = 0, Y = 0, Z = 0 });
            model.NodesDictionary.Add(2, new Node() { ID = 2, X = 5, Y = 0, Z = 0 });
            model.NodesDictionary.Add(3, new Node() { ID = 3, X = 10, Y = 0, Z = 0 });

            foreach (Node node in model.NodesDictionary.Values)
            {
                node.Constraints.Add(DOFType.Z);
                node.Constraints.Add(DOFType.RotX);
                node.Constraints.Add(DOFType.RotY);
            }
            model.NodesDictionary[1].Constraints.Add(DOFType.X);
            model.NodesDictionary[1].Constraints.Add(DOFType.Y);
            model.NodesDictionary[1].Constraints.Add(DOFType.RotZ);
            double b = 0.25;
            double h = 0.25;
            Element e;
            e = new Element()
            {
                ID = 1,
            ElementType = new FiberBeam3D(material, fibers, b, h)
            };
            e.NodesDictionary.Add(1, model.NodesDictionary[1]);
            e.NodesDictionary.Add(2, model.NodesDictionary[2]);
            model.ElementsDictionary.Add(e.ID, e);
            model.SubdomainsDictionary[1].ElementsDictionary.Add(e.ID, e);
            SteelFiberElementMaterial material1 = new SteelFiberElementMaterial(fibers, 210000000, 0.3, 0.1, 275000, 275000, -275000);
            //It is essential to define NEW material. :)
            e = new Element()
            {
                ID = 2,
                ElementType = new FiberBeam3D(material1,fibers,b,h)
            };
            e.NodesDictionary.Add(1, model.NodesDictionary[2]);
            e.NodesDictionary.Add(2, model.NodesDictionary[3]);
            model.ElementsDictionary.Add(e.ID, e);
            model.SubdomainsDictionary[1].ElementsDictionary.Add(e.ID, e);
            model.Loads.Add(new Load() { Node = model.NodesDictionary[3], DOF = DOFType.Y, Amount = -150 });
        }
        public static SteelFiberElementMaterial[] CreateMaterials(int fibers)
        {
            SteelFiberElementMaterial[] material = new SteelFiberElementMaterial[4];
            for (int jj = 0; jj < 4; jj++)
            {
                material[jj] = new SteelFiberElementMaterial(fibers, 210000000, 0.3, 0.1, 275000, 275000, -275000);
            }
            return material;
        }
    }
}
