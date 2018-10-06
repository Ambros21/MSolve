using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using ISAAR.MSolve.PreProcessor.Interfaces;

namespace ISAAR.MSolve.PreProcessor.Elements
{
    public class Fiber : IFiber
    {
        private readonly double b, h, x, y;

        public Fiber(double b, double h, double x, double y)
        {
            this.b = b;
            this.h = h;
            this.x = x;
            this.y = y;
        }

        public IFiberMaterial Material { get; set; }
        public double B
        {
            get { return b; }
        }

        public double H
        {
            get { return h; }
        }

        public double X
        {
            get { return x; }
        }

        public double Y
        {
            get { return y; }
        }
    }
}
