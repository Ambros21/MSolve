﻿using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace ISAAR.MSolve.PreProcessor
{
    public class ElementMassAccelerationLoad
    {
        public Element Element { get; set; }
        public DOFType DOF { get; set; }
        public double Amount { get; set; }
    }
}
