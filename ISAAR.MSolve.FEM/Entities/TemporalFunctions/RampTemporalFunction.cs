using System;
using System.Collections.Generic;
using System.Text;

namespace ISAAR.MSolve.FEM.Entities.TemporalFunctions
{
    public class RampTemporalFunction : ITemporalFunction
    {
        private readonly double maxValue;
        private readonly double totalDuration;
        private readonly double timeStepDuration;
        private readonly double constantSegmentStart;

        //TODO: Defining the total duration and timestep duration follows the logic from NewmarkAnalyzer.
        public RampTemporalFunction(double maxValue, double totalDuration, double timeStepDuration, double constantSegmentStart)
        {
            this.maxValue = maxValue;
            this.totalDuration = totalDuration;
            this.timeStepDuration = timeStepDuration;
            this.constantSegmentStart = constantSegmentStart;
        }

        public double CalculateValueAt(int timeStep)
        {
            double time = (timeStep+1) * timeStepDuration; //Since it follows the logic from Newmark it goes from 0---(numofsteps-1).
            if (time <= constantSegmentStart)
            {
                return maxValue * time / constantSegmentStart;
            }
            else
            {
                return maxValue;
            }
        }
    }
}
