/*
 * Copyright (C) 2019. Tim Vaughan
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <https://www.gnu.org/licenses/>.
 */

package pitchfork;

import beast.base.evolution.tree.TreeParser;
import beast.base.evolution.tree.coalescent.ConstantPopulation;
import beast.base.evolution.tree.coalescent.ExponentialGrowth;
import beast.base.evolution.tree.coalescent.PopulationFunction;
import beast.base.inference.parameter.RealParameter;

/**
 * Abstract class containing methods useful in testing Pitchfork code.
 */
public abstract class PitchforkTestClass {

    // Generate test tree.  Note that TreeParser will convert this to a binary
    // tree by replacing multifuractions with ladders containing zero-length
    // edges.  Conversely, CollapsedTreeIntervals ignores zero-length intervals.
    //protected TreeParser tree = new TreeParser("((A:1.0,B:1.0,C:1.0):0.5,D:0.7);",
    //        false, false, true, 0);
    //protected TreeParser tree = new TreeParser("(((((A:1,B:1):0.5,C:1.5,D:1.5,E:1.5):0.7,F:2.2):0.3,G:2.5,H:2.5):0.5,I:3);",
    //        false, false, true,0);
    protected TreeParser tree = new TreeParser("((A:1.0,B:1.0,C:1.0):0.5,D:0.7);",
            false, false, true,0);



    protected PopulationFunction getConstantPopulation(double size) {
        ConstantPopulation popFun = new ConstantPopulation();
        popFun.initByName("popSize", new RealParameter(String.valueOf(size)));

        return popFun;
    }

    protected PopulationFunction getExponentialPopulation(double N0, double lambda) {
        ExponentialGrowth popFun = new ExponentialGrowth();
        popFun.initByName("popSize", String.valueOf(N0),
                "growthRate", String.valueOf(lambda));

        return popFun;
    }

}


