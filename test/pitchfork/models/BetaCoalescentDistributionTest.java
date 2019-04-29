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

package pitchfork.models;

import beast.core.parameter.RealParameter;
import junit.framework.Assert;
import pitchfork.PitchforkTestClass;
import org.junit.Test;

public class BetaCoalescentDistributionTest extends PitchforkTestClass {

    @Test
    public void testDistribution() {
        BetaCoalescentModel bcModel = new BetaCoalescentModel();
        bcModel.initByName("alpha", new RealParameter("1.0"),
                "tree", tree);

        CollapsedTreeIntervals treeIntervals = new CollapsedTreeIntervals();
        treeIntervals.initByName("tree", tree);

        BetaCoalescentDistribution distribution = new BetaCoalescentDistribution();
        distribution.initByName("model", bcModel,
                "collapsedTreeIntervals", treeIntervals,
                "populationFunction", getConstantPopulation(1.0));

        double density = distribution.calculateLogP();
        System.out.println("Density: " + density);

        // TODO check this is actually right!
        Assert.assertEquals(-4.49175946928192, density, 1e-10);
    }
}
