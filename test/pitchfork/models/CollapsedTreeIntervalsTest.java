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

import beast.base.evolution.tree.IntervalType;
import pitchfork.PitchforkTestClass;
import org.junit.Assert;
import org.junit.Test;

public class CollapsedTreeIntervalsTest extends PitchforkTestClass {

    @Test
    public void testIntervalCalculation() {

        CollapsedTreeIntervals treeIntervals = new CollapsedTreeIntervals();
        treeIntervals.initByName("tree", tree);

        Assert.assertTrue(treeIntervals.getIntervalCount() == 3);
        Assert.assertTrue((treeIntervals.getSampleCount() == 4));
        Assert.assertTrue(treeIntervals.getIntervalType(0) == IntervalType.SAMPLE);
        Assert.assertTrue(treeIntervals.getIntervalType(1) == IntervalType.COALESCENT);
        Assert.assertTrue(treeIntervals.getIntervalType(2) == IntervalType.COALESCENT);
        Assert.assertEquals(treeIntervals.getInterval(0), 0.8, 1e-14);
        Assert.assertEquals(treeIntervals.getInterval(1), 0.2, 1e-14);
        Assert.assertEquals(treeIntervals.getInterval(2), 0.5, 1e-14);
    }
}
