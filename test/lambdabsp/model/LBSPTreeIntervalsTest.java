package lambdabsp.model;

import beast.evolution.tree.coalescent.IntervalType;
import beast.util.TreeParser;
import org.junit.Assert;
import org.junit.Test;

public class LBSPTreeIntervalsTest {

    @Test
    public void testIntervalCalculation() {

        // Generate test tree.  Note that TreeParser will convert this to a binary tree by replacing multifuractions
        // with ladders containing zero-length edges.  Conversely, LBSPTreeIntervals ignores zero-length intervals.
        TreeParser tree = new TreeParser("((A:1.0,B:1.0,C:1.0):0.5,D:0.7);",
                false, false, true,0);

        LBSPTreeIntervals treeIntervals = new LBSPTreeIntervals();
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
