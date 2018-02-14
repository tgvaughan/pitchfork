package pitchfork.models;

import beast.evolution.tree.coalescent.IntervalType;
import pitchfork.LBSPTest;
import org.junit.Assert;
import org.junit.Test;

public class CollapsedTreeIntervalsTest extends LBSPTest {

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
