package pitchfork.models;

import beast.base.inference.parameter.IntegerParameter;
import beast.base.inference.parameter.RealParameter;
import junit.framework.Assert;
import org.junit.Test;
import pitchfork.PitchforkTestClass;

public class BetaSkylineDistributionTest extends PitchforkTestClass {

    @Test
    public void testDistribution() {
        BetaCoalescentModel bcModel = new BetaCoalescentModel();
        bcModel.initByName("alpha", new RealParameter("1.0"),
                "tree", tree);

        CollapsedTreeIntervals treeIntervals = new CollapsedTreeIntervals();
        treeIntervals.initByName("tree", tree);

        BetaSkylineDistribution distribution = new BetaSkylineDistribution();
        distribution.initByName("model", bcModel,
                "collapsedTreeIntervals", treeIntervals,
                "groupSizes", new IntegerParameter("1 1 1"), //to do check total size
                "skylinePopulations", new RealParameter("1.0 1.0 1.0"));  //update

        double density = distribution.calculateLogP();
        System.out.println("Density: " + density);

        // TODO check this is actually right!
        Assert.assertEquals(-4.49175946928192, density, 1e-10);
    }
}
