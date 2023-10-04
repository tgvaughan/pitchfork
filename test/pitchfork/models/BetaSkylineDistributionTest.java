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

        AbstractBetaSkylineDistribution distribution = new AbstractBetaSkylineDistribution();
        distribution.initByName("model", bcModel,
                "collapsedTreeIntervals", treeIntervals,
                "groupSizes", new IntegerParameter("1 1 1"), //to do check total size
                "skylinePopulations", new RealParameter("1.0 1.0 1.0"));  //update

        double density = distribution.calculateLogP();
        System.out.println("Density: " + density);

        // TODO check this is actually right!
        Assert.assertEquals(-3.10546510810816, density, 1e-10);  //before -4.49175946928192
    }


   @Test
    public void testDistributionBig() {
        BetaCoalescentModel bcModel = new BetaCoalescentModel();
        bcModel.initByName("alpha", new RealParameter("1.0"),
                "tree", tree);

        CollapsedTreeIntervals treeIntervals = new CollapsedTreeIntervals();
        treeIntervals.initByName("tree", tree);

        AbstractBetaSkylineDistribution distribution = new AbstractBetaSkylineDistribution();
        distribution.initByName("model", bcModel,
                "collapsedTreeIntervals", treeIntervals,
                "groupSizes", new IntegerParameter("1 1 1 1 1 1 1 1"), //to do check total size of 8
                "skylinePopulations", new RealParameter("1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0"));  //update

        double density = distribution.calculateLogP();
        System.out.println("Density: " + density);

        // TODO check this is actually right!
        Assert.assertEquals(-14.09056208756590, density, 1e-8);
    }
}

