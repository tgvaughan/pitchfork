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
