package pitchfork.models;

import beast.core.parameter.RealParameter;
import junit.framework.Assert;
import pitchfork.LBSPTest;
import org.junit.Test;

public class LambdaCoalescentDistributionTest extends LBSPTest {

    @Test
    public void testDistribution() {
        LambdaCoalescentModel lcModel = new LambdaCoalescentModel();
        lcModel.initByName("alpha", new RealParameter("1.0"),
                "maxExtantLineages", tree.getLeafNodeCount());

        CollapsedTreeIntervals treeIntervals = new CollapsedTreeIntervals();
        treeIntervals.initByName("tree", tree);

        LambdaCoalescentDistribution distribution = new LambdaCoalescentDistribution();
        distribution.initByName("model", lcModel,
                "collapsedTreeIntervals", treeIntervals,
                "populationFunction", getConstantPopulation(1.0));

        double density = distribution.calculateLogP();
        System.out.println("Density: " + density);

        // TODO check this is actually right!
        Assert.assertEquals(-4.49175946928192, density, 1e-10);
    }
}
