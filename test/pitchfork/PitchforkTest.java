package pitchfork;

import beast.core.parameter.RealParameter;
import beast.evolution.tree.coalescent.ConstantPopulation;
import beast.evolution.tree.coalescent.ExponentialGrowth;
import beast.evolution.tree.coalescent.PopulationFunction;
import beast.util.TreeParser;

/**
 * Abstract class containing methods useful in testing Pitchfork code.
 */
public abstract class PitchforkTest {

    // Generate test tree.  Note that TreeParser will convert this to a binary
    // tree by replacing multifuractions with ladders containing zero-length
    // edges.  Conversely, CollapsedTreeIntervals ignores zero-length intervals.
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
