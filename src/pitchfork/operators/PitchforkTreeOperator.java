package pitchfork.operators;

import beast.core.Input;
import beast.core.parameter.RealParameter;
import beast.evolution.operators.TreeOperator;
import beast.util.Randomizer;
import pitchfork.Pitchforks;
import pitchfork.models.pop.SkylinePopulationFunction;

/**
 * Class of operators for traversing the space of multifurcating phylogenetic trees.
 */
public abstract class PitchforkTreeOperator extends TreeOperator {

    public Input<SkylinePopulationFunction> skylineInput = new Input<>("skylinePopFun",
            "Skyline population function. Required for skyline demographic models.");

    boolean isSkylineModel;
    SkylinePopulationFunction skyline;

    @Override
    public void initAndValidate() {
        skyline = skylineInput.get();
        isSkylineModel = (skyline != null);
    }

    @Override
    public final double proposal() {
        if (isSkylineModel)
            return pitchforkProposal();

        int initialIntervalCount = skyline.getTrueSkylineIntervalCount();

        double logHR = pitchforkProposal();

        if (isSkylineSafe() || logHR == Double.NEGATIVE_INFINITY)
            return logHR;

        // TODO Implement skyline function update.

        int finalIntervalCount = skyline.getTrueSkylineIntervalCount();

        RealParameter popSizeParameter = skyline.popSizesInput.get();

        if (finalIntervalCount > initialIntervalCount) {

            for (int i=initialIntervalCount; i<finalIntervalCount; i++) {
                popSizeParameter.setValue(i,
                        Randomizer.nextExponential(1.0/popSizeParameter.getValue(0)));

                logHR -= -popSizeParameter.getValue(i)/popSizeParameter.getValue(0)
                        + Math.log(1.0/popSizeParameter.getValue(0));
            }

        } else if (finalIntervalCount < initialIntervalCount) {
            for (int i=finalIntervalCount; i<initialIntervalCount; i++) {
                logHR += -popSizeParameter.getValue(i)/popSizeParameter.getValue(0)
                        + Math.log(1.0/popSizeParameter.getValue(0));
            }
        }

        return logHR;
    }

    protected abstract double pitchforkProposal();

    abstract boolean isSkylineSafe();
}
