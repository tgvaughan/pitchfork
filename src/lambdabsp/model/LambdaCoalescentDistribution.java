package lambdabsp.model;

import beast.core.Input;
import beast.core.parameter.RealParameter;
import beast.evolution.tree.TreeDistribution;
import beast.evolution.tree.coalescent.PopulationFunction;

public class LambdaCoalescentDistribution extends TreeDistribution {

    public Input<LBSPTreeIntervals> lambdaTreeIntervalsInput = new Input<>(
            "treeIntervals",
            "Lambda-coalescent tree intervals object.",
            Input.Validate.REQUIRED);

    public Input<RealParameter> alphaInput = new Input<>(
            "alpha",
            "Alpha parameter for Beta-coalescent process",
            Input.Validate.REQUIRED);

    public Input<PopulationFunction> populationFunctionInput = new Input<>(
            "populationFunction",
            "Population function object.",
            Input.Validate.REQUIRED);

    @Override
    public void initAndValidate() {

    }

    @Override
    public double calculateLogP() {
        logP = 0.0;

        return logP;
    }
}
