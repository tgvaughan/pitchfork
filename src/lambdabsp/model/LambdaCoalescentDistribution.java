package lambdabsp.model;

import beast.core.Input;
import beast.evolution.tree.TreeDistribution;
import beast.evolution.tree.coalescent.PopulationFunction;

public class LambdaCoalescentDistribution extends TreeDistribution {

    public Input<LBSPTreeIntervals> lambdaTreeIntervalsInput = new Input<>(
            "treeIntervals",
            "Lambda-coalescent tree intervals object.",
            Input.Validate.REQUIRED);

    public Input<LambdaCoalescentModel> lcModelInput = new Input<>(
            "model",
            "Lambda-coalescent model.",
            Input.Validate.REQUIRED);

    public Input<PopulationFunction> populationFunctionInput = new Input<>(
            "populationFunction",
            "Population function object.",
            Input.Validate.REQUIRED);

    LBSPTreeIntervals treeIntervals;
    LambdaCoalescentModel lcModel;

    double[][] cumulativeCoalRates;

    @Override
    public void initAndValidate() {
        treeIntervals = lambdaTreeIntervalsInput.get();
    }



    @Override
    public double calculateLogP() {
        logP = 0.0;

        return logP;
    }
}
