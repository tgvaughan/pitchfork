package pitchfork.models;

import beast.core.Input;
import beast.evolution.tree.TreeDistribution;
import beast.evolution.tree.coalescent.IntervalType;
import beast.evolution.tree.coalescent.PopulationFunction;

public class BetaCoalescentDistribution extends TreeDistribution {

    public Input<CollapsedTreeIntervals> collapsedTreeIntervalsInput = new Input<>(
            "collapsedTreeIntervals",
            "Collapsed tree intervals object.",
            Input.Validate.REQUIRED);

    public Input<BetaCoalescentModel> lcModelInput = new Input<>(
            "model",
            "Beta-coalescent model.",
            Input.Validate.REQUIRED);

    public Input<PopulationFunction> populationFunctionInput = new Input<>(
            "populationFunction",
            "Population function object.",
            Input.Validate.REQUIRED);

    private CollapsedTreeIntervals collapsedTreeIntervals;
    private BetaCoalescentModel lcModel;
    private PopulationFunction populationFunction;

    public BetaCoalescentDistribution() {
        treeIntervalsInput.setRule(Input.Validate.FORBIDDEN);
        treeInput.setRule(Input.Validate.FORBIDDEN);
    }

    @Override
    public void initAndValidate() {
        lcModel = lcModelInput.get();
        collapsedTreeIntervals = collapsedTreeIntervalsInput.get();
        populationFunction = populationFunctionInput.get();
    }

    @Override
    public double calculateLogP() {
        logP = 0.0;

        double t=0;
        for (int i = 0; i< collapsedTreeIntervals.getIntervalCount(); i++) {

            // Get interval details
            double dt = collapsedTreeIntervals.getInterval(i);
            int n = collapsedTreeIntervals.getLineageCount(i);

            // Waiting time contribution
            logP += -lcModel.getTotalCoalRate(n)*populationFunction.getIntegral(t, t+dt);

            // Increment time
            t += dt;

            if (collapsedTreeIntervals.getIntervalType(i) == IntervalType.COALESCENT) {
                // Beta-coalescent event contribution

                int k = collapsedTreeIntervals.getCoalescentEvents(i)+1;
                double N = populationFunction.getPopSize(t);

                double delta = lcModel.getLogLambda(n, k) - Math.log(N);

                logP += delta;
            }
        }

        return logP;
    }

    @Override
    protected boolean requiresRecalculation() {
        return true;
    }
}
