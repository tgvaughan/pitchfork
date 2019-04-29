/*
 * Copyright (C) 2019. Tim Vaughan
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <https://www.gnu.org/licenses/>.
 */

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
        treeInput.setRule(Input.Validate.OPTIONAL);
    }

    @Override
    public void initAndValidate() {
        lcModel = lcModelInput.get();
        collapsedTreeIntervals = collapsedTreeIntervalsInput.get();
        populationFunction = populationFunctionInput.get();

        treeInput.setValue(collapsedTreeIntervals.treeInput.get(), this);
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
