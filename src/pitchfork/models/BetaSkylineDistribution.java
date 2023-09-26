/*
 * Copyright (C) 2020. Tim Vaughan
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

import beast.base.core.Input;
import beast.base.evolution.tree.IntervalType;
import beast.base.evolution.tree.Tree;
import beast.base.evolution.tree.TreeDistribution;
import beast.base.inference.parameter.IntegerParameter;
import beast.base.inference.parameter.RealParameter;
import pitchfork.Pitchforks;

public class BetaSkylineDistribution extends TreeDistribution {

    public Input<CollapsedTreeIntervals> collapsedTreeIntervalsInput = new Input<>(
            "collapsedTreeIntervals",
            "Collapsed tree intervals object.",
            Input.Validate.REQUIRED);

    public Input<BetaCoalescentModel> betaCoalescentModelInput = new Input<>(
            "model",
            "Beta-coalescent model.",
            Input.Validate.REQUIRED);

    public Input<RealParameter> skylinePopulationsInput = new Input<>(
            "skylinePopulations",
            "Parameter containing skyline population sizes.",
            Input.Validate.REQUIRED);

    public Input<IntegerParameter> groupSizesInput = new Input<>(
            "groupSizes",
            "Parameter containing skyline group sizes.",
            Input.Validate.REQUIRED);

    public BetaSkylineDistribution() {
        treeIntervalsInput.setRule(Input.Validate.FORBIDDEN);
        treeInput.setRule(Input.Validate.OPTIONAL);
    }


    private CollapsedTreeIntervals collapsedTreeIntervals;
    private BetaCoalescentModel betaCoalescentModel;
    private RealParameter populationSizes;
    private IntegerParameter groupSizes;
    private Tree tree;

    @Override
    public void initAndValidate() {
        betaCoalescentModel = betaCoalescentModelInput.get();
        collapsedTreeIntervals = collapsedTreeIntervalsInput.get();
        populationSizes = skylinePopulationsInput.get();
        groupSizes = groupSizesInput.get();

        tree = collapsedTreeIntervals.treeInput.get();
        treeInput.setValue(tree, this);

        int groupSizesSum = 0;
        for (int groupSize : groupSizes.getValues())
            groupSizesSum += groupSize;

        if (groupSizesSum != tree.getLeafNodeCount() - 1)
            throw new IllegalArgumentException("Sum of elements of groupSizes input should equal number of leaves - 1.");
    }

    @Override
    public double calculateLogP() {
        logP = 0.0;

        double t = 0;
        double N = skylinePopulationsInput.get().getValue(0);
        int seenCoalescentEvents = 0;
        int group = 0;
        int totalCoalecents = 0;  //this variable is needed that N is not updated after the last coalescent events because there is no N value left
        for (int i = 0; i < collapsedTreeIntervals.getIntervalCount(); i++) {
                totalCoalecents += collapsedTreeIntervals.getCoalescentEvents(i);
        }
        int counter = 0;     //this variable counts all the coalescent events and is compared to totalCoalescents

        for (int i = 0; i < collapsedTreeIntervals.getIntervalCount(); i++) {
            // Get interval details
            double dt = collapsedTreeIntervals.getInterval(i);
            int n = collapsedTreeIntervals.getLineageCount(i);

            // Waiting time contribution
            logP += -betaCoalescentModel.getTotalCoalRate(n)*dt/N;

            // Increment time
            t += dt;

            if (collapsedTreeIntervals.getIntervalType(i) == IntervalType.COALESCENT) {
                // Beta-coalescent event contribution

                int k = collapsedTreeIntervals.getCoalescentEvents(i)+1;
                logP += betaCoalescentModel.getLogLambda(n, k) - Math.log(N);

                // while loop to update N until the seen coalescent events are smaller or equal than the actual group size
                seenCoalescentEvents += collapsedTreeIntervals.getCoalescentEvents(i);
                counter += collapsedTreeIntervals.getCoalescentEvents(i);

                if (counter < totalCoalecents) {           //because we don't want to update N after the last coalescent event
                    while (seenCoalescentEvents >= groupSizesInput.get().getValue(group)) {
                        seenCoalescentEvents -= groupSizesInput.get().getValue(group);
                        group += 1;
                    }
                    N = skylinePopulationsInput.get().getValue(group);
                }
                //logP += betaCoalescentModel.getLogLambda(n, k) - Math.log(N);
            }
        }
        return logP;
    }
}
