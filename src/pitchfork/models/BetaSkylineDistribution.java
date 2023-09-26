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

    //function to calculate sum of an array
    public static int findSum(int[] array) {
        int sum = 0;
        for (int value : array) {
            sum += value;
        }
        return sum;
    }

    //function to calculate cumulative sum array of other array
    public static int[] cumulativeSum(int[] array) {
        int[] vec = new int[array.length];
        int sum = 0;
        for (int i = 0; i < array.length; i++){
            sum += array[i];
            vec[i] = sum;
        }

        return vec;
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

        int nCoalescentIntrvals = Pitchforks.getTrueNodes(tree).size();
    }

    /*
    @Override
    public double calculateLogP2() {
        logP = 0.0;

        double t = 0;
        double N = skylinePopulationsInput.get().getValue(0);
        int seenCoalescentEvents = 0;
        int[] groupIndicator;  //creates cumulative sum array of group sizes to find out in which group we are after coalescent event
        int[] groupIndicator = cumulativeSum(groupSizes);

        for (int i = 0; i < collapsedTreeIntervals.getIntervalCount(); i++) {
            // Get interval details
            double dt = collapsedTreeIntervals.getInterval(i);
            int n = collapsedTreeIntervals.getLineageCount(i);

            // Waiting time contribution
            logP += -betaCoalescentModel.getTotalCoalRate(n)*N;

            // Increment time
            t += dt;

            if (collapsedTreeIntervals.getIntervalType(i) == IntervalType.COALESCENT) {
                // Beta-coalescent event contribution

                int k = collapsedTreeIntervals.getCoalescentEvents(i)+1;
                logP += betaCoalescentModel.getLogLambda(n, k) - Math.log(N);

                // while loop to update N until the cumulative sum of group sizes in group size vector is reached
                seenCoalescentEvents += collapsedTreeIntervals.getCoalescentEvents(i);

                for (int j = 0; j < groupIndicator.length; j++){
                    if (seenCoalescentEvents < groupIndicator[j]){
                        N = skylinePopulationsInput.get().getValue(j);
                    }
                }


            }
        }
        return logP;
    }

     */

    @Override
    public double calculateLogP() {
        logP = 0.0;

        double t = 0;
        double N = skylinePopulationsInput.get().getValue(0);
        int seenCoalescentEvents = 0;
        int group = 0;

        for (int i = 0; i < collapsedTreeIntervals.getIntervalCount(); i++) {
            // Get interval details
            double dt = collapsedTreeIntervals.getInterval(i);
            int n = collapsedTreeIntervals.getLineageCount(i);

            // Waiting time contribution
            logP += -betaCoalescentModel.getTotalCoalRate(n)*N;

            // Increment time
            t += dt;

            if (collapsedTreeIntervals.getIntervalType(i) == IntervalType.COALESCENT) {
                // Beta-coalescent event contribution

                int k = collapsedTreeIntervals.getCoalescentEvents(i)+1;
                logP += betaCoalescentModel.getLogLambda(n, k) - Math.log(N);

                // while loop to update N until the seen coalescent events are smaller or equal than the actual group size
                seenCoalescentEvents += collapsedTreeIntervals.getCoalescentEvents(i);

                while (seenCoalescentEvents > groupSizesInput.get().getValue(group)){
                    seenCoalescentEvents -= groupSizesInput.get().getValue(group);
                    group += 1;
                }
                N = skylinePopulationsInput.get().getValue(group);

            }
        }
        return logP;
    }
}
