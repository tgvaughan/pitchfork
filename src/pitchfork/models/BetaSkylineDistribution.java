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

    public BetaSkylineDistribution() {
        treeIntervalsInput.setRule(Input.Validate.FORBIDDEN);
        treeInput.setRule(Input.Validate.OPTIONAL);
    }

    private CollapsedTreeIntervals collapsedTreeIntervals;
    private BetaCoalescentModel betaCoalescentModel;
    private RealParameter skylinePopulations;
    private Tree tree;

    private int mmin;

    @Override
    public void initAndValidate() {
        betaCoalescentModel = betaCoalescentModelInput.get();
        collapsedTreeIntervals = collapsedTreeIntervalsInput.get();
        skylinePopulations = skylinePopulationsInput.get();

        tree = collapsedTreeIntervals.treeInput.get();
        treeInput.setValue(tree, this);

        int nCoalescentIntrvals = Pitchforks.getTrueNodes(tree).size();

        mmin = nCoalescentIntrvals/skylinePopulations.getDimension();
        if (nCoalescentIntrvals % skylinePopulationsInput.get().getDimension() > 0)
            mmin += 1;
    }

    @Override
    public double calculateLogP() {
        logP = 0.0;

        int nCoalescentIntervals = Pitchforks.getTrueInternalNodes(tree).size();

        int groupIdx=0, coalIntervalIdx = 0;

        for (int i=0; i<collapsedTreeIntervals.getIntervalCount(); i++) {

            double dt = collapsedTreeIntervals.getInterval(i);
            int n = collapsedTreeIntervals.getLineageCount(i);

            // Waiting time contribution
            logP += -dt*betaCoalescentModel.getTotalCoalRate(n)/skylinePopulations.getArrayValue(groupIdx);

            if (collapsedTreeIntervals.getIntervalType(i) == IntervalType.COALESCENT) {
                // Beta-coalescent event contribution
                int k = collapsedTreeIntervals.getCoalescentEvents(i)+1;
                logP += betaCoalescentModel.getLogLambda(n, k) - Math.log(skylinePopulations.getArrayValue(groupIdx));

                // Skyline interval adjustment
                if (coalIntervalIdx / mmin > groupIdx && (nCoalescentIntervals - (coalIntervalIdx + 1)) > mmin)
                    groupIdx += 1;

                coalIntervalIdx += 1;
            }
        }

        return logP;
    }
}
