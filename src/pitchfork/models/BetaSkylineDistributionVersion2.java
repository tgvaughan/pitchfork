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
import beast.base.util.Binomial;

public class BetaSkylineDistributionVersion2 extends AbstractBetaSkylineDistribution {

    public BetaSkylineDistributionVersion2() {
        super();
    }


    // This version ignores the number of binary coalescent events and updates N for polytomy as one coalescent event
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

            // Waiting time contribution   correct!
            logP += -betaCoalescentModel.getTotalCoalRate(n)*dt/N;


            // Increment time
            t += dt;

            if (collapsedTreeIntervals.getIntervalType(i) == IntervalType.COALESCENT) {
                // Beta-coalescent event contribution

                int k = collapsedTreeIntervals.getCoalescentEvents(i)+1;
                logP += betaCoalescentModel.getLogLambda(n, k) + Binomial.logChoose(n, k) - Math.log(N);

                // while loop to update N until the seen coalescent events are smaller or equal than the actual group size
                seenCoalescentEvents += 1;                    // independent of polytomy
                counter += collapsedTreeIntervals.getCoalescentEvents(i);

                if (counter < totalCoalecents) {           //because we don't want to update N after the last coalescent event
                    if (seenCoalescentEvents >= groupSizesInput.get().getValue(group)) {    //while loop is not necessary because we add seenCoalescentEvents + 1 even for polytomy and we can go at most one group further
                        seenCoalescentEvents -= groupSizesInput.get().getValue(group);
                        group += 1;
                    }
                    N = skylinePopulationsInput.get().getValue(group);
                }
            }
        }
        return logP;
    }
}