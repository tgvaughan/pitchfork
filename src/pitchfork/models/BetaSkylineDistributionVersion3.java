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

public class BetaSkylineDistributionVersion3 extends AbstractBetaSkylineDistribution {

    public BetaSkylineDistributionVersion3() {
        super();
    }


    // This version updates N after a certain time
    @Override
    public double calculateLogP() {
        logP = 0.0;

        double t = 0;
        double N = skylinePopulationsInput.get().getValue(0);
        int seenCoalescentEvents = 0;
        double dt_track = 0;                                        // this variable is first dt after every round and then t_update is deducted until dt_track < t_update
        double t_update = 0.2;                                      //this variable defines the length of the time intervals. Can easily be implemented as beast Input
        int counter = 0;                                             // counter for groups

        for (int i = 0; i < collapsedTreeIntervals.getIntervalCount(); i++) {
            // Get interval details
            double dt = collapsedTreeIntervals.getInterval(i);
            dt_track += dt;

            int n = collapsedTreeIntervals.getLineageCount(i);

            while (dt_track >= t_update){
                        logP += -betaCoalescentModel.getTotalCoalRate(n)*t_update/N;  // here we calculate the logP of the waiting time for the pre-defined interval t_update as long dt_track is bigger than t_update
                        dt_track -= t_update;
                        counter += 1;
                        N = skylinePopulationsInput.get().getValue(counter);
            }

            logP += -betaCoalescentModel.getTotalCoalRate(n)*dt_track/N;      // here we calculate the logP for the time interval which is smaller than t_update


            // Increment time
            t += dt;

            if (collapsedTreeIntervals.getIntervalType(i) == IntervalType.COALESCENT) {
                // Beta-coalescent event contribution

                int k = collapsedTreeIntervals.getCoalescentEvents(i)+1;
                logP += betaCoalescentModel.getLogLambda(n, k) + Binomial.logChoose(n, k) - Math.log(N);
            }
        }
        return logP;
    }

}