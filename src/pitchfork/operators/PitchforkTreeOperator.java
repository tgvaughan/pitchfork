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

package pitchfork.operators;

import beast.core.Input;
import beast.core.parameter.RealParameter;
import beast.evolution.operators.TreeOperator;
import beast.util.Randomizer;
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

        int initialIntervalCount = skyline.getSkylineIntervalCount();

        double logHR = pitchforkProposal();

        if (isSkylineSafe() || logHR == Double.NEGATIVE_INFINITY)
            return logHR;

        // Skyline function update.

        int finalIntervalCount = skyline.getSkylineIntervalCount();

        RealParameter deltaLogPopSizes = skyline.deltaLogPopSizesInput.get();

        // TODO Fix this!!
//        if (finalIntervalCount > initialIntervalCount) {
//
//            for (int i=initialIntervalCount; i<finalIntervalCount; i++) {
//                deltaLogPopSizes.setValue(i,
//                        Randomizer.nextExponential(1.0/deltaLogPopSizes.getValue(0)));
//
//                logHR -= -deltaLogPopSizes.getValue(i)/deltaLogPopSizes.getValue(0)
//                        + Math.log(1.0/deltaLogPopSizes.getValue(0));
//            }
//
//        } else if (finalIntervalCount < initialIntervalCount) {
//            for (int i=finalIntervalCount; i<initialIntervalCount; i++) {
//                logHR += -deltaLogPopSizes.getValue(i)/deltaLogPopSizes.getValue(0)
//                        + Math.log(1.0/deltaLogPopSizes.getValue(0));
//            }
//        }

        return logHR;
    }

    protected abstract double pitchforkProposal();

    abstract boolean isSkylineSafe();
}
