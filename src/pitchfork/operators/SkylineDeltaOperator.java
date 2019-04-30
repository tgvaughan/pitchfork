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
import beast.core.Operator;
import beast.core.parameter.RealParameter;
import beast.util.Randomizer;
import pitchfork.models.pop.SkylinePopulationFunction;

public class SkylineDeltaOperator extends Operator {

    public Input<SkylinePopulationFunction> skylineInput = new Input<>(
            "skyline",
            "Skyline population function",
            Input.Validate.REQUIRED);

    public Input<Double> windowSizeInput = new Input<>(
            "windowSize",
            "Size of uniform proposal distribution for new element values.",
            1.0);

    private SkylinePopulationFunction skyline;
    private RealParameter deltaLogPopSizes;
    private double windowSize;

    @Override
    public void initAndValidate() {
        skyline = skylineInput.get();
        deltaLogPopSizes = skyline.logNDeltasInput.get();
        windowSize = windowSizeInput.get();
    }

    @Override
    public double proposal() {

        if (skyline.getSkylineIntervalCount()==1)
            return Double.NEGATIVE_INFINITY;

        int idx = Randomizer.nextInt(skyline.getSkylineIntervalCount()-1);

        double newValue = deltaLogPopSizes.getValue(idx) + windowSize*(Randomizer.nextDouble()-0.5);

        if (newValue > deltaLogPopSizes.getUpper() || newValue < deltaLogPopSizes.getLower())
            return Double.NEGATIVE_INFINITY;

        deltaLogPopSizes.setValue(idx, newValue);

        return 0;
    }
}
