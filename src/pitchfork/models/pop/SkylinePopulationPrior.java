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

package pitchfork.models.pop;

import beast.core.Distribution;
import beast.core.Input;
import beast.core.State;
import beast.math.distributions.ParametricDistribution;

import java.util.List;
import java.util.Random;

public class SkylinePopulationPrior extends Distribution {

    public Input<SkylinePopulationFunction> skylineInput = new Input<>(
            "skyline",
            "Skyline population function to which to apply prior.",
            Input.Validate.REQUIRED);

    public Input<ParametricDistribution> priorOnN0Input = new Input<>(
            "priorOnN0",
            "Parametric distribution specifying prior on present-day " +
                    "population size.",
            Input.Validate.REQUIRED);

    public Input<ParametricDistribution> priorOnLogNDeltasInput = new Input<>(
            "priorOnLogNDeltas",
            "Parametric distribution specifying prior on differences " +
            "between logarithms of adjacent population sizes.",
            Input.Validate.REQUIRED);

    private SkylinePopulationFunction skyline;
    private ParametricDistribution priorOnN0, priorOnLogNDeltas;

    @Override
    public void initAndValidate() {

        skyline = skylineInput.get();
        priorOnN0 = priorOnN0Input.get();
        priorOnLogNDeltas = priorOnLogNDeltasInput.get();

        super.initAndValidate();
    }

    @Override
    public double calculateLogP() {
        logP = priorOnN0.logDensity(skyline.getIntervalPopSize(0));

        for (int i=0; i<skyline.getSkylineIntervalCount()-1; i++)
            logP += priorOnLogNDeltas.logDensity(skyline.getIntervalDelta(i));

        return logP;
    }

    @Override
    public List<String> getArguments() {
        return null;
    }

    @Override
    public List<String> getConditions() {
        return null;
    }

    @Override
    public void sample(State state, Random random) {

    }
}
