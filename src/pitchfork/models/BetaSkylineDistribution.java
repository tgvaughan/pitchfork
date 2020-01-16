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

import beast.core.Input;
import beast.core.parameter.RealParameter;
import beast.evolution.tree.TreeDistribution;

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

    @Override
    public void initAndValidate() {
    }
}
