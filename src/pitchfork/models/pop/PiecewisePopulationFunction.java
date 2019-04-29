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

import beast.core.Description;
import beast.core.Input;
import beast.core.Loggable;
import beast.core.parameter.RealParameter;
import beast.evolution.tree.Tree;
import beast.evolution.tree.coalescent.PopulationFunction;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;


@Description("Piecewise constant/linear population function.")
public class PiecewisePopulationFunction extends PopulationFunction.Abstract implements Loggable {

    public Input<RealParameter> popSizesInput = new Input<>("popSizes",
            "Population sizes in intervals", Input.Validate.REQUIRED);

    public Input<RealParameter> changeTimesInput = new Input<>("changeTimes",
            "Change times parameter.");

    public Input<Tree> treeInput = new Input<>("tree",
            "Tree used to select evenly-spaced interval times.",
            Input.Validate.XOR, changeTimesInput);

    RealParameter popSizes, changeTimes;
    int intervalCount;
    boolean evenlySpaced;

    Tree tree;

    private double[] intervalStartTimes, intervalStartIntensities;

    @Override
    public void initAndValidate() {
        popSizes = popSizesInput.get();
        changeTimes = changeTimesInput.get();
        tree = treeInput.get();
        evenlySpaced = (changeTimes == null);

        if (evenlySpaced)
            intervalCount = popSizes.getDimension();
        else
            intervalCount = changeTimes.getDimension();

        intervalStartTimes = new double[intervalCount];
        intervalStartIntensities = new double[intervalCount];

        super.initAndValidate();
    }

    @Override
    public List<String> getParameterIds() {
        List<String> ids = new ArrayList<>();
        ids.add(popSizesInput.get().getID());
        if (changeTimesInput.get() != null)
            ids.add(changeTimesInput.get().getID());

        return ids;
    }

    @Override
    public void prepare() {

        intervalStartTimes[0] = 0.0;

        if (evenlySpaced) {
            double intervalLength = tree.getRoot().getHeight()/intervalCount;

            for (int i=1; i<intervalStartTimes.length; i++) {
                intervalStartTimes[i] = i*intervalLength;
            }

        } else {
            for (int i = 1; i < intervalStartTimes.length; i++)
                intervalStartTimes[i] = changeTimes.getValue(i - 1);
        }

        intervalStartIntensities[0] = 0.0;

        for (int i = 1; i < intervalStartIntensities.length; i++) {
            intervalStartIntensities[i] = intervalStartIntensities[i - 1]
                    + (intervalStartTimes[i] - intervalStartTimes[i-1]) / popSizes.getValue(i - 1);
        }
    }

    @Override
    protected boolean requiresRecalculation() {
        return true;
    }

    @Override
    protected void store() {
        super.store();
    }

    @Override
    protected void restore() {
        super.restore();
    }

    @Override
    public double getPopSize(double t) {
        prepare();

        if (t <= 0)
            return popSizes.getValue(0);

        if (t >= intervalStartTimes[intervalStartTimes.length-1])
            return popSizes.getValue(popSizes.getDimension()-1);

        int interval = Arrays.binarySearch(intervalStartTimes, t);

        if (interval<0)
            interval = -(interval + 1) - 1;  // boundary to the left of time.

        return popSizes.getValue(interval);
    }

    @Override
    public double getIntensity(double t) {
        prepare();

        if (t <= 0 )
            return -t/popSizes.getValue(0);

        if (t >= intervalStartTimes[intervalStartTimes.length-1])
            return intervalStartIntensities[intervalStartIntensities.length-1]
                    + (t- intervalStartTimes[intervalStartIntensities.length-1])
                    /popSizes.getValue(popSizes.getDimension()-1);

        int interval = Arrays.binarySearch(intervalStartTimes, t);

        if (interval<0)
            interval = -(interval + 1) - 1; // boundary to the left of time.

        return intervalStartIntensities[interval] + (t- intervalStartTimes[interval])/popSizes.getValue(interval);
    }

    @Override
    public double getInverseIntensity(double x) {
        prepare();

        if (x<=0.0)
            return -x*popSizes.getValue(0);

        if (x >= intervalStartIntensities[intervalStartIntensities.length-1])
            return intervalStartTimes[intervalStartTimes.length-1]
                    + (x - intervalStartIntensities[intervalStartIntensities.length-1])
                    *popSizes.getValue(popSizes.getDimension()-1);

        int interval = Arrays.binarySearch(intervalStartIntensities, x);

        if (interval<0)
            interval = -(interval + 1) - 1; // boundary to the left of x

        return intervalStartTimes[interval]
                + (x- intervalStartIntensities[interval])*popSizes.getValue(interval);
    }

    // Loggable implementation:

    @Override
    public void init(PrintStream out) {
        prepare();

        for (int i=0; i<popSizes.getDimension(); i++) {
            out.print(getID() + ".t" + i + "\t");
            out.print(getID() + ".N" + i + "\t");
        }
    }

    @Override
    public void log(long nSample, PrintStream out) {
        prepare();

        for (int i=0; i<popSizes.getDimension(); i++) {
            out.print(intervalStartTimes[i] + "\t");
            out.print(popSizes.getValue(i) + "\t");
        }
    }

    @Override
    public void close(PrintStream out) {

    }
}
