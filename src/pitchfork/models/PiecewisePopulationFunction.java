package pitchfork.models;

import beast.core.Description;
import beast.core.Input;
import beast.core.Loggable;
import beast.core.parameter.RealParameter;
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
            "Change times parameter.", Input.Validate.REQUIRED);

    RealParameter popSizes, changeTimes;

    double[] intensities;
    double[] groupBoundaries;

    @Override
    public void initAndValidate() {
        popSizes = popSizesInput.get();
        changeTimes = changeTimesInput.get();

        groupBoundaries = new double[popSizes.getDimension()];
        intensities = new double[popSizes.getDimension()];

        super.initAndValidate();
    }

    @Override
    public List<String> getParameterIds() {
        List<String> ids = new ArrayList<>();
        ids.add(popSizesInput.get().getID());
        ids.add(changeTimesInput.get().getID());

        return ids;
    }

    @Override
    public void prepare() {

        groupBoundaries[0] = 0.0;
        for (int i=1; i<groupBoundaries.length; i++)
            groupBoundaries[i] = changeTimes.getValue(i-1);

        intensities[0] = 0.0;

        for (int i = 1; i < intensities.length; i++) {
            intensities[i] = intensities[i - 1]
                    + (groupBoundaries[i] - groupBoundaries[i-1]) / popSizes.getValue(i - 1);
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

        if (t >= groupBoundaries[groupBoundaries.length-1])
            return popSizes.getValue(popSizes.getDimension()-1);

        int interval = Arrays.binarySearch(groupBoundaries, t);

        if (interval<0)
            interval = -(interval + 1) - 1;  // boundary to the left of time.

        return popSizes.getValue(interval);
    }

    @Override
    public double getIntensity(double t) {
        prepare();

        if (t <= 0 )
            return -t/popSizes.getValue(0);

        if (t >= groupBoundaries[groupBoundaries.length-1])
            return intensities[intensities.length-1]
                    + (t-groupBoundaries[intensities.length-1])
                    /popSizes.getValue(popSizes.getDimension()-1);

        int interval = Arrays.binarySearch(groupBoundaries, t);

        if (interval<0)
            interval = -(interval + 1) - 1; // boundary to the left of time.

        return intensities[interval] + (t-groupBoundaries[interval])/popSizes.getValue(interval);
    }

    @Override
    public double getInverseIntensity(double x) {
        prepare();

        if (x<=0.0)
            return -x*popSizes.getValue(0);

        if (x >= intensities[intensities.length-1])
            return groupBoundaries[groupBoundaries.length-1]
                    + (x - intensities[intensities.length-1])
                    *popSizes.getValue(popSizes.getDimension()-1);

        int interval = Arrays.binarySearch(intensities, x);

        if (interval<0)
            interval = -(interval + 1) - 1; // boundary to the left of x

        return groupBoundaries[interval]
                + (x-intensities[interval])*popSizes.getValue(interval);
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
            out.print(groupBoundaries[i] + "\t");
            out.print(popSizes.getValue(i) + "\t");
        }
    }

    @Override
    public void close(PrintStream out) {

    }
}
