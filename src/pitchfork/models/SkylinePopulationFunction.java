package pitchfork.models;

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

/**
 * Skyline model population function fork Pitchfork.
 */
@Description("Piecewise constant/linear population function.")
public class SkylinePopulationFunction extends PopulationFunction.Abstract implements Loggable {

    public Input<Tree> treeInput = new Input<>("tree",
            "Tree on which skyline model is based.", Input.Validate.REQUIRED);

    public Input<RealParameter> popSizesInput = new Input<>("popSizes",
            "Population sizes in intervals", Input.Validate.REQUIRED);

    public Input<Boolean> piecewiseLinearInput = new Input<>("piecewiseLinear",
            "Use piecewise linear rather than piecewise constant " +
                    "population function.", false);

    public Input<RealParameter> epsilonInput = new Input<>("epsilon",
            "Epsilon parameter.", Input.Validate.REQUIRED);


    Tree tree;
    RealParameter popSizes, epsilon;
    boolean piecewiseLinear;

    private boolean dirty;

    double[] intensities;
    double[] intervalBoundaryTimes;
    boolean[] boundaryActive;

    int nIntervals;

    @Override
    public void initAndValidate() {
        super.initAndValidate();

        tree = treeInput.get();
        popSizes = popSizesInput.get();
        epsilon = epsilonInput.get();

        nIntervals = tree.getInternalNodeCount();

        intervalBoundaryTimes = new double[nIntervals+1];
        intensities = new double[nIntervals+1];
        boundaryActive = new boolean[nIntervals+1];

        dirty = true;
    }

    @Override
    public List<String> getParameterIds() {
        List<String> ids = new ArrayList<>();
        ids.add(treeInput.get().getID());
        ids.add(popSizesInput.get().getID());
        ids.add(epsilonInput.get().getID());

        return ids;
    }

    @Override
    public void prepare() {

        if (!dirty)
            return;

        intervalBoundaryTimes[0] = 0.0;

        for (int i = 0; i<tree.getInternalNodeCount(); i++)
            intervalBoundaryTimes[i+1] = tree.getNode(i+tree.getLeafNodeCount()).getHeight();
        Arrays.sort(intervalBoundaryTimes);

        intensities[0] = 0.0;

        if (!piecewiseLinearInput.get()) {
            for (int i = 1; i < intensities.length; i++) {
                intensities[i] = intensities[i - 1]
                        + (intervalBoundaryTimes[i] - intervalBoundaryTimes[i-1]) / popSizes.getValue(i - 1);
            }
        } else {
            for (int i = 1; i < intensities.length; i++) {

                if (!popSizes.getValue(i - 1).equals(popSizes.getValue(i)))
                    intensities[i] = intensities[i - 1]
                            + (intervalBoundaryTimes[i] - intervalBoundaryTimes[i - 1])
                            / (popSizes.getValue(i) - popSizes.getValue(i - 1))
                            * Math.log(popSizes.getValue(i) / popSizes.getValue(i - 1));
                else
                    intensities[i] = intensities[i-1]
                            + (intervalBoundaryTimes[i] - intervalBoundaryTimes[i-1])/popSizes.getValue(i-1);
            }
        }

        dirty = false;
    }

    @Override
    protected boolean requiresRecalculation() {
        dirty = true;
        return true;
    }

    @Override
    protected void store() {
        dirty = true;
        super.store();
    }

    @Override
    protected void restore() {
        dirty = true;
        super.restore();
    }



    @Override
    public double getPopSize(double t) {
        prepare();

        if (t <= 0)
            return popSizes.getValue(0);

        if (t >= intervalBoundaryTimes[intervalBoundaryTimes.length-1])
            return popSizes.getValue(popSizes.getDimension()-1);

        int interval = Arrays.binarySearch(intervalBoundaryTimes, t);

        if (interval<0)
            interval = -(interval + 1) - 1;  // boundary to the left of time.

        if (!piecewiseLinearInput.get())
            return popSizes.getValue(interval);
        else {
            double N0 = popSizes.getValue(interval);
            double N1 = popSizes.getValue(interval+1);

            double t0 = intervalBoundaryTimes[interval];
            double t1 = intervalBoundaryTimes[interval+1];

            return N0 + (t - t0)/(t1 - t0)*(N1 - N0);
        }
    }

    @Override
    public double getIntensity(double t) {
        prepare();

        if (t <= 0 )
            return -t/popSizes.getValue(0);

        if (t >= intervalBoundaryTimes[intervalBoundaryTimes.length-1])
            return intensities[intensities.length-1]
                    + (t- intervalBoundaryTimes[intensities.length-1])
                    /popSizes.getValue(popSizes.getDimension()-1);

        int interval = Arrays.binarySearch(intervalBoundaryTimes, t);

        if (interval<0)
            interval = -(interval + 1) - 1; // boundary to the left of time.

        if (!piecewiseLinearInput.get())
            return intensities[interval] + (t- intervalBoundaryTimes[interval])/popSizes.getValue(interval);
        else {
            double N0 = popSizes.getValue(interval);
            double N1 = popSizes.getValue(interval+1);

            double t0 = intervalBoundaryTimes[interval];
            double t1 = intervalBoundaryTimes[interval+1];

            double N = N0 + (t - t0)/(t1 - t0)*(N1 - N0);

            if (N1 != N0)
                return intensities[interval] + (t1-t0)/(N1-N0)*Math.log(N/N0);
            else
                return intensities[interval] + (t - t0)/N0;
        }
    }

    @Override
    public double getInverseIntensity(double x) {
        prepare();

        if (x<=0.0)
            return -x*popSizes.getValue(0);

        if (x >= intensities[intensities.length-1])
            return intervalBoundaryTimes[intervalBoundaryTimes.length-1]
                    + (x - intensities[intensities.length-1])
                    *popSizes.getValue(popSizes.getDimension()-1);

        int interval = Arrays.binarySearch(intensities, x);

        if (interval<0)
            interval = -(interval + 1) - 1; // boundary to the left of x

        if (!piecewiseLinearInput.get())
            return intervalBoundaryTimes[interval]
                    + (x-intensities[interval])*popSizes.getValue(interval);
        else {
            double N0 = popSizes.getValue(interval);
            double N1 = popSizes.getValue(interval+1);

            double t0 = intervalBoundaryTimes[interval];
            double t1 = intervalBoundaryTimes[interval+1];

            double a = N0 - t0*(N1-N0)/(t1-t0);
            double b = (N1-N0)/(t1-t0);

            if (N1 != N0)
                return (N0*Math.exp((N1-N0)/(t1-t0)*(x-intensities[interval])) - a)/b;
            else
                return intervalBoundaryTimes[interval] + N0*(x - intensities[interval]);
        }
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
            out.print(intervalBoundaryTimes[i] + "\t");
            out.print(popSizes.getValue(i) + "\t");
        }
    }

    @Override
    public void close(PrintStream out) {

    }
}
