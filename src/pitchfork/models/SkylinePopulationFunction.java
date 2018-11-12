package pitchfork.models;

import beast.core.Input;
import beast.core.Loggable;
import beast.core.parameter.RealParameter;
import beast.evolution.tree.Tree;
import beast.evolution.tree.coalescent.PopulationFunction;
import beast.util.TreeParser;

import java.io.FileNotFoundException;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;

/**
 * Skyline model population function for trees with polytomies.
 */
public class SkylinePopulationFunction extends PopulationFunction.Abstract implements Loggable {

    public Input<Tree> treeInput = new Input<>("tree",
            "Tree on which skyline model is based.", Input.Validate.REQUIRED);

    public Input<RealParameter> relativePopSizesInput = new Input<>("relativePopSizes",
            "Population sizes in intervals", Input.Validate.REQUIRED);

    public Input<RealParameter> popSizeScaleInput = new Input<>("popSizeScale",
            "Population size scale parameter.",
            Input.Validate.REQUIRED);

    public Input<Boolean> piecewiseLinearInput = new Input<>("piecewiseLinear",
            "Use piecewise linear rather than piecewise constant " +
                    "population function.", false);

    public Input<RealParameter> epsilonInput = new Input<>("epsilon",
            "Epsilon parameter.", Input.Validate.REQUIRED);

    public Input<Boolean> epsilonIsRelativeInput = new Input<>("epsilonIsRelative",
            "If true, epsilon parameter is relative to tree height.",
            false);

    public Input<Boolean> initializePopSizesInput = new Input<>("initializePopSizes",
            "If true, automatically initialize relative population size to all 1s.",
            false);


    Tree tree;
    RealParameter relativePopSizes, popSizeScale, epsilon;
    boolean piecewiseLinear, epsilonIsRelative;

    private boolean dirty;

    double[] intervalBoundaryTimes;

    List<Double> activeBoundaryTimes, activeBoundaryIntensities, activePopSizes;

    int nIntervals;

    @Override
    public void initAndValidate() {
        super.initAndValidate();

        tree = treeInput.get();
        relativePopSizes = relativePopSizesInput.get();
        popSizeScale = popSizeScaleInput.get();
        piecewiseLinear = piecewiseLinearInput.get();
        epsilon = epsilonInput.get();
        epsilonIsRelative = epsilonIsRelativeInput.get();

        nIntervals = tree.getInternalNodeCount();

        intervalBoundaryTimes = new double[nIntervals+1];
        activeBoundaryTimes = new ArrayList<>();
        activeBoundaryIntensities = new ArrayList<>();
        activePopSizes = new ArrayList<>();

        if (initializePopSizesInput.get()) {
            relativePopSizes.setDimension(nIntervals+1);
            for (int i=0; i<=nIntervals; i++)
                relativePopSizes.setValue(i, 1.0);
        }

        dirty = true;
    }

    @Override
    public List<String> getParameterIds() {
        List<String> ids = new ArrayList<>();
        ids.add(treeInput.get().getID());
        ids.add(relativePopSizesInput.get().getID());
        ids.add(epsilonInput.get().getID());

        return ids;
    }

    @Override
    public void prepare() {

        if (!dirty)
            return;

        // Set up interval boundary times
        intervalBoundaryTimes[0] = 0.0;
        for (int i = 0; i<tree.getInternalNodeCount(); i++)
            intervalBoundaryTimes[i+1] = tree.getNode(i+tree.getLeafNodeCount()).getHeight();
        Arrays.sort(intervalBoundaryTimes);

        double trueEpsilon = epsilonIsRelative
                ? epsilon.getValue()*tree.getRoot().getHeight()
                : epsilon.getValue();

        // Identify active interval boundaries and pop sizes
        activeBoundaryTimes.clear();
        activeBoundaryTimes.add(0.0);
        activePopSizes.clear();
        activePopSizes.add(relativePopSizes.getValue(0)*popSizeScale.getValue());
        double deltat = 0.0;
        for (int i=0; i<nIntervals; i++) {
            deltat += intervalBoundaryTimes[i+1]-intervalBoundaryTimes[i];

            if (deltat > trueEpsilon) {
                activeBoundaryTimes.add(intervalBoundaryTimes[i+1]);
                activePopSizes.add(relativePopSizes.getValue(i+1)*popSizeScale.getValue());
                deltat = 0.0;
            }
        }

        // Precompute intensities at active boundaries
        activeBoundaryIntensities.clear();
        activeBoundaryIntensities.add(0.0);

        if (!piecewiseLinearInput.get()) {
            for (int i = 1; i < activeBoundaryTimes.size(); i++) {
                activeBoundaryIntensities.add(activeBoundaryIntensities.get(i - 1)
                        + (activeBoundaryTimes.get(i) - activeBoundaryTimes.get(i-1)) / activePopSizes.get(i - 1));
            }
        } else {
            for (int i = 1; i < activeBoundaryTimes.size(); i++) {

                if (!relativePopSizes.getValue(i - 1).equals(relativePopSizes.getValue(i))) {
                    activeBoundaryIntensities.add(activeBoundaryIntensities.get(i - 1)
                            + (activeBoundaryTimes.get(i) - activeBoundaryTimes.get(i - 1))
                            / (activePopSizes.get(i) - activePopSizes.get(i - 1))
                            * Math.log(relativePopSizes.getValue(i) / relativePopSizes.getValue(i - 1)));
                } else {
                    activeBoundaryIntensities.add(activeBoundaryIntensities.get(i - 1)
                            + (activeBoundaryTimes.get(i) - activeBoundaryTimes.get(i - 1)) / activePopSizes.get(i - 1));
                }
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
            return relativePopSizes.getValue(0);

        if (t >= tree.getRoot().getHeight())
            return activePopSizes.get(activePopSizes.size()-1);

        int interval = Collections.binarySearch(activeBoundaryTimes, t);

        if (interval<0)
            interval = -(interval + 1) - 1;  // boundary to the left of time.

        if (!piecewiseLinearInput.get())
            return activePopSizes.get(interval);
        else {
            double N0 = activePopSizes.get(interval);
            double N1 = activePopSizes.get(interval+1);

            double t0 = activeBoundaryTimes.get(interval);
            double t1 = activeBoundaryTimes.get(interval+1);

            return N0 + (t - t0)/(t1 - t0)*(N1 - N0);
        }
    }

    @Override
    public double getIntensity(double t) {
        prepare();

        if (t <= 0 )
            return -t/ relativePopSizes.getValue(0);

        if (t >= tree.getRoot().getHeight())
            return activeBoundaryIntensities.get(activeBoundaryIntensities.size()-1)
                    + (t - activeBoundaryTimes.get(activeBoundaryTimes.size()-1))
                    /activePopSizes.get(activePopSizes.size()-1);

        int interval = Collections.binarySearch(activeBoundaryTimes, t);

        if (interval<0)
            interval = -(interval + 1) - 1; // boundary to the left of time.

        if (!piecewiseLinearInput.get())
            return activeBoundaryIntensities.get(interval) + (t - activeBoundaryTimes.get(interval))/activePopSizes.get(interval);
        else {
            double N0 = activePopSizes.get(interval);
            double N1 = activePopSizes.get(interval+1);

            double t0 = activeBoundaryTimes.get(interval);
            double t1 = activeBoundaryTimes.get(interval+1);

            double N = N0 + (t - t0)/(t1 - t0)*(N1 - N0);

            if (N1 != N0)
                return activeBoundaryIntensities.get(interval) + (t1-t0)/(N1-N0)*Math.log(N/N0);
            else
                return activeBoundaryIntensities.get(interval) + (t - t0)/N0;
        }
    }

    @Override
    public double getInverseIntensity(double x) {
        prepare();

        if (x<=0.0)
            return -x* relativePopSizes.getValue(0);

        if (x >= activeBoundaryIntensities.get(activeBoundaryIntensities.size()-1))
            return activeBoundaryTimes.get(activeBoundaryTimes.size()-1)
                    + (x - activeBoundaryIntensities.get(activeBoundaryIntensities.size()-1))
                    *activePopSizes.get(activePopSizes.size()-1);

        int interval = Collections.binarySearch(activeBoundaryIntensities, x);

        if (interval<0)
            interval = -(interval + 1) - 1; // boundary to the left of x

        if (!piecewiseLinearInput.get())
            return activeBoundaryTimes.get(interval)
                    + (x-activeBoundaryIntensities.get(interval))*activePopSizes.get(interval);
        else {
            double N0 = activePopSizes.get(interval);
            double N1 = activePopSizes.get(interval+1);

            double t0 = activeBoundaryTimes.get(interval);
            double t1 = activeBoundaryTimes.get(interval+1);

            double a = N0 - t0*(N1-N0)/(t1-t0);
            double b = (N1-N0)/(t1-t0);

            if (N1 != N0)
                return (N0*Math.exp((N1-N0)/(t1-t0)*(x-activeBoundaryIntensities.get(interval))) - a)/b;
            else
                return activeBoundaryTimes.get(interval) + N0*(x - activeBoundaryIntensities.get(interval));
        }
    }

    /*
     * Loggable implementation:
     */

    @Override
    public void init(PrintStream out) {
        prepare();

        for (int i = 0; i < nIntervals+1; i++) {
            out.print(getID() + ".t" + i + "\t");
            out.print(getID() + ".N" + i + "\t");
        }
    }

    @Override
    public void log(long nSample, PrintStream out) {
        prepare();

        for (int i = 0; i < nIntervals+1; i++) {
            out.print(intervalBoundaryTimes[i] + "\t");
            out.print(getPopSize(intervalBoundaryTimes[i]) + "\t");
        }
    }

    @Override
    public void close(PrintStream out) { }

    /**
     * Main method for testing.
     *
     * @param args
     */
    public static void main(String[] args) throws FileNotFoundException {

        String newickStr = "((t1:1.0,t2:1.0):1.0,t3:2.0):0.0;";

        TreeParser tree = new TreeParser(newickStr, false, false, true, 0);


        SkylinePopulationFunction skyline = new SkylinePopulationFunction();
        skyline.initByName(
                "tree", tree,
                "relativePopSizes", new RealParameter("1.0 0.5 2.0"),
                "popSizeScale", new RealParameter("1.0"),
                "epsilon", new RealParameter("0.1"),
                "epsilonIsRelative", true,
                "piecewiseLinear", true);

        try (PrintStream ps = new PrintStream("out.txt")){
            ps.println("t N intensity intensityInv");
            double t = 0.0;
            while (t<3) {
                ps.format("%g %g %g %g\n", t,
                        skyline.getPopSize(t), skyline.getIntensity(t),
                        skyline.getInverseIntensity(skyline.getIntensity(t)));
                t += 0.001;
            }
            ps.close();
        }
    }
}
