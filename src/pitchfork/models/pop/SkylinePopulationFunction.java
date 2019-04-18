package pitchfork.models.pop;

import beast.core.Input;
import beast.core.Loggable;
import beast.core.parameter.IntegerParameter;
import beast.core.parameter.RealParameter;
import beast.evolution.tree.Node;
import beast.evolution.tree.Tree;
import beast.evolution.tree.coalescent.PopulationFunction;
import beast.util.TreeParser;
import pitchfork.Pitchforks;

import java.io.FileNotFoundException;
import java.io.PrintStream;
import java.util.*;
import java.util.stream.Collectors;

/**
 * Skyline model population function for trees with polytomies.
 */
public class SkylinePopulationFunction extends PopulationFunction.Abstract implements Loggable {

    public Input<Tree> treeInput = new Input<>("tree",
            "Tree on which skyline model is based.", Input.Validate.REQUIRED);

    public Input<RealParameter> popSizesInput = new Input<>("relPopSizes",
            "Relative population sizes.", Input.Validate.REQUIRED);

    public Input<IntegerParameter> skylineIntervalCountInput = new Input<>("skylineIntervalCount",
            "Total number of skyline intervals.", Input.Validate.REQUIRED);

    public Input<Boolean> piecewiseLinearInput = new Input<>("piecewiseLinear",
            "Use piecewise linear rather than piecewise constant " +
                    "population function.", false);

    private Tree tree;
    private RealParameter popSizes;
    private IntegerParameter skylineIntervalCount;
    private boolean piecewiseLinear;

    private boolean dirty;

    private List<Double> intervalStartTimes, intervalStartIntensities;

    @Override
    public void initAndValidate() {
        tree = treeInput.get();
        popSizes = popSizesInput.get();
        skylineIntervalCount = skylineIntervalCountInput.get();
        piecewiseLinear = piecewiseLinearInput.get();

        int maxCoalCount = tree.getInternalNodeCount();

        popSizes.setDimension(maxCoalCount + 1);
        for (int i = 0; i<= maxCoalCount; i++)
            popSizes.setValue(i, 1.0);

        dirty = true;

        super.initAndValidate();
    }

    @Override
    public List<String> getParameterIds() {
        List<String> ids = new ArrayList<>();
        ids.add(treeInput.get().getID());
        ids.add(popSizesInput.get().getID());
        ids.add(skylineIntervalCountInput.get().getID());

        return ids;
    }

    /**
     * Compute and return the current true number of skyline intervals.
     * (This differs from the desired number, which is specified by
     * the skylineIntervalCount input, as the true number is constrained
     * by the actual number of coalescent intervals in the current tree.
     *
     * @return current true number  of skyline intervals.
     */
    public int getTrueSkylineIntervalCount() {
        int nCoal = Pitchforks.getTrueInternalNodes(tree).size();
        int coalsPerSkyline = (int)Math.ceil(nCoal/(double)skylineIntervalCount.getValue());

        return nCoal/coalsPerSkyline;
    }

    @Override
    public void prepare() {

        if (!dirty)
            return;

        List<Node> internalNodes = Pitchforks.getTrueInternalNodes(tree);
        internalNodes.sort(Comparator.comparingDouble(Node::getHeight));

        int nCoal = internalNodes.size();

        // Always at least one coalescent event per skyline interval:
        int coalsPerSkyline = (int)Math.ceil(nCoal/(double)skylineIntervalCount.getValue());

        // Set up interval start times
        intervalStartTimes = new ArrayList<>();
        intervalStartTimes.add(0.0);
        for (int i=0; i<nCoal; i++) {
            if (i % coalsPerSkyline == 0)
                intervalStartTimes.add(internalNodes.get(i).getHeight());
        }

        // Precompute intensities at interval start times
        intervalStartIntensities = new ArrayList<>();
        intervalStartIntensities.add(0.0);

        if (!piecewiseLinear) {
            for (int i = 1; i < intervalStartTimes.size(); i++) {
                intervalStartIntensities.add(intervalStartIntensities.get(i - 1)
                        + (intervalStartTimes.get(i) - intervalStartTimes.get(i-1)) / popSizes.getValue(i - 1));
            }
        } else {
            for (int i = 1; i < intervalStartTimes.size(); i++) {
                double N0 = popSizes.getValue(i-1);
                double N1 = popSizes.getValue(i);

                double t0 = intervalStartTimes.get(i-1);
                double t1 = intervalStartTimes.get(i);

                if (N0 != N1) {
                    intervalStartIntensities.add(intervalStartIntensities.get(i - 1)
                            + (t1 - t0) / (N1-N0) * Math.log(N1 / N0));
                } else {
                    intervalStartIntensities.add(intervalStartIntensities.get(i - 1)
                            + (t1 - t0) / N0);
                }
            }
        }

        dirty = false;
    }

    @Override
    public double getPopSize(double t) {
        prepare();

        if (t <= 0)
            return popSizes.getValue(0);

        if (t >= tree.getRoot().getHeight())
            return popSizes.getValue(intervalStartTimes.size()-1);

        int interval = Collections.binarySearch(intervalStartTimes, t);

        if (interval<0)
            interval = -(interval + 1) - 1;  // boundary to the left of time.

        if (!piecewiseLinear)
            return popSizes.getValue(interval);
        else {
            double N0 = popSizes.getValue(interval);
            double N1 = popSizes.getValue(interval+1);

            double t0 = intervalStartTimes.get(interval);
            double t1 = intervalStartTimes.get(interval+1);

            return N0 + (t - t0)/(t1 - t0)*(N1 - N0);
        }
    }

    @Override
    public double getIntensity(double t) {
        prepare();

        if (t <= 0 )
            return -t/ popSizes.getValue(0);

        if (t >= tree.getRoot().getHeight())
            return intervalStartIntensities.get(intervalStartIntensities.size()-1)
                    + (t - intervalStartTimes.get(intervalStartTimes.size()-1))
                    /popSizes.getValue(intervalStartTimes.size()-1);

        int interval = Collections.binarySearch(intervalStartTimes, t);

        if (interval<0)
            interval = -(interval + 1) - 1; // boundary to the left of time.

        if (!piecewiseLinear)
            return intervalStartIntensities.get(interval) + (t - intervalStartTimes.get(interval))/popSizes.getValue(interval);
        else {
            double N0 = popSizes.getValue(interval);
            double N1 = popSizes.getValue(interval+1);

            double t0 = intervalStartTimes.get(interval);
            double t1 = intervalStartTimes.get(interval+1);

            double N = N0 + (t - t0)/(t1 - t0)*(N1 - N0);

            if (N1 != N0)
                return intervalStartIntensities.get(interval) + (t1-t0)/(N1-N0)*Math.log(N/N0);
            else
                return intervalStartIntensities.get(interval) + (t - t0)/N0;
        }
    }

    @Override
    public double getInverseIntensity(double x) {
        prepare();

        if (x<=0.0)
            return x * popSizes.getValue(0);

        if (x >= intervalStartIntensities.get(intervalStartIntensities.size()-1))
            return intervalStartTimes.get(intervalStartTimes.size()-1)
                    + (x - intervalStartIntensities.get(intervalStartIntensities.size()-1))
                    * popSizes.getValue(intervalStartTimes.size()-1);

        int interval = Collections.binarySearch(intervalStartIntensities, x);

        if (interval<0)
            interval = -(interval + 1) - 1; // boundary to the left of x

        if (!piecewiseLinear)
            return intervalStartTimes.get(interval)
                    + (x-intervalStartIntensities.get(interval))*popSizes.getValue(interval);
        else {
            double N0 = popSizes.getValue(interval);
            double N1 = popSizes.getValue(interval+1);

            double t0 = intervalStartTimes.get(interval);
            double t1 = intervalStartTimes.get(interval+1);

            if (N1 != N0)
                return Math.exp((x-intervalStartIntensities.get(interval))*(N1-N0)/(t1-t0) - 1)*(t1-t0)*N0/(N1-N0);
            else
                return intervalStartTimes.get(interval) + N0*(x - intervalStartIntensities.get(interval));
        }
    }


    /*
     * StateNode implementation:
     */

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


    /*
     * Loggable implementation:
     */

    @Override
    public void init(PrintStream out) {
        for (int i = 0; i < tree.getInternalNodeCount()+1; i++) {
            out.print(getID() + ".t" + i + "\t");
            out.print(getID() + ".N" + i + "\t");
        }
    }

    @Override
    public void log(long nSample, PrintStream out) {
        double[] internalNodeTimes = tree.getInternalNodes().stream()
                .mapToDouble(Node::getHeight)
                .sorted()
                .toArray();

        out.print(0.0 + "\t" + getPopSize(0.0) + "\t");

        for (double internalNodeTime : internalNodeTimes) {
            out.print(internalNodeTime + "\t");
            out.print(getPopSize(internalNodeTime) + "\t");
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

        String newickStr = "((t1:1.0,t2:1.0,t3:1.0):1.0,t4:2.0):0.0;";

        TreeParser tree = new TreeParser(newickStr, false, false, true, 0);


        SkylinePopulationFunction skyline = new SkylinePopulationFunction();
        skyline.initByName(
                "tree", tree,
                "relativePopSizes", new RealParameter("1.0 0.5 5.0 2.0"),
                "popSizeScale", new RealParameter("1.0"),
                "epsilon", new RealParameter("0.1"),
                "epsilonIsRelative", true,
                "piecewiseLinear", true);

        try (PrintStream ps = new PrintStream("out.txt")){
            ps.println("t N intensity intensityInv");
            double t = -1.0;
            while (t<3) {
                ps.format("%g %g %g %g\n", t,
                        skyline.getPopSize(t), skyline.getIntensity(t),
                        skyline.getInverseIntensity(skyline.getIntensity(t)));
                t += 0.001;
            }
        }
    }
}
