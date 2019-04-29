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

/**
 * Skyline model population function for trees with polytomies.
 */
public class SkylinePopulationFunction extends PopulationFunction.Abstract implements Loggable {

    public Input<Tree> treeInput = new Input<>(
            "tree",
            "Tree on which skyline model is based.",
            Input.Validate.REQUIRED);

    public Input<RealParameter> N0Input = new Input<>(
            "N0",
            "Present-day population sizes.",
            Input.Validate.REQUIRED);

    public Input<RealParameter> deltaLogPopSizesInput = new Input<>(
            "deltaLogPopSizes",
            "Differences between logarithms of adjacent population sizes.",
            Input.Validate.REQUIRED);

    public Input<IntegerParameter> maxSkylineIntervalCountInput = new Input<>(
            "maxSkylineIntervalCount",
            "Total number of skyline intervals.",
            Input.Validate.REQUIRED);

    private Tree tree;
    private RealParameter N0, deltaLogPopSizes;
    private IntegerParameter maxSkylineIntervalCount;

    private boolean dirty;

    private List<Double> intervalStartTimes, intervalStartIntensities, intervalPopSizes;
    private List<Integer> skylineCountList;

    @Override
    public void initAndValidate() {
        tree = treeInput.get();
        N0 = N0Input.get();
        deltaLogPopSizes = deltaLogPopSizesInput.get();
        maxSkylineIntervalCount = maxSkylineIntervalCountInput.get();

        int maxCoalCount = tree.getInternalNodeCount();

        deltaLogPopSizes.setDimension(maxCoalCount-1);

        dirty = true;

        super.initAndValidate();
    }

    @Override
    public List<String> getParameterIds() {
        List<String> ids = new ArrayList<>();
        ids.add(treeInput.get().getID());
        ids.add(N0Input.get().getID());
        ids.add(deltaLogPopSizesInput.get().getID());
        ids.add(maxSkylineIntervalCountInput.get().getID());

        return ids;
    }


    @Override
    public void prepare() {

        if (!dirty)
            return;

        List<Node> internalNodes = Pitchforks.getTrueInternalNodes(tree);
        internalNodes.sort(Comparator.comparingDouble(Node::getHeight));

        int nCoal = internalNodes.size();

        // Always at least one coalescent event per skyline interval:
        int coalsPerSkyline = (int)Math.ceil(nCoal/(double) maxSkylineIntervalCount.getValue());

        // Set up interval start times
        intervalStartTimes = new ArrayList<>();
        intervalStartTimes.add(0.0);
        for (int i=0; i<nCoal-1; i++) {
            if ((i+1) % coalsPerSkyline == 0)
                intervalStartTimes.add(internalNodes.get(i).getHeight());
        }

        int intervalCount = intervalStartTimes.size();

        // Precompute population sizes:
        intervalPopSizes = new ArrayList<>();
        intervalPopSizes.add(N0.getValue());
        double prevLogPopSize = Math.log(N0.getValue());
        for (int i=1; i<intervalCount; i++) {
            double thisLogPopSize = prevLogPopSize + deltaLogPopSizes.getValue(i-1);

            intervalPopSizes.add(Math.exp(thisLogPopSize));

            prevLogPopSize = thisLogPopSize;
        }

        // Precompute intensities at interval start times
        intervalStartIntensities = new ArrayList<>();
        intervalStartIntensities.add(0.0);

        for (int i = 1; i < intervalStartTimes.size(); i++) {
            intervalStartIntensities.add(intervalStartIntensities.get(i - 1)
                    + (intervalStartTimes.get(i) - intervalStartTimes.get(i-1)) / intervalPopSizes.get(i-1));
        }

        dirty = false;
    }

    /**
     * Retrieve the actual true number of skyline intervals in the
     * skyline model for the current tree. This differs from the
     * maxSkylineIntervalCount input, as the true number is constrained
     * by the number of coalescent intervals in the current tree, which
     * may vary due to polytomies.
     *
     * @return current true number of skyline intervals.
     */
    public int getSkylineIntervalCount() {
        prepare();

        return intervalStartTimes.size();
    }

    /**
     * Retrieve the population size corresponding to a given skyline interval.
     *
     * @param interval index of interval
     * @return current population size
     */
    public double getIntervalPopSize(int interval) {
        prepare();

        return intervalPopSizes.get(interval);
    }

    @Override
    public double getPopSize(double t) {
        prepare();

        if (t <= 0)
            return intervalPopSizes.get(0);

        if (t >= intervalStartTimes.get(intervalStartTimes.size()-1))
            return intervalPopSizes.get(intervalStartTimes.size()-1);

        int interval = Collections.binarySearch(intervalStartTimes, t);

        if (interval<0)
            interval = -(interval + 1) - 1;  // boundary to the left of time.

        return intervalPopSizes.get(interval);
    }

    @Override
    public double getIntensity(double t) {
        prepare();

        if (t <= intervalStartTimes.get(0))
            return (t-intervalStartTimes.get(0))/intervalPopSizes.get(0);

        if (t >= intervalStartTimes.get(intervalStartTimes.size()-1))
            return intervalStartIntensities.get(intervalStartIntensities.size()-1)
                    + (t - intervalStartTimes.get(intervalStartTimes.size()-1))
                    /intervalPopSizes.get(intervalStartTimes.size()-1);

        int interval = Collections.binarySearch(intervalStartTimes, t);

        if (interval<0)
            interval = -(interval + 1) - 1; // boundary to the left of time.

        return intervalStartIntensities.get(interval)
                + (t - intervalStartTimes.get(interval))
                /intervalPopSizes.get(interval);
    }

    @Override
    public double getInverseIntensity(double x) {
        prepare();

        if (x<=0.0)
            return x * intervalPopSizes.get(0);

        if (x >= intervalStartIntensities.get(intervalStartIntensities.size()-1))
            return intervalStartTimes.get(intervalStartTimes.size()-1)
                    + (x - intervalStartIntensities.get(intervalStartIntensities.size()-1))
                    * intervalPopSizes.get(intervalStartTimes.size()-1);

        int interval = Collections.binarySearch(intervalStartIntensities, x);

        if (interval<0)
            interval = -(interval + 1) - 1; // boundary to the left of x

        return intervalStartTimes.get(interval)
                + (x-intervalStartIntensities.get(interval))*intervalPopSizes.get(interval);
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

//        String newickStr = "((t1:1.0,t2:1.0,t3:1.0):1.0,t4:2.0):0.0;";
        String newickStr =
                "((1:0.17778950960078332,2:0.020499889752980893):0.942815088997381,(((((3:3.404349842587745,((((4:11.459500475823369,(5:0.9956462693336272,6:0.1015638355924926):0.21331396292499072):0.2601671924876241,((7:0.38518837312929044,8:2.594347394361007):0.5877112703291463,9:0.08091626659444806):2.0214263671134773):0.5705748530674519,(10:0.7945143039368698,11:9.092528886522512):0.1550077634275313):0.2654286942261468,12:5.593871735238542):0.7470866910164045):0.2320344962945624,(((13:1.8707417094108614,(14:2.8626951643328447,15:0.3179919681682062):1.3332860012736418):1.5241547190825635,(16:1.2275279595874116,(17:1.8771216890283764,18:20.978191289788278):6.032716448889871):0.39898067867950626):1.0267316806819702,19:4.099463379464418):0.33633044560204794):0.15633271301423246,(((20:3.5878985089939075,21:1.389403088889674):3.8891761424526266,(22:0.1681547577930198,(23:1.6216812894426358,24:1.2143640464890897):1.9915257185493918):0.45363466893641124):0.10093856621934538,(25:16.538882045187492,((26:3.3005303040398246,(27:8.300198065812248,28:19.252294836710753):0.2734490423717091):3.9608845537528747,29:1.2402467058373885):0.3269923070167078):0.20672211809075858):0.39482167726375206):0.5424960630268938,(30:0.07552685557481986,(31:0.4992190270298451,(((32:1.5064082096688765,33:17.348565328095027):1.9293316559686975,34:0.025890660448214753):0.02558760628041501,35:0.23032681135826127):1.517620578922382):1.485552033991313):0.2892609674994997):0.40363639059171064,(((((36:2.6715757112267315,37:1.1456463889778785):0.2628387034858113,(38:4.682592562005885,39:3.3426052910814743):3.02063389485524):0.5528548078579885,(((40:0.8476545078883033,41:0.45377155746685194):0.06535907593951062,(42:8.289685092033858,(43:16.105981702435876,44:3.871853871744751):2.154733134884422):0.6265712814522093):0.7958272921351153,(45:2.3068362725965823,46:0.37886020621205674):0.5426352372931778):0.11578229047289579):0.12472943033921524,47:2.1078543572447144):0.439217527162703,((48:6.982557611557448,49:6.810157913220159):2.416090828632739,50:0.6489097723645534):0.4350069513773063):1.499013117287216):4.074000881516865):0.08292494732618373;";

        TreeParser tree = new TreeParser(newickStr, false, false, true, 0);


        SkylinePopulationFunction skyline = new SkylinePopulationFunction();
        skyline.initByName(
                "tree", tree,
                "N0", new RealParameter("1.0"),
                "deltaLogPopSizes", new RealParameter("5.0 -5.0"),
                "maxSkylineIntervalCount", new IntegerParameter("5"));

        try (PrintStream ps = new PrintStream("out.txt")){
            ps.println("t N intensity intensityInv");
            double t = -1.0;
            while (t<40) {
                ps.format("%g %g %g %g\n", t,
                        skyline.getPopSize(t), skyline.getIntensity(t),
                        skyline.getInverseIntensity(skyline.getIntensity(t)));
                t += 0.1;
            }
        }
    }
}
