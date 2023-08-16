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

package pitchfork.models;

import beast.base.core.Input;
import beast.base.core.Log;
import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.Tree;
import beast.base.inference.distribution.ParametricDistribution;
import beast.base.inference.parameter.RealParameter;
import beast.base.util.Randomizer;
import org.apache.commons.math.MathException;

import java.io.FileNotFoundException;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Comparator;
import java.util.List;

/**
 * Simulator for multifurcating skyline trees.
 */
public class SimulatedBetaSkylineTree extends Tree {

    public Input<BetaCoalescentModel> lcModelInput = new Input<>(
            "model",
            "Beta coalescent model.",
            Input.Validate.REQUIRED);

    public Input<RealParameter> skylinePopulationsInput = new Input<>(
            "skylinePopulations",
            "Parameter containing skyline population sizes.",
            Input.Validate.REQUIRED);

    public Input<ParametricDistribution> skylinePopDistrInput = new Input<>(
            "skylinePopDistr",
            "Skyline population distribution. (Must be proper.)",
            Input.Validate.REQUIRED);

    public Input<String> fileNameInput = new Input<>(
            "fileName",
            "Name of file to save Newick representation of tree to.");

    private double[] leafAges;
    private String[] leafNames;
    private int nLeaves;

    private ParametricDistribution skylinePopDistr;
    private BetaCoalescentModel betaCoalescentModel;
    private RealParameter skylinePopulations;

    public SimulatedBetaSkylineTree() { }

    @Override
    public void initAndValidate() {
        super.initAndValidate();

        skylinePopulations = skylinePopulationsInput.get();
        skylinePopDistr = skylinePopDistrInput.get();

        if (!hasDateTrait()) {
            throw new RuntimeException("Must specify tip dates!");
        }

        nLeaves = getDateTrait().taxaInput.get().getTaxonCount();
        leafAges = new double[nLeaves];
        leafNames = new String[nLeaves];

        for (int nodeNr=0; nodeNr<nLeaves; nodeNr++) {
            String taxonName = getDateTrait().taxaInput.get().getTaxonId(nodeNr);

            leafAges[nodeNr] = getDateTrait().getValue(taxonName);
            leafNames[nodeNr] = taxonName;
        }

        betaCoalescentModel = lcModelInput.get();

        List<Event> treeEvents = simulateTreeEvents();

        initArrays();

        simulateTree(treeEvents);

        // Write output file
        if (fileNameInput.get() != null) {
            try (PrintStream ps = new PrintStream(fileNameInput.get())) {
                ps.println(getRoot().toNewick().concat(";"));
            } catch (FileNotFoundException ex) {
                Log.err.println("Could not write to output file.");
                System.exit(1);
            }
        }
    }

    private abstract class Event {
        double time;
        int linageCount;

        abstract int getMultiplicity();

        Event(double time, int lineageCount) {
            this.time = time;
            this.linageCount = lineageCount;
        }
    }

    private class SampEvent extends Event {
        List<Node> leafNodes;

        SampEvent(double time, int lineageCount, List<Node> leafNodes) {
            super(time, lineageCount);
            this.leafNodes = leafNodes;
        }

        @Override
        int getMultiplicity() {
            return leafNodes.size();
        }
    }

    private class CoalEvent extends Event {
        int multiplicity;

        CoalEvent(double time, int lineageCount, int multiplicity) {
            super(time, lineageCount);

            this.multiplicity = multiplicity;
        }

        @Override
        public int getMultiplicity() {
            return multiplicity;
        }
    }

    /**
     * Simulate event sequence
     */
    private List<Event> simulateTreeEvents() {
        List<Node> sampleNodes = new ArrayList<>();

        // Create leaf nodes:
        for (int i=0; i<nLeaves; i++)  {
            Node leaf = new Node(leafNames[i]);
            leaf.setNr(i);
            leaf.setHeight(leafAges[i]);
            sampleNodes.add(leaf);
        }

        sampleNodes.sort(Comparator.comparingDouble(Node::getHeight));

        List<SampEvent> sampleEvents = new ArrayList<>();

        SampEvent currentEvent = null;
        for (Node sampleNode : sampleNodes) {
            if (currentEvent == null || sampleNode.getHeight()>currentEvent.time) {
                currentEvent = new SampEvent(sampleNode.getHeight(), -1, new ArrayList<>());
                currentEvent.leafNodes.add(sampleNode);
                sampleEvents.add(currentEvent);
            } else {
                currentEvent.leafNodes.add(sampleNode);
            }
        }

        List<Event> treeEvents = new ArrayList<>();

        int n = 0;

        double t = 0;
        while (!sampleEvents.isEmpty() || n > 1) {

            // Compute propensities

            double totalPropensity = n>=2 ? betaCoalescentModel.getTotalCoalRate(n) : 0.0 ;

            // Increment time
            t += Randomizer.nextExponential(totalPropensity);

            // Check whether next sample time exceeded.
            if (!sampleEvents.isEmpty() && t>sampleEvents.get(0).time) {
                SampEvent sampEvent = sampleEvents.get(0);
                t = sampEvent.time;
                sampEvent.linageCount = n;
                n += sampEvent.getMultiplicity();
                treeEvents.add(sampEvent);
                sampleEvents.remove(0);

                continue;
            }

            // Choose reaction
            double u = Randomizer.nextDouble()*totalPropensity;
            int k = 2 + -(Arrays.binarySearch(betaCoalescentModel.getCumulativeCoalRateArray(n), u) + 1);

            if (k>n || k<2)
                throw new IllegalStateException("Numerical error in furcation degree sampler.");


            // Implement coalescence
            CoalEvent coalEvent = new CoalEvent(t, n, k-1);
            n -= k;

            treeEvents.add(coalEvent);
        }

        return treeEvents;
    }

    private void simulateTree(List<Event> treeEvents) {

        int nextInternalNr = nLeaves;

        // Count total number of coalescent intervals
        int nCoalescentIntervals = (int)treeEvents.stream().filter(e -> e instanceof CoalEvent).count();

        int maxSkylineIntervals = skylinePopulations.getDimension();

        int mmin = nCoalescentIntervals/maxSkylineIntervals;
        if (getInternalNodeCount() % maxSkylineIntervals > 0)
            mmin += 1;

        int groupIdx=0;
        drawPopulationSize(groupIdx);


        int coalIntervalIdx = 0;

        List<Node> activeLineages = new ArrayList<>();
        List<Node> coalescingLineages = new ArrayList<>();

        double t = 0.0;

        for (int i=0; i<treeEvents.size(); i++) {

            double dt = i>0 ? treeEvents.get(i).time-treeEvents.get(i-1).time : 0.0;
            t += dt*skylinePopulations.getArrayValue(groupIdx);

            if (treeEvents.get(i) instanceof SampEvent) {
                SampEvent sampEvent = (SampEvent)treeEvents.get(i);

                for (Node leafNode : sampEvent.leafNodes) {
                    leafNode.setHeight(t);
                    activeLineages.add(leafNode);
                }

            } else {
                CoalEvent coalEvent = (CoalEvent)treeEvents.get(i);
                int k = coalEvent.getMultiplicity() + 1;

                Node newParent = new Node(String.valueOf(nextInternalNr));
                newParent.setNr(nextInternalNr++);
                newParent.setHeight(t);

                int[] indices = Randomizer.shuffled(activeLineages.size());

                coalescingLineages.clear();
                for (int l=0; l<k; l++)
                    coalescingLineages.add(activeLineages.get(indices[l]));

                newParent.addChild(coalescingLineages.get(0));
                newParent.addChild(coalescingLineages.get(1));

                for (int l=2; l<k; l++) {
                    Node newNewParent = new Node(String.valueOf(nextInternalNr));
                    newNewParent.setNr(nextInternalNr++);
                    newNewParent.setHeight(t);

                    newNewParent.addChild(newParent);
                    newNewParent.addChild(coalescingLineages.get(l));

                    newParent = newNewParent;
                }

                activeLineages.removeAll(coalescingLineages);

                activeLineages.add(newParent);

                // Switch to next pop size group if necessary

                if (coalIntervalIdx/mmin > groupIdx && (nCoalescentIntervals - (coalIntervalIdx + 1)) > mmin) {
                    groupIdx += 1;
                    drawPopulationSize(groupIdx);
                }

                coalIntervalIdx += 1;
            }
        }

        assignFromWithoutID(new Tree(activeLineages.get(0)));
    }

    /**
     * Sample and set a new coalescent interval group population size.
     *
     * @param groupIdx index of coalescent interval group
     */
    private void drawPopulationSize(int groupIdx) {
        try {
            skylinePopulations.setValue(groupIdx, skylinePopDistr.sample(1)[0][0]);
        } catch (MathException e) {
            e.printStackTrace();
        }
    }
}
