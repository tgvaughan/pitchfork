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

import beast.core.Input;
import beast.core.parameter.RealParameter;
import beast.core.util.Log;
import beast.evolution.tree.Node;
import beast.evolution.tree.Tree;
import beast.math.distributions.ParametricDistribution;
import beast.util.Randomizer;

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

    public Input<Integer> maxSkylineIntervalsInput = new Input<>(
            "maxSkylineIntervals",
            "Maximum number of skyline intervals.",
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

    private int maxSkylineIntervals, mmin;
    private ParametricDistribution skylinePopDistr;
    private BetaCoalescentModel betaCoalescentModel;

    public SimulatedBetaSkylineTree() { }

    @Override
    public void initAndValidate() {
        super.initAndValidate();

        maxSkylineIntervals = maxSkylineIntervalsInput.get();
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

        initArrays();
        simulate();

        // Stretch time in coalescent intervals according to population sizes

        mmin = getInternalNodeCount()/maxSkylineIntervals;
        if (getInternalNodeCount() % maxSkylineIntervals > 0)
            mmin += 1;

        // TODO: implement time stretch

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

    /**
     * Perform coalescent simulation.
     */
    private void simulate() {

        List<Node> activeLineages = new ArrayList<>();
        List<Node> unusedLineages = new ArrayList<>();
        List<Node> coalescingLineages = new ArrayList<>();

        // Create leaf nodes:
        for (int i=0; i<nLeaves; i++)  {
            Node leaf = new Node(leafNames[i]);
            leaf.setNr(i);
            leaf.setHeight(leafAges[i]);
            unusedLineages.add(leaf);
        }

        int nextInternalNr = nLeaves;

        unusedLineages.sort(Comparator.comparingDouble(Node::getHeight));

        double t = 0;
        while (!unusedLineages.isEmpty() || activeLineages.size() > 1) {

            // Compute propensities

            int n = activeLineages.size();
            double totalPropensity = n>=2 ? betaCoalescentModel.getTotalCoalRate(n) : 0.0 ;

            // Increment time
            t += Randomizer.nextExponential(totalPropensity);

            // Check whether next sample time exceeded.
            if (!unusedLineages.isEmpty() && t>unusedLineages.get(0).getHeight()) {
                t = unusedLineages.get(0).getHeight();
                activeLineages.add(unusedLineages.remove(0));
                continue;
            }

            // Choose reaction
            double u = Randomizer.nextDouble()*totalPropensity;
            int k = 2 + -(Arrays.binarySearch(betaCoalescentModel.getCumulativeCoalRateArray(n), u) + 1);

            if (k>n || k<2)
                throw new IllegalStateException("Numerical error in furcation degree sampler.");


            // Implement coalescence
            // (Note: BEAST really only deals with binary trees, so have to
            // fake multifurcations using zero-length edges.)

            Node newParent = new Node(String.valueOf(nextInternalNr));
            newParent.setNr(nextInternalNr++);
            newParent.setHeight(t);

            int[] indices = Randomizer.shuffled(n);

            coalescingLineages.clear();
            for (int i=0; i<k; i++)
                coalescingLineages.add(activeLineages.get(indices[i]));

            newParent.addChild(coalescingLineages.get(0));
            newParent.addChild(coalescingLineages.get(1));

            for (int i=2; i<k; i++) {
                Node newNewParent = new Node(String.valueOf(nextInternalNr));
                newNewParent.setNr(nextInternalNr++);
                newNewParent.setHeight(t);

                newNewParent.addChild(newParent);
                newNewParent.addChild(coalescingLineages.get(i));

                newParent = newNewParent;
            }

            activeLineages.removeAll(coalescingLineages);

            activeLineages.add(newParent);
        }

        assignFromWithoutID(new Tree(activeLineages.get(0)));
    }
}
