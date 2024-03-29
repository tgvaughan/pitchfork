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
import beast.base.evolution.tree.IntervalList;
import beast.base.evolution.tree.IntervalType;
import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.Tree;
import beast.base.inference.CalculationNode;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Comparator;
import java.util.List;

/**
 * Computes tree intervals, skipping zero-length intervals.
 */
public class CollapsedTreeIntervals extends CalculationNode implements IntervalList {

    public Input<Tree> treeInput = new Input<>("tree",
            "Tree from which to compute coalescent intervals.",
            Input.Validate.REQUIRED);

    protected boolean isDirty;

    protected List<Integer> lineageCounts;
    protected List<Double> intervalDurations;
    protected int nIntervals, nSamples;
    protected boolean isBinaryTree;

    protected Tree tree;

    int nLeaves;

    @Override
    public void initAndValidate() {

        tree = treeInput.get();
        nLeaves = treeInput.get().getLeafNodeCount();

        lineageCounts = new ArrayList<>();
        intervalDurations = new ArrayList<>();
        isDirty = true;
    }

    /**
     * @return number of leaf nodes in parent tree.  Necessary because this
     * is also the maximum number of extant lineages possible.
     */
    public int getNLeaves() {
        return nLeaves;
    }

    void update() {
        if (!isDirty)
            return;

        // Do update

        lineageCounts.clear();
        intervalDurations.clear();

        /* The following looks redundant, but Arrays.asList() returns a list
        backed by the original array, which causes the sorting to mess up the
        tree's internal node array! */
        List<Node> sortedNodeList = new ArrayList<>(Arrays.asList(tree.getNodesAsArray()));
        sortedNodeList.sort(Comparator.comparingDouble(Node::getHeight));

        int lineages = 0;
        double prevHeight = 0.0;
        nIntervals = 0;
        nSamples = 0;
        isBinaryTree = true;

        boolean isFirst = true;

        for (Node node : sortedNodeList) {

            if (!isFirst) {
                double thisDuration = node.getHeight() - prevHeight;

                if (thisDuration > 0.0) {
                    lineageCounts.add(lineages);
                    intervalDurations.add(thisDuration);
                    prevHeight = node.getHeight();
                }
            } else
                isFirst = false;

            if (node.isLeaf()) {
                lineages += 1;
                nSamples += 1;
            } else
                lineages -= 1;
        }

        // Add number of lineages above root (explicitly including this makes other calculations neater)
        lineageCounts.add(1);

        nIntervals = intervalDurations.size();

        isDirty = false;
    }

    @Override
    public int getIntervalCount() {
        update();

        return nIntervals;
    }

    @Override
    public int getSampleCount() {
        update();

        return nSamples;
    }

    @Override
    public double getInterval(int i) {
        update();

        return intervalDurations.get(i);
    }

    @Override
    public int getLineageCount(int i) {
        update();

        return lineageCounts.get(i);
    }

    @Override
    public int getCoalescentEvents(int i) {
        update();

        return Math.max(0, lineageCounts.get(i) - lineageCounts.get(i+1));
    }

    @Override
    public IntervalType getIntervalType(int i) {
        update();

        if (lineageCounts.get(i+1) > lineageCounts.get(i))
            return IntervalType.SAMPLE;
        else
            return IntervalType.COALESCENT;
    }

    @Override
    public double getTotalDuration() {
        return tree.getRoot().getHeight();
    }

    @Override
    public boolean isBinaryCoalescent() {
        update();

        return isBinaryTree;
    }

    @Override
    public boolean isCoalescentOnly() {
        update();

        return nSamples == 0;
    }

    @Override
    protected void restore() {
        isDirty = true;

        super.restore();
    }

    @Override
    protected boolean requiresRecalculation() {
        isDirty = true;

        return true;
    }
}
