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

package pitchfork.util;

import beast.base.core.Function;
import beast.base.core.Input;
import beast.base.core.Loggable;
import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.Tree;
import beast.base.inference.CalculationNode;
import pitchfork.Pitchforks;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.List;

public class PolytomyStatsLogger extends CalculationNode implements Loggable, Function {

    public Input<Tree> treeInput = new Input<>("tree",
            "Tree whose polytomies to count.",
            Input.Validate.REQUIRED);

    public Input<Boolean> polytomyCountOnlyInput = new Input<>("polytomyCountOnly",
            "Only log polytomy count: no node order histogram.",
            false);

    public Input<Integer> maxOrderInput = new Input<>("maxOrder",
            "Maximum node order to include in histogram.");

    private Tree tree;
    private boolean polytomyCountOnly;
    private int[] nodeOrderHist;
    private int maxOrder;

    public PolytomyStatsLogger() { }

    @Override
    public void initAndValidate() {
        tree = treeInput.get();
        polytomyCountOnly = polytomyCountOnlyInput.get();

        if (!polytomyCountOnly) {
            maxOrder = tree.getLeafNodeCount();
            if (maxOrderInput.get() != null)
                maxOrder = Math.min(maxOrder, maxOrderInput.get());

            nodeOrderHist = new int[maxOrder - 1];
        }
    }


    @Override
    public void init(PrintStream out) {
        String prefix = getID() == null ? "" : getID() + ".";

        out.print(prefix + "PolytomyCount\t");

        if (!polytomyCountOnly) {
            for (int i = 2; i <= maxOrder; i++)
                out.print(prefix + "NodesOrder" + i + "\t");
        }
    }

    private int getPolytomyCount() {
        int count = 0;

        List<Node> trueNodes = new ArrayList<>();
        for (Node node : tree.getNodesAsArray())
            if (node.isRoot() || node.getParent().getHeight()>node.getHeight())
                trueNodes.add(node);

        for (Node node : trueNodes) {
            if (!node.isLeaf() && (node.getChildren().get(0).getHeight()==node.getHeight()
                    || node.getChildren().get(1).getHeight()==node.getHeight()))
                count += 1;
        }

        return count;
    }

    private void updateNodeOrderHist() {
        // Zero entries
        for (int i=0; i<nodeOrderHist.length; i++)
            nodeOrderHist[i] = 0;

        // Compute histogram
        List<Node> trueNodes = Pitchforks.getTrueInternalNodes(tree);

        for (Node node : trueNodes) {
            int order = Pitchforks.getLogicalChildren(node).size();
            if (order <= maxOrder)
                nodeOrderHist[order-2] += 1;
        }
    }

    @Override
    public void log(long sample, PrintStream out) {
        out.print(getPolytomyCount() + "\t");

        if (!polytomyCountOnly) {
            updateNodeOrderHist();
            for (int count : nodeOrderHist)
                out.print(count + "\t");
        }
    }

    @Override
    public void close(PrintStream out) { }

    @Override
    public int getDimension() {
        return 1 + nodeOrderHist.length;
    }

    @Override
    public double getArrayValue() {
        return getPolytomyCount();
    }

    @Override
    public double getArrayValue(int dim) {
        if (dim == 1)
            return getPolytomyCount();
        else if (dim < 1 + nodeOrderHist.length)
            return nodeOrderHist[dim-1];
        else
            return Double.NaN;
    }
}
