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

package pitchfork.operators;

import beast.core.Description;
import beast.core.Input;
import beast.evolution.tree.Node;
import beast.evolution.tree.Tree;
import beast.util.Randomizer;

import java.util.ArrayList;
import java.util.List;

import static pitchfork.Pitchforks.getGroupAndLogicalChildren;
import static pitchfork.Pitchforks.getTrueInternalNodes;

@Description("Uniform node height operator compatible with trees having polytomies.")
public class UniformOperator extends PitchforkTreeOperator {

    public Input<Boolean> scaleRootInput = new Input<>(
            "scaleRoot",
            "Whether to scale the age of the root node.",
            true);

    public Input<Double> scaleFactorInput = new Input<>(
            "scaleFactor",
            "Tuning parameter for scaling root.",
            0.8);

    Tree tree;

    boolean scaleRoot;
    double scaleFactor;

    @Override
    public void initAndValidate() {
        tree = treeInput.get();
        scaleRoot = scaleRootInput.get();
        scaleFactor = scaleFactorInput.get();

        super.initAndValidate();
    }

    @Override
    boolean isSkylineSafe() {
        return true;
    }

    @Override
    public double pitchforkProposal() {

        double logHR = 0.0;

        List<Node> trueNodes = getTrueInternalNodes(tree);

        if (trueNodes.size() == 1 && !scaleRoot)
            return Double.NEGATIVE_INFINITY;

        Node logicalNode;
        do {
            logicalNode = trueNodes.get(Randomizer.nextInt(trueNodes.size()));
        } while (!scaleRoot && logicalNode.isRoot());


        List<Node> nodesInLogicalGroup = new ArrayList<>();
        List<Node> logicalChildren = new ArrayList<>();
        getGroupAndLogicalChildren(logicalNode, nodesInLogicalGroup, logicalChildren);
        double minHeight = logicalChildren.stream().mapToDouble(Node::getHeight).max().getAsDouble();

        double newHeight;
        if (logicalNode.isRoot()) {

            double minf = Math.min(scaleFactor, 1.0/scaleFactor);
            double maxf = 1.0/minf;
            double f = minf + Randomizer.nextDouble()*(maxf - minf);

            newHeight = logicalNode.getHeight() * f;

            if (newHeight < minHeight)
                return Double.NEGATIVE_INFINITY;

            logHR = -Math.log(f);

        } else {
            Node parent = logicalNode.getParent();
            double maxHeight = parent.getHeight();

            newHeight = minHeight + (maxHeight - minHeight) * Randomizer.nextDouble();
        }

        logicalNode.setHeight(newHeight);
        for (Node node : nodesInLogicalGroup)
            node.setHeight(newHeight);

        return logHR;
    }
}
