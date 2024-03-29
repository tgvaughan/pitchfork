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

import beast.base.core.Input;
import beast.base.evolution.operator.TreeOperator;
import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.Tree;
import beast.base.util.Randomizer;
import pitchfork.Pitchforks;

import java.util.ArrayList;
import java.util.List;

public class ExpandCollapseOperator extends TreeOperator {

    public Input<Double> rootAttachLambdaInput = new Input<>("rootAttachLambda",
            "Mean of exponential (relative to tree height) from which " +
                    "expanded node height is drawn if expanded from a " +
                    "root polytomy.",
            0.1);

    Tree tree;
    double lambda;

    @Override
    public void initAndValidate() {
        tree = treeInput.get();
        lambda = rootAttachLambdaInput.get();
    }

    @Override
    public double proposal() {

        double logHR = 0.0;

        if (Randomizer.nextBoolean()) {
            // Collapse

            List<Node> collapsableEdges = getCollapsableEdges(tree);

            if (collapsableEdges.isEmpty())
                return Double.NEGATIVE_INFINITY;

            Node edgeToCollapse = collapsableEdges.get(Randomizer.nextInt(collapsableEdges.size()));
            logHR -= Math.log(1.0/collapsableEdges.size());

            Node edgeParent = edgeToCollapse.getParent();
            Node sister = getOtherChild(edgeParent, edgeToCollapse);

            if (edgeParent.isRoot()) {
                double expRate = 1.0/(lambda*sister.getHeight());
                logHR += -expRate*(edgeParent.getHeight() - sister.getHeight()) + Math.log(expRate);
            } else {
                double L = edgeParent.getParent().getHeight() - sister.getHeight();
                logHR += Math.log(1.0/L);
            }

            edgeParent.setHeight(sister.getHeight());

            logHR += Math.log(1.0/getExpandableEdges(tree).size());

        } else {
            // Expand

            List<Node> expandableEdges = getExpandableEdges(tree);

            if (expandableEdges.isEmpty())
                return Double.NEGATIVE_INFINITY;

            Node edgeToExpand = expandableEdges.get(Randomizer.nextInt(expandableEdges.size()));
            logHR -= Math.log(1.0/expandableEdges.size());

            Node logicalParent = Pitchforks.getLogicalParent(edgeToExpand);
            assert logicalParent != null;

            double newHeight;
            if (logicalParent.isRoot()) {
                double expRate = 1.0/(lambda*logicalParent.getHeight());
                newHeight = logicalParent.getHeight() + Randomizer.nextExponential(expRate);
                logHR -= -expRate*(newHeight - logicalParent.getHeight()) + Math.log(expRate);
            } else {
                double L = logicalParent.getParent().getHeight() - logicalParent.getHeight();
                newHeight = logicalParent.getHeight() + Randomizer.nextDouble()*L;
                logHR -= Math.log(1.0/L);
            }

            // Disconnect edge

            Node nodeToMove = edgeToExpand.getParent();
            Node sisterNode = getOtherChild(nodeToMove, edgeToExpand);
            Node nodeToMoveParent = nodeToMove.getParent();

            nodeToMove.removeChild(sisterNode);
            if (nodeToMoveParent != null) {
                nodeToMoveParent.removeChild(nodeToMove);
                nodeToMoveParent.addChild(sisterNode);
            } else {
                sisterNode.setParent(null);
            }
            nodeToMove.setParent(null);

            if (nodeToMove == logicalParent)
                logicalParent = sisterNode;

            // Attach edge

            if (sisterNode.isRoot()) {
                nodeToMove.addChild(sisterNode);
                tree.setRoot(nodeToMove);
             } else {
                if (logicalParent.isRoot()) {
                    nodeToMove.addChild(logicalParent);
                    tree.setRoot(nodeToMove);
                } else {
                    Node logicalParentParent = logicalParent.getParent();
                    logicalParentParent.removeChild(logicalParent);
                    logicalParentParent.addChild(nodeToMove);
                    nodeToMove.addChild(logicalParent);
                }
            }

            // Set new node height
            nodeToMove.setHeight(newHeight);

            // Complete HR calculation

            logHR += Math.log(1.0/getCollapsableEdges(tree).size());
        }


        return logHR;
    }

    private List<Node> getCollapsableEdges(Tree tree) {

        List<Node> trueNodes = Pitchforks.getTrueNodes(tree);
        List<Node> collapsableEdges = new ArrayList<>();

        for (Node node : trueNodes) {
            if (node.isRoot() || Pitchforks.isPolytomy(node.getParent()))
                    continue;

            Node sister = getOtherChild(node.getParent(), node);
            if (sister.getHeight() < node.getHeight() || sister.isLeaf())
                continue;

            collapsableEdges.add(node);
        }

        return collapsableEdges;
    }

    private List<Node> getExpandableEdges(Tree tree) {

        List<Node> trueNodes = Pitchforks.getTrueNodes(tree);
        List<Node> expandableEdges = new ArrayList<>();

        for (Node node : trueNodes) {
            if (!node.isRoot() && Pitchforks.isPolytomy(node.getParent()))
                expandableEdges.add(node);
        }

        return expandableEdges;
    }
}
