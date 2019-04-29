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

import static pitchfork.Pitchforks.getTrueNodes;
import static pitchfork.Pitchforks.isPolytomy;

@Description("SPR operator for trees with polytomies.")
public class SPROperator extends PitchforkTreeOperator {

    public Input<Double> rootAttachLambdaInput = new Input<>(
            "rootAttachLambda",
            "Mean of exponential distribution (relative to tree height)" +
                    "used to position attachments above the root.",
            2.0);

    public Input<Double> probCoalAttachInput = new Input<>(
            "probCoalAttach",
            "Probability of attaching to existing coalescent event.",
            0.1);

    Tree tree;
    Double rootAttachLambda, probCoalAttach;

    @Override
    public void initAndValidate() {
        rootAttachLambda = rootAttachLambdaInput.get();
        tree = treeInput.get();
        probCoalAttach = probCoalAttachInput.get();

        super.initAndValidate();
    }

    @Override
    boolean isSkylineSafe() {
        return false;
    }

    @Override
    public double pitchforkProposal() {

        double logHR = 0.0;

        // Get list of nodes below finite-length edges
        List<Node> trueNodes = getTrueNodes(tree);

        // Record number of (true) edges in original tree:
        int nEdges = trueNodes.size() - 1;

        // Select non-root subtree at random

        Node srcNode, srcNodeParent;
        do {
            srcNode = trueNodes.get(Randomizer.nextInt(trueNodes.size()));
            srcNodeParent = srcNode.getParent();
        } while (srcNodeParent == null);

        Node srcNodeSister = getOtherChild(srcNodeParent, srcNode);

        // Record whether the the original attachment was a polytomy
        boolean origAttachWasPolytomy = isPolytomy(srcNodeParent);

        // Incorporate probability of existing attachment point into HR

        if (origAttachWasPolytomy) {
            logHR += Math.log(probCoalAttach);
        } else {
            if (!srcNodeSister.isLeaf() && srcNodeSister.getHeight() > srcNode.getHeight())
                logHR += Math.log(1 - probCoalAttach);

            if (srcNodeParent.isRoot()) {
                double offset = Math.max(srcNodeSister.getHeight(), srcNode.getHeight());
                double expRate = 1.0/(rootAttachLambda*offset);
                logHR += -expRate*(srcNodeParent.getHeight() - offset) + Math.log(expRate);
            } else {
                double L = srcNodeParent.getParent().getHeight() - Math.max(srcNodeSister.getHeight(), srcNode.getHeight());
                logHR += Math.log(1.0/L);
            }
        }

        // Disconnect subtree

        srcNodeParent.removeChild(srcNodeSister);

        if (srcNodeParent.isRoot()) {
            srcNodeSister.setParent(null);
        } else {
            Node srcNodeGrandparent = srcNodeParent.getParent();
            srcNodeGrandparent.removeChild(srcNodeParent);
            srcNodeGrandparent.addChild(srcNodeSister);
        }

        srcNodeParent.setParent(null);

        // Select new attachment node

        Node remainingSubtreeRoot;
        if (srcNodeSister.isRoot())
            remainingSubtreeRoot = srcNodeSister;
        else
            remainingSubtreeRoot = tree.getRoot();

        List<Node> subtreeNodes = getNodesInSubtree(remainingSubtreeRoot, srcNode.getHeight());
        Node attachmentNode = subtreeNodes.get(Randomizer.nextInt(subtreeNodes.size()));


        // Determine whether polytomy is to be created

        boolean newAttachIsPolytomy;

        if (attachmentNode.isLeaf() || attachmentNode.getHeight() < srcNode.getHeight()) {
            newAttachIsPolytomy = false;
        } else {
            if (Randomizer.nextDouble() < probCoalAttach) {
                newAttachIsPolytomy = true;
                logHR -= Math.log(probCoalAttach);
            } else {
                newAttachIsPolytomy = false;
                logHR -= Math.log(1 - probCoalAttach);
            }
        }

        // Select new attachment height

        double attachmentHeight;

        if (newAttachIsPolytomy) {
            attachmentHeight = attachmentNode.getHeight();
        } else {

            if (attachmentNode.isRoot()) {
                double offset = Math.max(srcNode.getHeight(), attachmentNode.getHeight());
                double expRate = 1.0/(rootAttachLambda*offset);
                attachmentHeight = offset + Randomizer.nextExponential(expRate);

                logHR -= -expRate*(attachmentHeight-offset)
                        + Math.log(expRate);
            } else {
                double L = attachmentNode.getParent().getHeight() -
                        Math.max(srcNode.getHeight(), attachmentNode.getHeight());
                attachmentHeight = Randomizer.nextDouble()*L +
                        Math.max(srcNode.getHeight(), attachmentNode.getHeight());

                logHR -= Math.log(1.0/L);
            }
        }

        // Reconnect subtree

        srcNodeParent.setHeight(attachmentHeight);

        if (attachmentNode.isRoot()) {
            srcNodeParent.addChild(attachmentNode);
        } else {
            Node oldParent = attachmentNode.getParent();
            oldParent.removeChild(attachmentNode);
            oldParent.addChild(srcNodeParent);
            srcNodeParent.addChild(attachmentNode);
        }

        // Ensure correct root if set if this has been modified:
        if (srcNodeSister.isRoot())
            tree.setRoot(srcNodeSister);
        else if (srcNodeParent.isRoot())
            tree.setRoot(srcNodeParent);

        // Account for edge selection probability in HR:
        if (origAttachWasPolytomy != newAttachIsPolytomy) {
            if (origAttachWasPolytomy) {
                logHR += Math.log(nEdges/(nEdges+1.0));
            } else {
                logHR += Math.log(nEdges/(nEdges-1.0));
            }
        }

        return logHR;
    }

    private List<Node> getNodesInSubtree(Node subtreeRoot, double minAge) {
        List<Node> nodeList = new ArrayList<>();

        if (subtreeRoot.isRoot() || subtreeRoot.getParent().getHeight()>subtreeRoot.getHeight())
            nodeList.add(subtreeRoot);

        if (subtreeRoot.getHeight()>minAge) {
            for (Node child : subtreeRoot.getChildren())
                nodeList.addAll(getNodesInSubtree(child, minAge));
        }

        return nodeList;
    }
}
