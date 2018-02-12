package lambdabsp.operators;

import beast.core.Input;
import beast.evolution.tree.Node;
import beast.evolution.tree.Tree;
import beast.util.Randomizer;

import java.util.ArrayList;
import java.util.List;

public class SPROperator extends LambdaTreeOperator {

    public Input<Double> rootAttachLambdaInput = new Input<>(
            "rootAttachLambda",
            "Rate of exponential distribution " +
                    "used to position attachements above the root.",
            Input.Validate.REQUIRED);

    public Input<Double> probCoalAttachInput = new Input<>(
            "probCoalAttach",
            "Probability of attaching to existing coalescent event.",
            Input.Validate.REQUIRED);

    Tree tree;
    Double rootAttachLambda, probCoalAttach;

    @Override
    public void initAndValidate() {
        rootAttachLambda = rootAttachLambdaInput.get();
        tree = treeInput.get();
        probCoalAttach = probCoalAttachInput.get();
    }

    @Override
    public double proposal() {

        double logHR = 0.0;

        // Get list of nodes below finite-length edges
        List<Node> trueNodes = getTrueNodes();

        // Record number of (true) edges in original tree:
        int nEdges = trueNodes.size() - 1;

        // Select non-root subtree at random

        Node i, ip;
        do {
            i = trueNodes.get(Randomizer.nextInt(trueNodes.size()));
            ip = i.getParent();
        } while (ip == null);

        Node is = getOtherChild(ip, i);

        // Record whether the the original attachment was a polytomy
        boolean origAttachWasPolytomy = isPolytomy(ip);

        // Incorporate probability of existing attachment point into HR

        if (origAttachWasPolytomy) {
            logHR += Math.log(probCoalAttach);
        } else {
            if (!is.isLeaf() && is.getHeight() > i.getHeight())
                logHR += Math.log(1 - probCoalAttach);

            if (ip.isRoot()) {
                logHR += -rootAttachLambda*(ip.getHeight() - Math.max(is.getHeight(),i.getHeight()))
                        + Math.log(rootAttachLambda);
            } else {
                double L = ip.getParent().getHeight() - Math.max(is.getHeight(), i.getHeight());
                logHR += Math.log(1.0/L);
            }
        }

        // Disconnect subtree

        ip.removeChild(is);

        if (ip.isRoot()) {
            is.setParent(null);
        } else {
            Node ipp = ip.getParent();
            ipp.removeChild(ip);
            ipp.addChild(is);
        }

        ip.setParent(null);

        // Select new attachment node

        Node remainingSubtreeRoot;
        if (is.isRoot())
            remainingSubtreeRoot = is;
        else
            remainingSubtreeRoot = tree.getRoot();

        List<Node> subtreeNodes = getNodesInSubtree(remainingSubtreeRoot, i.getHeight());
        Node attachmentNode = subtreeNodes.get(Randomizer.nextInt(subtreeNodes.size()));


        // Determine whether polytomy is to be created

        boolean newAttachIsPolytomy;

        if (attachmentNode.isLeaf() || attachmentNode.getHeight()>i.getHeight()) {
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
                double offset = Math.max(i.getHeight(), attachmentNode.getHeight());
                attachmentHeight = offset + Randomizer.nextExponential(rootAttachLambda);

                logHR -= -rootAttachLambda*(attachmentHeight-offset)
                        + Math.log(rootAttachLambda);
            } else {
                double L = attachmentNode.getParent().getHeight() -
                        Math.max(i.getHeight(), attachmentNode.getHeight());
                attachmentHeight = Randomizer.nextDouble()*L +
                        Math.max(i.getHeight(), attachmentNode.getHeight());

                logHR -= Math.log(1.0/L);
            }
        }

        // Reconnect subtree

        ip.setHeight(attachmentHeight);

        if (attachmentNode.isRoot()) {
            ip.addChild(attachmentNode);
        } else {
            Node oldParent = attachmentNode.getParent();
            oldParent.removeChild(attachmentNode);
            oldParent.addChild(ip);
            ip.addChild(attachmentNode);
        }

        // Ensure correct root if set if this has been modified:
        if (is.isRoot())
            tree.setRoot(is);
        else if (ip.isRoot())
            tree.setRoot(ip);

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
