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

    class EdgePosition {
        Node baseNode;
        double heightAboveBase;
    }

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

        // Select non-root subtree at random

        Node i, ip;
        do {
            i = trueNodes.get(Randomizer.nextInt(trueNodes.size()));
            ip = i.getParent();
        } while (ip == null);

        Node is = getOtherChild(ip, i);

        // Incorporate probability of attachment point into HR

        if (isPolytomy(ip)) {
            logHR += Math.log(probCoalAttach);
        } else {
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


        // Select attachment height:

        double attachmentHeight;
        if (!attachmentNode.isLeaf()
                && attachmentNode.getHeight()>i.getHeight()
                && Randomizer.nextDouble() < probCoalAttach) {
            attachmentHeight = attachmentNode.getHeight();

            logHR -= Math.log(probCoalAttach);

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

            logHR -= Math.log(1.0 - probCoalAttach);
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
