package lambdabsp.operators;

import beast.core.Input;
import beast.core.parameter.RealParameter;
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

    int count = 0;

    @Override
    public double proposal() {

        count += 1;
        System.out.println("\nCount: " + count);

        // Select non-root subtree at random

        Node i, ip;
        do {
            i = tree.getNode(Randomizer.nextInt(tree.getNodeCount()-1));
            ip = i.getParent();
        } while (ip.getHeight() == i.getHeight());

        // Disconnect subtree
        tree.startEditing(this);

        Node is = getOtherChild(ip, i);
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

        List<Node> subtreeNodes = getNodesInSubtree(tree.getRoot(), i.getHeight());
        Node attachmentNode = subtreeNodes.get(Randomizer.nextInt(subtreeNodes.size()));


        // Select attachment height:

        double attachmentHeight;
        if (attachmentNode.isRoot()) {
            System.out.println("Rerooting");
            attachmentHeight = Math.max(i.getHeight(), attachmentNode.getHeight())
                    + Randomizer.nextExponential(rootAttachLambda);
        } else {
            if (!attachmentNode.isLeaf()
                    && attachmentNode.getHeight()>i.getHeight()
                    && Randomizer.nextDouble() < probCoalAttach) {
                System.out.println("Creating polytomy");
                attachmentHeight = attachmentNode.getHeight();
            } else {
                System.out.println("Regular edge attachment");
                double branchLength = attachmentNode.getParent().getHeight() -
                        Math.max(i.getHeight(), attachmentNode.getHeight());
                attachmentHeight = Randomizer.nextDouble()*branchLength +
                        Math.max(i.getHeight(), attachmentNode.getHeight());
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

        return 0;
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
