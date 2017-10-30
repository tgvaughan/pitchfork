package lambdabsp.operators;

import beast.core.Input;
import beast.evolution.tree.Node;
import beast.evolution.tree.Tree;
import beast.util.Randomizer;

import java.util.ArrayList;
import java.util.List;

public class EdgeCollapseOperator extends LambdaTreeOperator {

    public Input<Double> alphaInput = new Input<>(
            "alpha",
            "Mean of exponentially distributed length of expanded edge, " +
                    "relative to child node height.",
            Input.Validate.REQUIRED);

    double alpha;

    @Override
    public void initAndValidate() {
        alpha = alphaInput.get();
    }

    @Override
    public double proposal() {
        double logHR = 0.0;

        Tree tree = treeInput.get();

        // Choose anchor node
        List<Node> logicalNodes = getLogicalNodesInSubtree(treeInput.get().getRoot(), false);
        int nInternalNodes = logicalNodes.size();
        logHR -= -Math.log(nInternalNodes-1);

        Node anchorNode;
        do {
            anchorNode = logicalNodes.get(Randomizer.nextInt(nInternalNodes));
        } while (anchorNode.isRoot());

        // Choose whether to expand or contract
        if (Randomizer.nextBoolean()) {
            // Expand

            logHR += doExpand(anchorNode);
            logHR += -Math.log(nInternalNodes+1);

        } else {
            // Contract

            logHR += doContract(anchorNode);
            logHR += Math.log(nInternalNodes-1);
        }

        return logHR;
    }

    double doExpand(Node anchorNode) {
        double logHR = 0.0;

        List<Node> anchorSibs = getLogicalChildren(anchorNode.getParent());
        anchorSibs.remove(anchorNode);

        List<Node> upNodes = new ArrayList<>();
        List<Node> downNodes = new ArrayList<>();

        // Choose edges to move up or leave down.

        for (Node child : anchorSibs) {
            if (Randomizer.nextBoolean())
                upNodes.add(child);
            else
                downNodes.add(child);
        }

        logHR -= anchorSibs.size()*Math.log(0.5);

        if (upNodes.isEmpty() || downNodes.isEmpty())
            return Double.NEGATIVE_INFINITY;

        // Choose length of new edge:
        double newLength = Randomizer.nextExponential(1.0/(anchorNode.getHeight()*alpha));
        logHR -= -newLength/(anchorNode.getHeight()*alpha);

        // Disconnect up nodes
        for (Node upNode : upNodes) {
            disconnectSubtree(upNode, anchorNode);
        }

        // Find new representative of logical node formed by remaining down nodes
        Node logicalDownRoot = getLogicalNodeRoot(downNodes.get(0).getParent());

        if (logicalDownRoot.getHeight()+newLength>logicalDownRoot.getParent().getHeight())
            return Double.NEGATIVE_INFINITY;

        // Reattach up nodes at chosen height
        for (Node upNode : upNodes) {
            connectSubtree(upNode, logicalDownRoot, logicalDownRoot.getHeight()+newLength);
        }

        return logHR;
    }

    void disconnectSubtree(Node logicalChild, Node anchorNode) {
       Node actualParent = logicalChild.getParent();
       Node actualGrandParent = actualParent.getParent();
       Node actualSib = getOtherChild(actualParent, logicalChild);

       actualGrandParent.removeChild(actualParent);
       actualParent.removeChild(actualSib);
       actualGrandParent.addChild(actualSib);
    }

    void connectSubtree(Node logicalChild, Node anchorNode, double height) {

        Node actualParent = logicalChild.getParent();
        Node anchorNodeParent = anchorNode.getParent();

        anchorNodeParent.removeChild(anchorNode);
        anchorNodeParent.addChild(actualParent);
        actualParent.addChild(anchorNode);

        actualParent.setHeight(height);
    }

    double doContract(Node anchorNode) {
        double logHR = 0.0;

        

        return logHR;
    }
}
