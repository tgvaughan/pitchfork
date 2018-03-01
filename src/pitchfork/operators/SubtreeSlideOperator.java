package pitchfork.operators;

import beast.core.Description;
import beast.core.Input;
import beast.evolution.tree.Node;
import beast.evolution.tree.Tree;
import beast.util.Randomizer;
import pitchfork.util.Pitchforks;

import java.util.ArrayList;
import java.util.List;

import static java.lang.Math.PI;

@Description("Implements a version of BEAST's subtree slide operator which " +
        "is applicable to trees with hard polytomies.")
public class SubtreeSlideOperator extends PitchforkTreeOperator {

    public Input<Double> relSizeInput = new Input<>("relSize",
            "Size of slide window, relative to tree height.",
            0.15);

    public Input<Double> probCoalAttachInput = new Input<>("probCoalAttach",
            "Probability of attaching to the nearest coalescent node following slide.",
            0.1);

    Tree tree;
    double probCoalAttach;

    @Override
    public void initAndValidate() {
        tree = treeInput.get();
        probCoalAttach = probCoalAttachInput.get();
    }


    @Override
    public double proposal() {

        double logHR = 0.0;

        // Select base node of edge to move:

        List<Node> logicalNodes = Pitchforks.getTrueNodes(treeInput.get());
        Node edgeBaseNode;
        do {
            edgeBaseNode = logicalNodes.get(Randomizer.nextInt(logicalNodes.size()));
        } while (edgeBaseNode.isRoot());

        Node edgeParentNode = edgeBaseNode.getParent();
        Node edgeBaseSister = getOtherChild(edgeParentNode, edgeBaseNode);

        boolean wasPolytomy = Pitchforks.isPolytomy(edgeParentNode);

        // Choose new edge attachment height:

        double oldAttachmentHeight = edgeParentNode.getHeight();

        double window = relSizeInput.get()*tree.getRoot().getHeight();
        double deltaHeight = Randomizer.nextExponential(1.0/window);
        if (Randomizer.nextBoolean())
            deltaHeight = -deltaHeight;

        double newAttachmentHeight = oldAttachmentHeight + deltaHeight;

        // Avoid illegal height changes:

        if (newAttachmentHeight < edgeBaseNode.getHeight())
            return Double.NEGATIVE_INFINITY;

        List<Node> intersectingEdges = new ArrayList<>();
        List<Node> coalescentNodes = new ArrayList<>();

        getIntersections(edgeBaseNode, edgeParentNode, newAttachmentHeight, intersectingEdges);

        if (intersectingEdges.isEmpty())
            return Double.NEGATIVE_INFINITY;

        Node newEdgeBaseSister = intersectingEdges.get(Randomizer.nextInt(intersectingEdges.size()));

        logHR -= Math.log(1.0/intersectingEdges.size());

        if (Randomizer.nextDouble()<probCoalAttach) {

            if (deltaHeight>0) {
                newAttachmentHeight = newEdgeBaseSister.getParent().getHeight();

                double maxHeight = newEdgeBaseSister.getParent().getHeight();
                double minHeight = newEdgeBaseSister.getHeight();

                logHR -= Math.log(Math.exp(-minHeight/window) - Math.exp(-maxHeight/window));

            } else {
                newAttachmentHeight = newEdgeBaseSister.getHeight();
            }

            logHR -= Math.log(probCoalAttach);
        } else {
            logHR -= Math.log(1-probCoalAttach)
                    - deltaHeight/window + Math.log(1.0/window);
        }

        // Perform required topology adjustments

        if (newEdgeBaseSister != edgeBaseSister && newEdgeBaseSister != edgeParentNode) {
            edgeParentNode.removeChild(edgeBaseSister);

            if (!edgeParentNode.isRoot()) {
                Node grandParent = edgeParentNode.getParent();
                grandParent.removeChild(edgeParentNode);
                edgeParentNode.setParent(null);
                grandParent.addChild(edgeBaseSister);
            } else {
                edgeBaseSister.setParent(null);
            }

            if (!newEdgeBaseSister.isRoot()) {
                Node newEdgeBaseSisterParent = newEdgeBaseSister.getParent();
                newEdgeBaseSisterParent.removeChild(newEdgeBaseSister);
                newEdgeBaseSisterParent.addChild(edgeParentNode);
            }
            edgeParentNode.addChild(newEdgeBaseSister);

            if (edgeParentNode.isRoot())
                tree.setRoot(edgeParentNode);
            else if (edgeBaseSister.isRoot())
                tree.setRoot(edgeBaseSister);
        }

        edgeParentNode.setHeight(newAttachmentHeight);

        // Incorporate probability of reverse move into HR

        intersectingEdges.clear();
        coalescentNodes.clear();
        getIntersections(edgeBaseNode, edgeParentNode,
                oldAttachmentHeight, intersectingEdges);

        logHR += Math.log(1.0/intersectingEdges.size());

        double reverseWindow = relSizeInput.get()*tree.getRoot().getHeight();
        logHR += -deltaHeight/reverseWindow + Math.log(1.0/reverseWindow);

        return logHR;
    }

    private class AttachmentPoint {
        Node attachmentEdgeBase;
        double attachmentHeight;
        double logProb = 0;
    }

    private class AttachmentException extends Exception { };

    private AttachmentPoint getOlderAttachmentPoint(Node startNode, double lambda) {

        AttachmentPoint ap = new AttachmentPoint();

        ap.attachmentEdgeBase = startNode;
        while(true) {
            Node logicalParent = Pitchforks.getLogicalParent(ap.attachmentEdgeBase);

            if (logicalParent != null && Randomizer.nextDouble() < probCoalAttach) {
                ap.attachmentHeight = logicalParent.getHeight();
                ap.logProb += Math.log(probCoalAttach);
                break;
            }

            ap.logProb += Math.log(1-probCoalAttach);
            double delta = Randomizer.nextExponential(lambda);

            if (logicalParent == null || delta < ap.attachmentEdgeBase.getLength()) {
                ap.logProb += -lambda*delta + Math.log(lambda);
                ap.attachmentHeight = ap.attachmentEdgeBase.getHeight() + delta;
                break;
            }

            ap.logProb += -lambda*ap.attachmentEdgeBase.getLength();
            ap.attachmentEdgeBase = logicalParent;
        }

        return ap;
    }

    private double computeOlderAttachmentPointProb(AttachmentPoint ap, Node startNode, double lambda) {

        double logProb = 0.0;

        Node currentEdgeBase = startNode;
        Node logicalParent;
        while(true) {
            logicalParent = Pitchforks.getLogicalParent(currentEdgeBase);

            if (ap.attachmentEdgeBase == currentEdgeBase)
                break;

            if (logicalParent == null)
                return Double.NEGATIVE_INFINITY;

            logProb += Math.log(1-probCoalAttach) - lambda*currentEdgeBase.getLength();
            currentEdgeBase = logicalParent;
        }

        if (ap.attachmentHeight == logicalParent.getHeight()) {
            logProb += Math.log(probCoalAttach);
        } else {
            logProb += Math.log(1-probCoalAttach)
                    - lambda*(ap.attachmentHeight-currentEdgeBase.getHeight())
                    + Math.log(lambda);
        }

        return logProb;
    }

    private AttachmentPoint getYoungerAttachmentPoint(Node edgeBaseNode, Node startNode, double lambda) throws AttachmentException {

        AttachmentPoint ap = new AttachmentPoint();

        ap.attachmentEdgeBase = startNode;
        while(true) {
            List<Node> logicalChildren = Pitchforks.getLogicalChildren(ap.attachmentEdgeBase);
            if (ap.attachmentEdgeBase == startNode)
                logicalChildren.remove(edgeBaseNode);

            ap.attachmentEdgeBase = Pitchforks.randomChoice(logicalChildren);
            ap.logProb += Math.log(1.0/logicalChildren.size());

            if (Randomizer.nextDouble() < probCoalAttach) {
                ap.attachmentHeight = ap.attachmentEdgeBase.getHeight();
                ap.logProb += Math.log(probCoalAttach);
                break;
            }

            ap.logProb += Math.log(1-probCoalAttach);
            double delta = Randomizer.nextExponential(lambda);

            if (delta < ap.attachmentEdgeBase.getLength()) {
                ap.logProb += -lambda*delta + Math.log(lambda);
                ap.attachmentHeight = ap.attachmentEdgeBase.getHeight() + (ap.attachmentEdgeBase.getLength() - delta);
                break;
            }

            if (ap.attachmentEdgeBase.isLeaf())
                throw new AttachmentException();

            ap.logProb += -lambda*ap.attachmentEdgeBase.getLength();
        }

        return ap;
    }
}
