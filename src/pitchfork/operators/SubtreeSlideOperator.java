package pitchfork.operators;

import beast.core.Description;
import beast.core.Input;
import beast.evolution.tree.Node;
import beast.evolution.tree.Tree;
import beast.util.Randomizer;
import pitchfork.util.Pitchforks;

import java.util.List;

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
    double probCoalAttach, relSize;

    @Override
    public void initAndValidate() {
        tree = treeInput.get();
        probCoalAttach = probCoalAttachInput.get();
        relSize = relSizeInput.get();
    }

//    int count = 0;

    @Override
    public double proposal() {

        double logHR = 0.0;

//        System.out.println("Count: " + (++count));

        // Select base node of edge to move:

        List<Node> logicalNodes = Pitchforks.getTrueNodes(tree);
        Node edgeBaseNode;
        do {
            edgeBaseNode = logicalNodes.get(Randomizer.nextInt(logicalNodes.size()));
        } while (edgeBaseNode.isRoot());

        // Forward HR contribution of edge node selection:

        logHR -= Math.log(1.0/(logicalNodes.size()-1.0));

        // Slide edge in randomly chosen direction:

        boolean isSlideUp = Randomizer.nextBoolean();

        if (isSlideUp)
            logHR += slideUp(edgeBaseNode);
        else
            logHR += slideDown(edgeBaseNode);

        // Reverse HR contribution of edge node selection:

        logHR += Math.log(1.0/(Pitchforks.getTrueNodes(tree).size()-1.0));

//        if (!isSlideUp)
//            System.out.println(logHR);

        return logHR;

    }

    double getCurrentLambda() {
        return 1.0/(relSize*tree.getRoot().getHeight());
    }

    double slideUp(Node edgeBaseNode) {

        Node edgeParentNode = edgeBaseNode.getParent();
        Node edgeSisterNode = getOtherChild(edgeParentNode, edgeBaseNode);

        double lambda = getCurrentLambda();

        AttachmentPoint newAttachmentPoint = getOlderAttachmentPoint(Pitchforks.getLogicalNode(edgeParentNode), lambda);

        AttachmentPoint oldAttachmentPoint = new AttachmentPoint();
        oldAttachmentPoint.attachmentHeight = edgeParentNode.getHeight();
        oldAttachmentPoint.attachmentEdgeBase = Pitchforks.getLogicalNode(edgeSisterNode);

        // Topology modification:

        if (edgeParentNode != newAttachmentPoint.attachmentEdgeBase) {

            Node grandParent = edgeParentNode.getParent();
            grandParent.removeChild(edgeParentNode);
            edgeParentNode.removeChild(edgeSisterNode);
            grandParent.addChild(edgeSisterNode);
            edgeParentNode.setParent(null);

            if (!newAttachmentPoint.attachmentEdgeBase.isRoot()) {
                Node newGrandParent = newAttachmentPoint.attachmentEdgeBase.getParent();
                newGrandParent.removeChild(newAttachmentPoint.attachmentEdgeBase);
                newGrandParent.addChild(edgeParentNode);
            }
            edgeParentNode.addChild(newAttachmentPoint.attachmentEdgeBase);

            if (edgeParentNode.isRoot())
                tree.setRoot(edgeParentNode);
        }
        edgeParentNode.setHeight(newAttachmentPoint.attachmentHeight);

        // Probability of reverse move:

        double lambdaPrime = getCurrentLambda();
        computeYoungerAttachmentPointProb(oldAttachmentPoint,
                edgeParentNode, lambdaPrime);

        return oldAttachmentPoint.logProb - newAttachmentPoint.logProb;
    }

    double slideDown(Node edgeBaseNode) {
        Node edgeParentNode = edgeBaseNode.getParent();
        Node edgeSisterNode = getOtherChild(edgeParentNode, edgeBaseNode);

        double lambda = getCurrentLambda();

        AttachmentPoint newAttachmentPoint;
        try {
            newAttachmentPoint = getYoungerAttachmentPoint(edgeBaseNode,
                    Pitchforks.getLogicalNode(edgeParentNode), lambda);
        } catch (AttachmentException ex) {
            return Double.NEGATIVE_INFINITY;
        }

        AttachmentPoint oldAttachmentPoint = new AttachmentPoint();
        oldAttachmentPoint.attachmentHeight = edgeParentNode.getHeight();
        oldAttachmentPoint.attachmentEdgeBase = Pitchforks.getLogicalNode(edgeSisterNode);

        // Topology modification

        if (edgeSisterNode != newAttachmentPoint.attachmentEdgeBase) {
            if (!edgeParentNode.isRoot()) {
                Node grandParent = edgeParentNode.getParent();
                grandParent.removeChild(edgeParentNode);
                edgeParentNode.removeChild(edgeSisterNode);
                grandParent.addChild(edgeSisterNode);
                edgeParentNode.setParent(null);
            } else {
                edgeParentNode.removeChild(edgeSisterNode);
                edgeSisterNode.setParent(null);
            }

            Node newGrandParent = newAttachmentPoint.attachmentEdgeBase.getParent();
            newGrandParent.removeChild(newAttachmentPoint.attachmentEdgeBase);
            newGrandParent.addChild(edgeParentNode);
            edgeParentNode.addChild(newAttachmentPoint.attachmentEdgeBase);

            if (edgeSisterNode.isRoot())
                tree.setRoot(edgeSisterNode);
        } else {
            // If topology is unchanged, node below edge supporting original
            // attachment will be the original edge parent node:

            oldAttachmentPoint.attachmentEdgeBase = edgeParentNode;
        }
        edgeParentNode.setHeight(newAttachmentPoint.attachmentHeight);

        // Probability of reverse move

        double lambdaPrime = getCurrentLambda();
        computeOlderAttachmentPointProb(oldAttachmentPoint, edgeParentNode, lambdaPrime);

        return oldAttachmentPoint.logProb - newAttachmentPoint.logProb;
    }


    class AttachmentPoint {
        Node attachmentEdgeBase;
        double attachmentHeight;
        double logProb = 0;

        @Override
        public String toString() {
            return "attachmentEdgeBase: " + attachmentEdgeBase.getNr() + ", " +
                    "attachmentHeight: " + attachmentHeight + ", " +
                    "logProb: " + logProb;
        }
    }

    class AttachmentException extends Exception { }

    AttachmentPoint getOlderAttachmentPoint(Node startNode, double lambda) {

        AttachmentPoint ap = new AttachmentPoint();

        ap.attachmentEdgeBase = startNode;
        while(true) {
            Node logicalParent = Pitchforks.getLogicalParent(ap.attachmentEdgeBase);

            if (logicalParent != null) {
                if (Randomizer.nextDouble() < probCoalAttach) {
                    ap.attachmentEdgeBase = logicalParent;
                    ap.attachmentHeight = logicalParent.getHeight();
                    ap.logProb += Math.log(probCoalAttach);
                    break;
                } else {
                    ap.logProb += Math.log(1-probCoalAttach);
                }
            }

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

    void computeOlderAttachmentPointProb(AttachmentPoint ap, Node startNode, double lambda) {

        ap.logProb = 0;

        Node currentEdgeBase = startNode;
        Node logicalParent;
        while(true) {
            logicalParent = Pitchforks.getLogicalParent(currentEdgeBase);

            if (ap.attachmentEdgeBase == currentEdgeBase)
                break;

            if (logicalParent == null) {
                ap.logProb = Double.NEGATIVE_INFINITY;
                return;
            }

            if (ap.attachmentHeight != logicalParent.getHeight())
                ap.logProb += Math.log(1-probCoalAttach) - lambda*currentEdgeBase.getLength();

            currentEdgeBase = logicalParent;
        }

        if (!currentEdgeBase.isRoot()) {
            if (ap.attachmentHeight == currentEdgeBase.getHeight()) {
                ap.logProb += Math.log(probCoalAttach);
            } else {
                ap.logProb += Math.log(1 - probCoalAttach);
            }

        }

        if (ap.attachmentHeight > currentEdgeBase.getHeight()) {
            ap.logProb += -lambda * (ap.attachmentHeight - currentEdgeBase.getHeight())
                    + Math.log(lambda);
        }
    }

    AttachmentPoint getYoungerAttachmentPoint(Node edgeBaseNode,
                                                      Node startNode,
                                                      double lambda) throws AttachmentException {

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
                ap.attachmentHeight = ap.attachmentEdgeBase.getHeight()
                        + (ap.attachmentEdgeBase.getLength() - delta);
                break;
            }

            if (ap.attachmentEdgeBase.isLeaf())
                throw new AttachmentException();

            ap.logProb += -lambda*ap.attachmentEdgeBase.getLength();
        }

        if (ap.attachmentHeight < edgeBaseNode.getHeight())
            throw new AttachmentException();

        return ap;
    }

    void computeYoungerAttachmentPointProb(AttachmentPoint ap,
                                                   Node startNode,
                                                   double lambda) {
        ap.logProb = 0.0;

        Node currentEdgeBase = ap.attachmentEdgeBase;

        do {
            if (currentEdgeBase == null)
                throw new IllegalStateException("Probability calculation loop failed to find startNode.");

            if (currentEdgeBase.getHeight() <= ap.attachmentHeight) {
                if (ap.attachmentHeight == currentEdgeBase.getHeight()) {
                    ap.logProb += Math.log(probCoalAttach);
                } else {
                    ap.logProb += Math.log(1.0 - probCoalAttach)
                            - lambda * (currentEdgeBase.getParent().getHeight() - ap.attachmentHeight)
                            + Math.log(lambda);
                }
            } else {
                ap.logProb += Math.log(1.0 - probCoalAttach)
                        - lambda*currentEdgeBase.getLength();
            }

            currentEdgeBase = Pitchforks.getLogicalParent(currentEdgeBase);

            List<Node> logicalChildren = Pitchforks.getLogicalChildren(currentEdgeBase);

            if (currentEdgeBase == startNode)
                ap.logProb += Math.log(1.0 / (logicalChildren.size() - 1));
            else
                ap.logProb += Math.log(1.0 / logicalChildren.size());

        } while (currentEdgeBase != startNode);
    }


}
