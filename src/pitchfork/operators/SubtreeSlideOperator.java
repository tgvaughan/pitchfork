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
            0.1);

    public Input<Double> probCoalAttachInput = new Input<>("probCoalAttach",
            "Probability of attaching to any given coalescent node within" +
                    "slide window.",
            0.1);

    Tree tree;

    @Override
    public void initAndValidate() {
        tree = treeInput.get();
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

        // Avoid polytomy deletion:

        if (Pitchforks.isPolytomy(edgeParentNode))
            return Double.NEGATIVE_INFINITY;

        // Choose new edge attachment height:

        double oldAttachmentHeight = edgeParentNode.getHeight();
        double window = relSizeInput.get()*tree.getRoot().getHeight();
        double deltaHeight = window*Randomizer.nextGaussian();
        double newAttachmentHeight = oldAttachmentHeight + deltaHeight;
        logHR -= -deltaHeight*deltaHeight/(2*window*window)
                - 0.5*Math.log(2*PI*window*window);

        // Avoid illegal height changes:

        if (newAttachmentHeight < edgeBaseNode.getHeight())
            return Double.NEGATIVE_INFINITY;

        List<Node> intersectingEdges = new ArrayList<>();
        List<Node> coalescentNodes = new ArrayList<>();

        getIntersectionsAndCoalescences(edgeBaseNode, edgeParentNode,
                newAttachmentHeight, intersectingEdges, coalescentNodes);

        if (intersectingEdges.isEmpty())
            return Double.NEGATIVE_INFINITY;

        Node newEdgeBaseSister = intersectingEdges.get(Randomizer.nextInt(intersectingEdges.size()));

        logHR -= Math.log(1.0/intersectingEdges.size());

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
        getIntersectionsAndCoalescences(edgeBaseNode, edgeParentNode,
                oldAttachmentHeight, intersectingEdges, coalescentNodes);

        logHR += Math.log(1.0/intersectingEdges.size());

        double reverseWindow = relSizeInput.get()*tree.getRoot().getHeight();
        logHR += -deltaHeight*deltaHeight/(2*reverseWindow*reverseWindow)
                - 0.5*Math.log(2*PI*reverseWindow*reverseWindow);

        return logHR;
    }


    /**
     * Given a starting node and a height, assemble a list of (nodes below)
     * directly ancestral or descendant edges which intersect that height
     * and a list of all directly ancestral/descendant coalescent nodes
     * passed along the way.
     *
     * @param baseNode node to avoid when traversing downwards
     * @param logicalNode starting node
     * @param height height at which to detect intersections.
     * @param intersectingEdges list to populate with (nodes below) intersecting edges
     */
    private void getIntersectionsAndCoalescences(Node baseNode, Node logicalNode, double height,
                                                List<Node> intersectingEdges,
                                                 List<Node> coalescentNodes) {

        if (height > logicalNode.getHeight()) {

            Node logicalParent = Pitchforks.getLogicalParent(logicalNode);
            if (logicalParent == null || logicalNode.getParent().getHeight() > height) {
                intersectingEdges.add(logicalNode);
            } else {
                coalescentNodes.add(logicalParent);
                getIntersectionsAndCoalescences(baseNode, logicalParent, height,
                        intersectingEdges, coalescentNodes);
            }

        } else {

            List<Node> logicalChildren = Pitchforks.getLogicalChildren(logicalNode);

            for (Node logicalChild : logicalChildren) {
                if (logicalChild == baseNode)
                    continue;

                if (logicalChild.getHeight() < height) {
                    intersectingEdges.add(logicalChild);
                } else if (!logicalChild.isLeaf()) {
                    coalescentNodes.add(logicalChild);
                    getIntersectionsAndCoalescences(baseNode, logicalChild, height,
                            intersectingEdges, coalescentNodes);
                }
            }
        }
    }
}
