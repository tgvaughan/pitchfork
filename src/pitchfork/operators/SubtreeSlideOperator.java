package pitchfork.operators;

import beast.core.Description;
import beast.core.Input;
import beast.evolution.tree.Node;
import beast.evolution.tree.Tree;
import beast.util.Randomizer;
import pitchfork.util.Pitchforks;

import java.util.ArrayList;
import java.util.List;

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


        List<Node> logicalNodes = Pitchforks.getTrueNodes(treeInput.get());

        Node edgeBaseNode;
        do {
            edgeBaseNode = logicalNodes.get(Randomizer.nextInt(logicalNodes.size()));
        } while (edgeBaseNode.isRoot());

        Node edgeParentNode = Pitchforks.getLogicalParent(edgeBaseNode);
        double oldAttachmentHeight = edgeParentNode.getHeight();

        double window = relSizeInput.get()*tree.getRoot().getHeight();
        double newAttachmentHeight = oldAttachmentHeight + window*Randomizer.nextGaussian();

        if (newAttachmentHeight < edgeBaseNode.getHeight())
            return Double.NEGATIVE_INFINITY;



        List<Node> coalescenceNodes = new ArrayList<>();
        List<Node> intersectingEdges = new ArrayList<>();

        getIntersectionsAndCoalescences(edgeBaseNode, edgeParentNode, newAttachmentHeight,
                coalescenceNodes, intersectingEdges);



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
     * @param coalescenceNodes list to populate with seen coalescent nodes
     * @param intersectingEdges list to populate with (nodes below) intersecting edges
     */
    private void getIntersectionsAndCoalescences(Node baseNode, Node logicalNode, double height,
                                                List<Node> coalescenceNodes,
                                                List<Node> intersectingEdges) {

        if (height > logicalNode.getHeight()) {

            Node logicalParent = Pitchforks.getLogicalParent(logicalNode);
            if (logicalParent == null || logicalNode.getParent().getHeight() > height) {
                intersectingEdges.add(logicalNode);
            } else {
                coalescenceNodes.add(logicalParent);
                getIntersectionsAndCoalescences(baseNode, logicalParent, height,
                        coalescenceNodes, intersectingEdges);
            }

        } else {

            List<Node> logicalChildren = Pitchforks.getLogicalChildren(logicalNode);

            for (Node logicalChild : logicalChildren) {
                if (logicalChild == baseNode)
                    continue;

                if (logicalChild.getHeight() < height) {
                    intersectingEdges.add(logicalChild);
                } else {
                    coalescenceNodes.add(logicalChild);
                    getIntersectionsAndCoalescences(baseNode, logicalChild, height,
                            coalescenceNodes, intersectingEdges);
                }
            }
        }

    }

}
