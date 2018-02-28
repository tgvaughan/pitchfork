package pitchfork.operators;

import beast.core.Description;
import beast.core.Input;
import beast.evolution.tree.Node;
import beast.evolution.tree.Tree;
import beast.util.Randomizer;
import org.apache.commons.math.special.Erf;
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

//                logHR -= Math.log(Erf.erf(maxHeight/Math.sqrt(2)/window));

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
    private void getIntersections(Node baseNode, Node logicalNode, double height,
                                  List<Node> intersectingEdges) {

        if (height > logicalNode.getHeight()) {

            Node logicalParent = Pitchforks.getLogicalParent(logicalNode);
            if (logicalParent == null || logicalNode.getParent().getHeight() > height) {
                intersectingEdges.add(logicalNode);
            } else {
                getIntersections(baseNode, logicalParent, height, intersectingEdges);
            }

        } else {

            List<Node> logicalChildren = Pitchforks.getLogicalChildren(logicalNode);

            for (Node logicalChild : logicalChildren) {
                if (logicalChild == baseNode)
                    continue;

                if (logicalChild.getHeight() < height) {
                    intersectingEdges.add(logicalChild);
                } else if (!logicalChild.isLeaf()) {
                    getIntersections(baseNode, logicalChild, height, intersectingEdges);
                }
            }
        }
    }
}
