package lambdabsp.operators;

import beast.core.Input;
import beast.evolution.operators.SubtreeSlide;
import beast.evolution.operators.TreeOperator;
import beast.evolution.tree.Node;
import beast.evolution.tree.Tree;
import beast.util.Randomizer;
import beast.util.TreeParser;
import lambdabsp.model.CollapsedTreeIntervals;

import java.util.ArrayList;
import java.util.List;
import java.util.Random;

/**
 * Class of operators for traversing the space of multifurcating phylogenetic trees.
 */
public abstract class LambdaTreeOperator extends TreeOperator {


    protected void computeItersectionsAndTraversedNodes(Node node, double height,
                                                     List<Node> intersectingEdges,
                                                     List<Node> nodesEncountered) {

        if (!node.isRoot() && node.getParent().getHeight() < height)
            return;

        if (node.getHeight() < height) {
            intersectingEdges.add(node);
            return;
        }

        if (node.isRoot() || node.getHeight()!=node.getParent().getHeight())
            nodesEncountered.add(node);

        for (Node child : node.getChildren()) {
            computeItersectionsAndTraversedNodes(child, height,
                    intersectingEdges, nodesEncountered);
        }
    }

    /**
     * Find root of logical node tree.
     *
     * @param node member of logical node tree.
     * @return root of logical node tree.
     */
    protected Node getLogicalNodeRoot(Node node) {
        while (!node.isRoot() && node.getParent().getHeight()==node.getHeight())
            node = node.getParent();

        return node;
    }

    /**
     * Retrieve the list of children descending from this node,
     * omitting any node objects below zero-length edges.
     *
     * @param node Root of tree representing logical node
     * @return List of logical children
     */
    protected List<Node> getLogicalChildren(Node node) {

        List<Node> logicalChildren = new ArrayList<>();

        for (Node child : node.getChildren()) {
            if (child.getHeight()==node.getHeight())
                logicalChildren.addAll(getLogicalChildren(child));
            else
                logicalChildren.add(child);
        }

        return logicalChildren;
    }

    /**
     * Retrieve a list of all parents of logical nodes in the subtree
     * descending from node (which must itself be the parent of a logical node).
     *
     * @param node parent of a logical node
     * @return all parents of logical nodes descending from node (including node)
     */
    protected List<Node> getLogicalNodesInSubtree(Node node) {
        List<Node> logicalDescendents = new ArrayList<>();
        logicalDescendents.add(node);

        for (Node child : getLogicalChildren(node)) {
            logicalDescendents.addAll(getLogicalNodesInSubtree(child));
        }

        return logicalDescendents;
    }

    public static void main(String[] args) {

        TreeParser tree = new TreeParser("((A:1,B:1,C:1):1,(D:0.4,E:0.4):1.6):0.0;",
                false, false, true,1, true);

        System.out.println(tree.getRoot().toString());
        System.out.println(tree.getRoot().toNewick());


        Node root = tree.getRoot();

        List<Node> intersectingEdges = new ArrayList<>();
        List<Node> nodesEncountered = new ArrayList<>();

        LambdaTreeOperator op = new LambdaTreeOperator() {
            @Override
            public double proposal() { return 0; }

            @Override
            public void initAndValidate() { }
        };

        op.computeItersectionsAndTraversedNodes(root, 1.2,
                intersectingEdges, nodesEncountered);

        System.out.print("Intersecting edges: ");
        intersectingEdges.stream().mapToInt(Node::getNr).forEach(i -> System.out.print(i + " "));
        System.out.print("\nNodes encountered: ");
        nodesEncountered.stream().mapToInt(Node::getNr).forEach(i -> System.out.print(i + " "));

    }
}
