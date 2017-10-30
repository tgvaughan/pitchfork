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
     * Find all nodes descending from node which are part of the same logical
     * node.
     *
     * @param node start of traversal
     * @return List of nodes belonging to same logical node as "node"
     */
    protected List<Node> getIndividualNodesInLogicalNode(Node node) {
        List<Node> nodeList = new ArrayList<>();
        nodeList.add(node);

        for (Node child : node.getChildren()) {
            if (child.getHeight() == node.getHeight())
                nodeList.addAll(getIndividualNodesInLogicalNode(child));
        }

        return nodeList;
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
     * @param internalOnly if true, only internal logical nodes are included
     * @return all parents of logical nodes descending from node (including node)
     */
    protected List<Node> getLogicalNodesInSubtree(Node node, boolean internalOnly) {
        List<Node> logicalDescendents = new ArrayList<>();
        logicalDescendents.add(node);

        for (Node child : getLogicalChildren(node)) {
            if (!internalOnly || child.isLeaf())
                logicalDescendents.addAll(getLogicalNodesInSubtree(child, internalOnly));
        }

        return logicalDescendents;
    }
}
