package lambdabsp.operators;

import beast.core.Input;
import beast.evolution.operators.SubtreeSlide;
import beast.evolution.operators.TreeOperator;
import beast.evolution.tree.Node;
import beast.evolution.tree.Tree;
import beast.util.Randomizer;
import beast.util.TreeParser;
import lambdabsp.model.CollapsedTreeIntervals;

import java.lang.reflect.Array;
import java.util.ArrayList;
import java.util.List;
import java.util.Random;

/**
 * Class of operators for traversing the space of multifurcating phylogenetic trees.
 */
public abstract class LambdaTreeOperator extends TreeOperator {

    public List<Node> getTrueNodes() {
        List<Node> trueNodes = new ArrayList<>();

        for (Node node : treeInput.get().getNodesAsArray())
            if (node.isRoot() || node.getParent().getHeight() > node.getHeight())
                trueNodes.add(node);

        return trueNodes;
    }

    public Node getLogicalNode(Node node) {
        while (!node.isRoot() && node.getParent().getHeight()==node.getHeight())
            node = node.getParent();

        return node;
    }

    public boolean isPolytomy(Node node) {
        return (!node.isRoot() && node.getParent().getHeight() == node.getHeight())
                || (!node.isLeaf() && (node.getChildren().get(0).getHeight() == node.getHeight()
                || node.getChildren().get(1).getHeight() == node.getHeight()));

    }

}
