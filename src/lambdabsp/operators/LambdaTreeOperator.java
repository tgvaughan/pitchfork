package lambdabsp.operators;

import beast.evolution.operators.TreeOperator;
import beast.evolution.tree.Node;

import java.util.ArrayList;
import java.util.List;

/**
 * Class of operators for traversing the space of multifurcating phylogenetic trees.
 */
public abstract class LambdaTreeOperator extends TreeOperator {

    public List<Node> getTrueNodes() {
        List<Node> trueNodes = new ArrayList<>(getTrueInternalNodes());

        for (int nodeNr=0; nodeNr < treeInput.get().getLeafNodeCount(); nodeNr++) {
            trueNodes.add(treeInput.get().getNode(nodeNr));
        }

        return trueNodes;
    }

    public List<Node> getTrueInternalNodes() {
         List<Node> trueNodes = new ArrayList<>();

         for (int nodeNr = treeInput.get().getLeafNodeCount();
                 nodeNr < treeInput.get().getNodeCount();
                 nodeNr += 1) {
             Node node = treeInput.get().getNode(nodeNr);

             if (node.isRoot() || node.getParent().getHeight() > node.getHeight())
                 trueNodes.add(node);
         }

        return trueNodes;
    }

    public Node getLogicalNode(Node node) {
        while (!node.isRoot() && node.getParent().getHeight()==node.getHeight())
            node = node.getParent();

        return node;
    }

    public void getGroupAndLogicalChildren(Node node, List<Node> group, List<Node> logicalChildren) {
        for (Node child : node.getChildren()) {
            if (child.getHeight() == node.getHeight()) {
                group.add(child);
                getGroupAndLogicalChildren(child, group, logicalChildren);
            } else {
                logicalChildren.add(child);
            }
        }
    }

    public boolean isPolytomy(Node node) {
        return (!node.isRoot() && node.getParent().getHeight() == node.getHeight())
                || (!node.isLeaf() && (node.getChildren().get(0).getHeight() == node.getHeight()
                || node.getChildren().get(1).getHeight() == node.getHeight()));

    }

}
