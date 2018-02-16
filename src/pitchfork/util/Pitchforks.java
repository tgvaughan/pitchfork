package pitchfork.util;

import beast.evolution.tree.Node;
import beast.evolution.tree.Tree;

import java.util.ArrayList;
import java.util.List;

public class Pitchforks {

    public static List<Node> getTrueNodes(Tree tree) {
        List<Node> trueNodes = new ArrayList<>(getTrueInternalNodes(tree));

        for (int nodeNr=0; nodeNr < tree.getLeafNodeCount(); nodeNr++) {
            trueNodes.add(tree.getNode(nodeNr));
        }

        return trueNodes;
    }

    public static List<Node> getTrueInternalNodes(Tree tree) {
        List<Node> trueNodes = new ArrayList<>();

        for (int nodeNr = tree.getLeafNodeCount();
             nodeNr < tree.getNodeCount();
             nodeNr += 1) {
            Node node = tree.getNode(nodeNr);

            if (node.isRoot() || node.getParent().getHeight() > node.getHeight())
                trueNodes.add(node);
        }

        return trueNodes;
    }

    public static Node getLogicalNode(Node node) {
        while (!node.isRoot() && node.getParent().getHeight()==node.getHeight())
            node = node.getParent();

        return node;
    }

    public static Node getLogicalParent(Node logicalNode) {
        if (logicalNode.isRoot())
            return null;

        return getLogicalNode(logicalNode.getParent());
    }

    public static void getGroupAndLogicalChildren(Node node, List<Node> group, List<Node> logicalChildren) {
        for (Node child : node.getChildren()) {
            if (child.getHeight() == node.getHeight()) {
                if (group != null)
                    group.add(child);
                getGroupAndLogicalChildren(child, group, logicalChildren);
            } else {
                if (logicalChildren != null)
                    logicalChildren.add(child);
            }
        }
    }

    /**
     * Get nodes included in same logical node as groupRoot which are descendants
     * of groupRoot.  (Excludes groupRoot itself.)
     * @param groupRoot root of logical node group.
     * @return list of nodes comprising group which descend from groupRoot.
     */
    public static List<Node> getGroup(Node groupRoot) {
        List<Node> group = new ArrayList<>();

        getGroupAndLogicalChildren(groupRoot, group, null);

        return group;
    }

    public static List<Node> getLogicalChildren(Node groupRoot) {
        List<Node> logicalChildren = new ArrayList<>();

        getGroupAndLogicalChildren(groupRoot, null, logicalChildren);

        return logicalChildren;
    }
    public static boolean isPolytomy(Node node) {
        return (!node.isRoot() && node.getParent().getHeight() == node.getHeight())
                || (!node.isLeaf() && (node.getChildren().get(0).getHeight() == node.getHeight()
                || node.getChildren().get(1).getHeight() == node.getHeight()));
    }
}
