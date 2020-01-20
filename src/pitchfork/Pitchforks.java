/*
 * Copyright (C) 2019. Tim Vaughan
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <https://www.gnu.org/licenses/>.
 */

package pitchfork;

import beast.evolution.tree.Node;
import beast.evolution.tree.Tree;
import beast.util.Randomizer;

import java.util.ArrayList;
import java.util.List;

/**
 * Class of static methods useful for traversing and manipulating
 * pitchfork's style of polytomy trees.
 */
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

    public static boolean isLogicalNode(Node node) {
        return node.isRoot() || node.getParent().getHeight()>node.getHeight();
    }

    public static Node getLogicalNode(Node node) {
        while (!isLogicalNode(node))
            node = node.getParent();

        return node;
    }

    public static Node getLogicalParent(Node logicalNode) {
        if (logicalNode.isRoot())
            return null;

        return getLogicalNode(logicalNode.getParent());
    }

    /**
     * Set the height of node in polytomy tree.  (May involve setting
     * the height of multiple real nodes in underlying BEAST tree.)
     *
     * @param node logical node or group member.
     * @param height new height of logical node and nodes in group.
     */
    public static void setLogicalNodeHeight(Node node, double height) {
        Node logicalNode = getLogicalNode(node);
        logicalNode.setHeight(height);

        for (Node gNode : getGroup(logicalNode))
            gNode.setHeight(height);
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

    /**
     * Get logical children descending from logical node groupRoot.
     *
     * @param groupRoot logical node in Tree
     * @return newly created list of child nodes.
     */
    public static List<Node> getLogicalChildren(Node groupRoot) {
        List<Node> logicalChildren = new ArrayList<>();

        getGroupAndLogicalChildren(groupRoot, null, logicalChildren);

        return logicalChildren;
    }

    /**
     * Returns true iff node is a member of the set of Tree nodes representing
     * a polytomy.
     *
     * @param node Tree node to determine polytomy status of
     * @return true iff node is polytomy.
     */
    public static boolean isPolytomy(Node node) {
        return (!node.isRoot() && node.getParent().getHeight() == node.getHeight())
                || (!node.isLeaf() && (node.getChildren().get(0).getHeight() == node.getHeight()
                || node.getChildren().get(1).getHeight() == node.getHeight()));
    }

    /**
     * Obtain a randomly chosen element of the given list.  (Something like
     * this should be in Randomizer.)
     *
     * @param list input list
     * @param <T> type of objects in list
     * @return object of type T.
     */
    public static <T> T randomChoice(List<T> list) {
        return list.get(Randomizer.nextInt(list.size()));
    }
}
