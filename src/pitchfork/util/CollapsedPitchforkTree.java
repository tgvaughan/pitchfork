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

package pitchfork.util;

import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.Tree;
import pitchfork.Pitchforks;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

@Description("Tree with true polytomies for logging only.")
public class CollapsedPitchforkTree extends Tree {

    public Input<Tree> treeInput = new Input<>("tree",
            "Tree to collapse.",
            Input.Validate.REQUIRED);

    public CollapsedPitchforkTree() {};

    @Override
    public void initAndValidate() {
        update();
    }

    private int nextNodeNr;

    public void update() {
        nextNodeNr = treeInput.get().getLeafNodeCount();
        setRoot(getCollapsedTree(treeInput.get().getRoot()));

        // This breaks for some reason to do with varying node numbers.
        // initArrays();
        // (Not necessary for logging, so leaving out for now.)
    }

    public Node getCollapsedTree(Node root) {
        Node newRoot = new Node();
        newRoot.setHeight(root.getHeight());

        for (String key : root.getMetaDataNames())
            newRoot.setMetaData(key, root.getMetaData(key));

        newRoot.metaDataString = root.metaDataString;

        newRoot.setID(root.getID());

        if (root.isLeaf()) {
            newRoot.setNr(root.getNr());
        } else {
            newRoot.setNr(nextNodeNr++);
        }

        List<Node> logicalChildren = Pitchforks.getLogicalChildren(root);

        for (Node child : logicalChildren)
            newRoot.addChild(getCollapsedTree(child));

        return newRoot;
    }

    /* Need to override these methods because the Tree implementation
       assumes a binary tree. */

    /**
     * print translate block for NEXUS beast.tree file
     */
    public static void printTranslate(final Node node, final PrintStream out, final int nodeCount) {
        final List<String> translateLines = new ArrayList<>();
        printTranslate(node, translateLines, nodeCount);
        Collections.sort(translateLines);
        for (final String line : translateLines) {
            out.println(line);
        }
    }

    private static void printTranslate(Node node, List<String> translateLines, int nodeCount) {
        if (node.isLeaf()) {
            final String nr = (node.getNr() + taxaTranslationOffset) + "";
            String line = "\t\t" + "    ".substring(nr.length()) + nr + " " + node.getID();
            if (node.getNr() < nodeCount - 1) {
                line += ",";
            }
            translateLines.add(line);
        } else {
            for (Node child : node.getChildren())
                printTranslate(child, translateLines, nodeCount);
        }
    }

    public static void printTaxa(final Node node, final PrintStream out, final int nodeCount) {
        final List<String> translateLines = new ArrayList<>();
        printTranslate(node, translateLines, nodeCount);
        Collections.sort(translateLines);
        for (String line : translateLines) {
            line = line.split("\\s+")[2];
            out.println("\t\t\t" + line.replace(',', ' '));
        }
    }


    @Override
	public void init(PrintStream out) {
        update();

        Node node = getRoot();
        out.println("#NEXUS\n");
        out.println("Begin taxa;");
        out.println("\tDimensions ntax=" + root.getLeafNodeCount() + ";");
        out.println("\t\tTaxlabels");
        printTaxa(node, out, root.getLeafNodeCount());
        out.println("\t\t\t;");
        out.println("End;");

        out.println("Begin trees;");
        out.println("\tTranslate");
        printTranslate(node, out, root.getLeafNodeCount());
        out.print(";");
    }


    @Override
    public void log(long sample, PrintStream out) {
        update();
        out.print("tree STATE_" + sample + " = ");
        final int[] dummy = new int[1];
        final String newick = getRoot().toSortedNewick(dummy, false);
        out.print(newick);
        out.print(";");
    }

}
