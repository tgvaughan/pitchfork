package pitchfork.util;

import beast.core.Description;
import beast.core.Input;
import beast.evolution.tree.Node;
import beast.evolution.tree.Tree;

import java.io.PrintStream;
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

    @Override
    public void log(int sample, PrintStream out) {
        update();
        out.print("tree STATE_" + sample + " = ");
        final int[] dummy = new int[1];
        final String newick = getRoot().toSortedNewick(dummy, false);
        out.print(newick);
        out.print(";");
    }

}
