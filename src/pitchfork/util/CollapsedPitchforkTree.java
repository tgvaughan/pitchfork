package pitchfork.util;

import beast.core.Input;
import beast.evolution.tree.Node;
import beast.evolution.tree.Tree;

public class CollapsedPitchforkTree extends Tree {

    public Input<Tree> treeInput = new Input<>("tree",
            "Tree to collapse.",
            Input.Validate.REQUIRED);

    public CollapsedPitchforkTree() {};

    @Override
    public void initAndValidate() {
        setRoot(getCollapsedTree(treeInput.get().getRoot()));
        initArrays();
    }

    public Node getCollapsedTree(Node root) {
        Node newRoot = new Node();



        return newRoot;
    }
}
