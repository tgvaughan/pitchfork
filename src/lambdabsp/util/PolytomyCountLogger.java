package lambdabsp.util;

import beast.core.BEASTObject;
import beast.core.Input;
import beast.core.Loggable;
import beast.evolution.tree.Node;
import beast.evolution.tree.Tree;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.List;

public class PolytomyCountLogger extends BEASTObject implements Loggable {

    public Input<Tree> treeInput = new Input<>("tree",
            "Tree whose polytomies to count.",
            Input.Validate.REQUIRED);

    public PolytomyCountLogger() { }

    @Override
    public void initAndValidate() { }

    @Override
    public void init(PrintStream out) {
        if (getID() == null)
            out.print("PolytomyCount\t");
        else
            out.print(getID() + "\t");
    }

    public int getPolytomyCount() {
        int count = 0;

        List<Node> trueNodes = new ArrayList<>();
        for (Node node : treeInput.get().getNodesAsArray())
            if (node.isRoot() || node.getParent().getHeight()>node.getHeight())
                trueNodes.add(node);

        for (Node node : trueNodes) {
            if (!node.isLeaf() && (node.getChildren().get(0).getHeight()==node.getHeight()
                    || node.getChildren().get(1).getHeight()==node.getHeight()))
                count += 1;
        }

        return count;
    }

    @Override
    public void log(int sample, PrintStream out) {
        out.print(getPolytomyCount() + "\t");
    }

    @Override
    public void close(PrintStream out) {
    }
}
