package lambdabsp.operators;

import beast.core.Description;
import beast.evolution.tree.Node;
import beast.evolution.tree.Tree;
import beast.util.Randomizer;

import java.util.ArrayList;
import java.util.List;

@Description("Uniform node height operator compatible with trees having polytomies.")
public class UniformOperator extends LambdaTreeOperator {

    Tree tree;

    @Override
    public void initAndValidate() {
        tree = treeInput.get();
    }

    @Override
    public double proposal() {

        List<Node> trueNodes = getTrueInternalNodes();

        Node logicalNode;
        do {
            logicalNode = trueNodes.get(Randomizer.nextInt(trueNodes.size()));
        } while (logicalNode.isRoot());


        List<Node> nodesInLogicalGroup = new ArrayList<>();
        List<Node> logicalChildren = new ArrayList<>();
        getGroupAndLogicalChildren(logicalNode, nodesInLogicalGroup, logicalChildren);

        Node parent = logicalNode.getParent();

        double maxHeight = parent.getHeight();
        double minHeight = logicalChildren.stream().mapToDouble(Node::getHeight).max().getAsDouble();

        double newHeight = minHeight + (maxHeight-minHeight)*Randomizer.nextDouble();

        logicalNode.setHeight(newHeight);
        for (Node node : nodesInLogicalGroup)
            node.setHeight(newHeight);

        return 0;
    }
}
