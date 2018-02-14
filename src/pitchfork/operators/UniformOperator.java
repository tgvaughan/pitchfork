package pitchfork.operators;

import beast.core.Description;
import beast.core.Input;
import beast.evolution.tree.Node;
import beast.evolution.tree.Tree;
import beast.util.Randomizer;

import java.util.ArrayList;
import java.util.List;

@Description("Uniform node height operator compatible with trees having polytomies.")
public class UniformOperator extends PitchforkTreeOperator {

    public Input<Boolean> scaleRootInput = new Input<>(
            "scaleRoot",
            "Whether to scale the age of the root node.",
            true);

    public Input<Double> scaleFactorInput = new Input<>(
            "scaleFactor",
            "Tuning parameter for scaling root.",
            0.8);

    Tree tree;

    boolean scaleRoot;
    double scaleFactor;

    @Override
    public void initAndValidate() {
        tree = treeInput.get();
        scaleRoot = scaleRootInput.get();
        scaleFactor = scaleFactorInput.get();
    }

    @Override
    public double proposal() {

        double logHR = 0.0;

        List<Node> trueNodes = getTrueInternalNodes();

        if (trueNodes.size() == 1 && !scaleRoot)
            return Double.NEGATIVE_INFINITY;

        Node logicalNode;
        do {
            logicalNode = trueNodes.get(Randomizer.nextInt(trueNodes.size()));
        } while (!scaleRoot && logicalNode.isRoot());


        List<Node> nodesInLogicalGroup = new ArrayList<>();
        List<Node> logicalChildren = new ArrayList<>();
        getGroupAndLogicalChildren(logicalNode, nodesInLogicalGroup, logicalChildren);
        double minHeight = logicalChildren.stream().mapToDouble(Node::getHeight).max().getAsDouble();

        double newHeight;
        if (logicalNode.isRoot()) {

            double minf = Math.min(scaleFactor, 1.0/scaleFactor);
            double maxf = 1.0/minf;
            double f = minf + Randomizer.nextDouble()*(maxf - minf);

            newHeight = logicalNode.getHeight() * f;

            if (newHeight < minHeight)
                return Double.NEGATIVE_INFINITY;

            logHR = -Math.log(f);

        } else {
            Node parent = logicalNode.getParent();
            double maxHeight = parent.getHeight();

            newHeight = minHeight + (maxHeight - minHeight) * Randomizer.nextDouble();
        }

        logicalNode.setHeight(newHeight);
        for (Node node : nodesInLogicalGroup)
            node.setHeight(newHeight);

        return logHR;
    }
}
