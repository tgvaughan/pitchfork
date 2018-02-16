package pitchfork.operators;

import beast.core.Description;
import beast.core.Input;
import beast.evolution.tree.Node;
import beast.evolution.tree.Tree;
import beast.util.Randomizer;
import pitchfork.util.Pitchforks;

import java.util.List;

@Description("Exchange operator compatible with pitchfork trees.")
public class ExchangeOperator extends PitchforkTreeOperator {

    public Input<Boolean> isNarrowInput = new Input<>(
            "isNarrow",
            "Whether narrow exchange is used.",
            true);

    boolean isNarrow;
    Tree tree;

    @Override
    public void initAndValidate() {
        isNarrow = isNarrowInput.get();
        tree = treeInput.get();
    }

//    int count = 0;

    @Override
    public double proposal() {

//        System.out.println("Count: " + (++count));

        if (isNarrow) {
            List<Node> trueNodes = Pitchforks.getTrueNodes(tree);

            if (trueNodes.size() - tree.getLeafNodeCount() <= 1)
                return Double.NEGATIVE_INFINITY;

            Node srcNode, srcNodeParent, destNode, destNodeParent;
            Node srcNodeLogicalParent, srcNodeLogicalGrandparent = null;

            do {
                srcNode = trueNodes.get(Randomizer.nextInt(trueNodes.size()));
                srcNodeLogicalParent = Pitchforks.getLogicalParent(srcNode);

                if (srcNodeLogicalParent != null)
                    srcNodeLogicalGrandparent = Pitchforks.getLogicalParent(srcNodeLogicalParent);

            } while (srcNodeLogicalParent == null || srcNodeLogicalGrandparent == null);
            srcNodeParent = srcNode.getParent();

            List<Node> possibleDestNodes = Pitchforks.getLogicalChildren(srcNodeLogicalGrandparent);

            do {
                destNode = possibleDestNodes.get(Randomizer.nextInt(possibleDestNodes.size()));
            } while (destNode == srcNodeLogicalParent);
            destNodeParent = destNode.getParent();

            // Reject if substitution would result in negative edge length:
            if (destNode.getHeight() > srcNodeParent.getHeight()
                || srcNode.getHeight() > destNodeParent.getHeight())
                return Double.NEGATIVE_INFINITY;

            srcNodeParent.removeChild(srcNode);
            destNodeParent.removeChild(destNode);
            srcNodeParent.addChild(destNode);
            destNodeParent.addChild(srcNode);

        } else {
            throw new UnsupportedOperationException("Wide exchange for " +
                    "pitchfork trees is not yet supported.");
        }


        return 0;
    }
}
