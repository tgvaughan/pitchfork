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

package pitchfork.operators;

import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.evolution.operator.TreeOperator;
import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.Tree;
import beast.base.util.Randomizer;
import pitchfork.Pitchforks;

import java.util.List;

@Description("Exchange operator compatible with pitchfork trees.")
public class ExchangeOperator extends TreeOperator {

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

    @Override
    public double proposal() {

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
