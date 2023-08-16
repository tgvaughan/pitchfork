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
import beast.base.inference.parameter.RealParameter;
import beast.base.util.Randomizer;
import pitchfork.Pitchforks;

import java.util.ArrayList;
import java.util.List;

/**
 * This scale operator is largely redundant now that the beast core ScaleOperator
 * correctly computes the HR when zero-length internal edges exist.  However,
 * it is still useful as it can correctly perform root-only scale proposals
 * when the root is a polytomy, while the core operator cannot.
 */
@Description("Scale operator for pitchfork trees.")
public class ScaleOperator extends TreeOperator {

    public Input<Boolean> rootOnlyInput = new Input<>("rootOnly",
            "Scale only age of root node node.",
            false);

    public Input<List<RealParameter>> parametersInput =
            new Input<>("parameter",
                    "Scale this scalar parameter by the same amount as tree.",
                    new ArrayList<RealParameter>());

    public Input<List<RealParameter>> parametersInverseInput =
            new Input<>("parameterInverse",
                    "Scale this scalar parameter inversely.",
                    new ArrayList<RealParameter>());


    public Input<Double> scaleFactorInput = new Input<>("scaleFactor",
            "Scaling is restricted to the range [1/scaleFactor, scaleFactor]",
            0.8);

    private Tree tree;
    private boolean rootOnly;


    @Override
    public void initAndValidate() {
        tree = treeInput.get();
        rootOnly = rootOnlyInput.get();
    }

    @Override
    public double proposal() {

        if (tree.getRoot().isLeaf())
            return Double.NEGATIVE_INFINITY;

        // Choose scale factor:
        double minf = Math.min(scaleFactorInput.get(), 1.0/scaleFactorInput.get());
        double maxf = 1.0/minf;
        double f = minf + (maxf - minf)* Randomizer.nextDouble();

        // Keep track of Hastings factor:
        double logf = Math.log(f);
        double logHR = -2*logf;

        // Scale tree
        if (rootOnly) {
            List<Node> rootGroup = new ArrayList<>();
            List<Node> logicalChildren = new ArrayList<>();
            Pitchforks.getGroupAndLogicalChildren(tree.getRoot(), rootGroup, logicalChildren);
            rootGroup.add(tree.getRoot());

            double newHeight = f*tree.getRoot().getHeight();

            if (f<1.0) {
                for (Node child : logicalChildren)
                    if (newHeight < child.getHeight())
                        return Double.NEGATIVE_INFINITY;
            }

            for (Node node : rootGroup)
                node.setHeight(newHeight);

            logHR += logf;

        } else {
            try {
                tree.getRoot().scale(f);
            } catch (IllegalArgumentException ex) {
                return Double.NEGATIVE_INFINITY;
            }
            for (int nodeNr = tree.getLeafNodeCount(); nodeNr < tree.getNodeCount(); nodeNr++) {
                Node node = tree.getNode(nodeNr);
                if (node.isRoot() || node.getHeight() < node.getParent().getHeight())
                    logHR += logf;
            }
        }

        // Scale parameters
        for (RealParameter param : parametersInput.get()) {
            try {
                param.startEditing(this);
                int nScaled = param.scale(f);
                logHR += nScaled*logf;
            } catch (IllegalArgumentException ex) {

                // Parameter scaled out of range
                return Double.NEGATIVE_INFINITY;
            }
        }

        // Scale inverse parameters
        for (RealParameter param : parametersInverseInput.get()) {
            try {
                param.startEditing(this);
                int nScaled = param.scale(1.0/f);
                logHR -= nScaled*logf;
            } catch (IllegalArgumentException ex) {

                // Parameter scaled out of range
                return Double.NEGATIVE_INFINITY;
            }
        }

        return logHR;
    }

}
