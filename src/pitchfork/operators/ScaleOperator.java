package pitchfork.operators;

import beast.core.Input;
import beast.core.parameter.RealParameter;
import beast.evolution.tree.Node;
import beast.evolution.tree.Tree;
import beast.util.Randomizer;

import java.util.ArrayList;
import java.util.List;

public class ScaleOperator extends PitchforkTreeOperator {

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


    @Override
    public void initAndValidate() {
        tree = treeInput.get();
    }

    @Override
    public double proposal() {

        // Choose scale factor:
        double minf = Math.min(scaleFactorInput.get(), 1.0/scaleFactorInput.get());
        double maxf = 1.0/minf;
        double f = minf + (maxf - minf)* Randomizer.nextDouble();

        // Keep track of Hastings factor:
        double logf = Math.log(f);
        double logHR = -2*logf;

        // Scale tree
        try {
            tree.getRoot().scale(f);
        } catch (IllegalArgumentException ex) {
            return Double.NEGATIVE_INFINITY;
        }
        for (int nodeNr = tree.getLeafNodeCount(); nodeNr<tree.getNodeCount(); nodeNr++) {
            Node node = tree.getNode(nodeNr);
            if (node.isRoot() || node.getHeight()<node.getParent().getHeight())
                logHR += logf;
        }

        // Scale parameters
        for (RealParameter param : parametersInput.get()) {
            try {
                int nScaled = param.scale(f);
                logHR += param.getDimension()*nScaled;
            } catch (IllegalArgumentException ex) {

                // Parameter scaled out of range
                return Double.NEGATIVE_INFINITY;
            }
        }

        // Scale inverse parameters
        for (RealParameter param : parametersInverseInput.get()) {
            try {
                int nScaled = param.scale(1.0/f);
                logHR -= param.getDimension()*nScaled;
            } catch (IllegalArgumentException ex) {

                // Parameter scaled out of range
                return Double.NEGATIVE_INFINITY;
            }
        }

        return logHR;
    }

}
