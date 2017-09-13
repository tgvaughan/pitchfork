package lambdabsp.model;

import beast.core.Input;
import beast.core.Logger;
import beast.core.parameter.RealParameter;
import beast.core.util.Log;
import beast.evolution.tree.Node;
import beast.evolution.tree.Tree;
import beast.evolution.tree.coalescent.PopulationFunction;
import beast.math.Binomial;
import beast.util.Randomizer;
import org.apache.commons.math.special.Beta;

import java.io.FileNotFoundException;
import java.io.IOError;
import java.io.PrintStream;
import java.util.*;

/**
 * Simulator for Lambda coalescent trees.
 */
public class SimulatedLambdaCoalescentTree extends Tree {

    public Input<RealParameter> alphaInput = new Input<>(
            "alpha",
            "Parameter determining shape of Lambda distribution",
            Input.Validate.REQUIRED);

    public Input<PopulationFunction> populationFunctionInput = new Input<>(
            "populationFunction",
            "Population function used in coalescent simulation.",
            Input.Validate.REQUIRED);

    public Input<String> fileNameInput = new Input<>(
            "fileName",
            "Name of file to save Newick representation of tree to.");

    private double[] leafAges;
    private String[] leafNames;
    private int nLeaves;

    private PopulationFunction populationFunction;

    @Override
    public void initAndValidate() {
        super.initAndValidate();

        nLeaves = getDateTrait().taxaInput.get().getTaxonCount();
        leafAges = new double[nLeaves];
        leafNames = new String[nLeaves];

        for (int nodeNr=0; nodeNr<nLeaves; nodeNr++) {
            String taxonName = getDateTrait().taxaInput.get().getTaxonId(nodeNr);

            leafAges[nodeNr] = getDateTrait().getValue(taxonName);
            leafNames[nodeNr] = taxonName;
        }

        populationFunction = populationFunctionInput.get();

        simulate();

        // Write output file
        if (fileNameInput.get() != null) {
            try (PrintStream ps = new PrintStream(fileNameInput.get())) {
                ps.println(getRoot().toNewick().concat(";"));
            } catch (FileNotFoundException ex) {
                Log.err.println("Could not write to output file.");
                System.exit(1);
            }
        }
    }

    private double getLogLambda(int n, int k) {

        double alpha = alphaInput.get().getValue();

        return Binomial.logChoose(nLeaves, k)
                + Beta.logBeta(k - alpha, nLeaves - k + alpha)
                - Beta.logBeta(2-alpha, alpha);
    }

    private double[][] cumulativeCoalRates;
    private void computeCoalRateDistribs() {
        cumulativeCoalRates = new double[nLeaves-1][nLeaves-1];

        for (int n=2; n<=nLeaves; n++) {
            cumulativeCoalRates[n][0] = Math.exp(getLogLambda(n, 2) + Binomial.choose(n, 2));

            for (int k=3; k<=n; k++) {
                cumulativeCoalRates[n-2][k-2] = cumulativeCoalRates[n][k-3]
                        + Math.exp(getLogLambda(n, k) + Binomial.choose(n, k));
            }
        }
    }

    /**
     * Perform simulation.
     */
    public void simulate() {

        computeCoalRateDistribs();

        List<Node> activeLineages = new ArrayList<>();
        List<Node> unusedLineages = new ArrayList<>();
        List<Node> coalescingLineages = new ArrayList<>();

        // Create leaf nodes:
        for (int i=0; i<nLeaves; i++)  {
            Node leaf = new Node(leafNames[i]);
            leaf.setNr(i);
            leaf.setHeight(leafAges[i]);
            unusedLineages.add(leaf);
        }

        int nextInternalNr = nLeaves;

        // TODO Check sort direction!
        unusedLineages.sort(Comparator.comparingDouble(Node::getHeight));

        double t = 0;
        while (!unusedLineages.isEmpty() || activeLineages.size() > 1) {

            // Compute propensities

            int n = activeLineages.size();
            double totalPropensity = cumulativeCoalRates[n-2][n-2];

            // Increment time
            t += Randomizer.nextExponential(totalPropensity);

            // Check whether next sample time exceeded.
            if (!unusedLineages.isEmpty() && t>unusedLineages.get(0).getHeight()) {
                t = unusedLineages.get(0).getHeight();
                activeLineages.add(unusedLineages.remove(0));
                continue;
            }

            // Choose reaction
            double u = Randomizer.nextDouble()*totalPropensity;
            int k = -1;
            for (k=2; k<=n; k++) {
                u -= cumulativeCoalRates[n-2][k-2];

                if (u < 0)
                    break;
            }

            if (k<0)
                throw new IllegalStateException("Numerical error: loop in coalescent simulator fell through.");

            // Implement coalescence
            Node newParent = new Node();
            newParent.setNr(nextInternalNr++);
            newParent.setHeight(t);

            int[] indices = Randomizer.shuffled(n);

            for (int i=0; i<k; i++)
                newParent.addChild(activeLineages.get(indices[i]));

            for (Node child : newParent.getChildren())
                activeLineages.remove(child);

            activeLineages.add(newParent);
        }

        assignFromWithoutID(new Tree(activeLineages.get(0)));
    }
}
