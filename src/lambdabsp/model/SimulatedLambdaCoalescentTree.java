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

    public SimulatedLambdaCoalescentTree() { }

    @Override
    public void initAndValidate() {
        super.initAndValidate();

        if (!hasDateTrait()) {
            throw new RuntimeException("Must specify tip dates!");
        }

        nLeaves = getDateTrait().taxaInput.get().getTaxonCount();
        leafAges = new double[nLeaves];
        leafNames = new String[nLeaves];


        for (int nodeNr=0; nodeNr<nLeaves; nodeNr++) {
            String taxonName = getDateTrait().taxaInput.get().getTaxonId(nodeNr);

            leafAges[nodeNr] = getDateTrait().getValue(taxonName);
            leafNames[nodeNr] = taxonName;
        }

        populationFunction = populationFunctionInput.get();

        initArrays();
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

        return Beta.logBeta(k - alpha, n - k + alpha)
                - Beta.logBeta(2-alpha, alpha);
    }

    private double[][] cumulativeCoalRates;
    private void computeCoalRateDistribs() {
        cumulativeCoalRates = new double[nLeaves-1][nLeaves-1];

        for (int n=2; n<=nLeaves; n++) {
            cumulativeCoalRates[n-2][0] = Math.exp(getLogLambda(n, 2) + Binomial.logChoose(n, 2));

            for (int k=3; k<=n; k++) {
                cumulativeCoalRates[n-2][k-2] = cumulativeCoalRates[n-2][k-3]
                        + Math.exp(getLogLambda(n, k) + Binomial.logChoose(n, k));
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
            double totalPropensity = n>=2 ? cumulativeCoalRates[n-2][n-2] : 0.0 ;

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
            int k;
            for (k=2; k<=n; k++) {
                u -= cumulativeCoalRates[n-2][k-2];

                if (u < 0)
                    break;
            }

            if (k>n)
                throw new IllegalStateException("Numerical error: loop in coalescent simulator fell through.");

            // Implement coalescence
            // (Note: BEAST really only deals with binary trees, so have to
            // fake multifurcations using zero-length edges.)

            Node newParent = new Node(String.valueOf(nextInternalNr));
            newParent.setNr(nextInternalNr++);
            newParent.setHeight(t);

            int[] indices = Randomizer.shuffled(n);

            coalescingLineages.clear();
            for (int i=0; i<k; i++)
                coalescingLineages.add(activeLineages.get(indices[i]));

            newParent.addChild(coalescingLineages.get(0));
            newParent.addChild(coalescingLineages.get(1));

            for (int i=2; i<k; i++) {
                Node newNewParent = new Node(String.valueOf(nextInternalNr));
                newNewParent.setNr(nextInternalNr++);
                newNewParent.setHeight(t);

                newNewParent.addChild(newParent);
                newNewParent.addChild(coalescingLineages.get(i));

                newParent = newNewParent;
            }

            activeLineages.removeAll(coalescingLineages);

            activeLineages.add(newParent);
        }

//        setRoot(activeLineages.get(0));
        assignFromWithoutID(new Tree(activeLineages.get(0)));
    }
}
