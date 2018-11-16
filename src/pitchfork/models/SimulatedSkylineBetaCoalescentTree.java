package pitchfork.models;

import beast.core.Input;
import beast.core.parameter.IntegerParameter;
import beast.core.parameter.RealParameter;
import beast.core.util.Log;
import beast.evolution.tree.Node;
import beast.evolution.tree.Tree;
import beast.math.distributions.ParametricDistribution;
import beast.util.Randomizer;

import java.io.FileNotFoundException;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Comparator;
import java.util.List;

/**
 * Simulator for Skyline-Beta coalescent trees.
 */
public class SimulatedSkylineBetaCoalescentTree extends Tree {

    public Input<BetaCoalescentModel> lcModelInput = new Input<>(
            "model",
            "Beta coalescent model.",
            Input.Validate.REQUIRED);

    public Input<RealParameter> initialPopSizeInput = new Input<>(
            "initialPopSize",
            "Initial population size.",
            new RealParameter("1.0"));

    public Input<IntegerParameter> intervalsPerGroupInput = new Input<>(
            "intervalsPerGroup",
            "Intervals per group in skyline plot.",
            Input.Validate.REQUIRED);

    public Input<String> fileNameInput = new Input<>(
            "fileName",
            "Name of file to save Newick representation of tree to.");

    public Input<RealParameter> populationSizeParameterInput = new Input<>(
            "popSizeParameter",
            "If provided, this parameter is initialized with the simulated " +
                    "population sizes.");

    public Input<RealParameter> changeTimesParameterInput = new Input<>(
            "changeTimesParameter",
            "If provided, this parameter is initialized with the simulated " +
                    "population change times.");

    private double[] leafAges;
    private String[] leafNames;
    private int nLeaves;

    private RealParameter initialPopSize;
    private IntegerParameter intervalsPerGroup;
    private BetaCoalescentModel lcModel;

    private List<Double> popSizes, changeTimes;

    public SimulatedSkylineBetaCoalescentTree() { }

    @Override
    public void initAndValidate() {
        super.initAndValidate();

        if (!hasDateTrait()) {
            throw new RuntimeException("Must specify tip dates!");
        }

        nLeaves = getDateTrait().taxaInput.get().getTaxonCount();
        leafAges = new double[nLeaves];
        leafNames = new String[nLeaves];

        initialPopSize = initialPopSizeInput.get();
        intervalsPerGroup = intervalsPerGroupInput.get();

        popSizes = new ArrayList<>();
        changeTimes = new ArrayList<>();

        popSizes.add(initialPopSize.getValue());
        changeTimes.add(0.0);

        for (int nodeNr=0; nodeNr<nLeaves; nodeNr++) {
            String taxonName = getDateTrait().taxaInput.get().getTaxonId(nodeNr);

            leafAges[nodeNr] = getDateTrait().getValue(taxonName);
            leafNames[nodeNr] = taxonName;
        }

        lcModel = lcModelInput.get();

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

    /**
     * Perform simulation.
     */
    private void simulate() {

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

        unusedLineages.sort(Comparator.comparingDouble(Node::getHeight));

        double t = 0.0;
        int interval = 0;
        double N = initialPopSize.getValue();
        while (!unusedLineages.isEmpty() || activeLineages.size() > 1) {

            // Compute propensities
            int n = activeLineages.size();
            double totalPropensity = n>=2 ? lcModel.getTotalCoalRate(n) : 0.0 ;

            // Increment time
            t += Randomizer.nextExponential(totalPropensity/N);

            // Check whether next sample time exceeded.
            if (!unusedLineages.isEmpty() && t>unusedLineages.get(0).getHeight()) {
                t = unusedLineages.get(0).getHeight();
                activeLineages.add(unusedLineages.remove(0));
                continue;
            }

            // Choose reaction
            double u = Randomizer.nextDouble()*totalPropensity;
            int k = 2 + -(Arrays.binarySearch(lcModel.getCumulativeCoalRateArray(n), u) + 1);

            if (k>n || k<2)
                throw new IllegalStateException("Numerical error in furcation degree sampler.");


            // Implement coalescence
            // (Note: BEAST really only deals with binary trees, so have to
            // fake multifurcations using zero-length edges.)

            Node newParent = new Node(String.valueOf(nextInternalNr));
            newParent.setNr(nextInternalNr++);
            newParent.setHeight(t);

            popSizes.add(N);
            changeTimes.add(t);

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

                popSizes.add(N);
                changeTimes.add(t);

                newNewParent.addChild(newParent);
                newNewParent.addChild(coalescingLineages.get(i));

                newParent = newNewParent;
            }

            activeLineages.removeAll(coalescingLineages);

            activeLineages.add(newParent);

            // Draw new population size if necessary

            interval += 1;

            if (interval % intervalsPerGroup.getValue() == 0)
                N = Randomizer.nextExponential(1.0/N);
        }

        assignFromWithoutID(new Tree(activeLineages.get(0)));

        if (populationSizeParameterInput.get() != null) {
            populationSizeParameterInput.get().assignFromWithoutID(
                    new RealParameter(popSizes.toArray(new Double[popSizes.size()])));
        }

        if (changeTimesParameterInput.get() != null) {
            changeTimesParameterInput.get().assignFromWithoutID(
                    new RealParameter(changeTimes.toArray(new Double[changeTimes.size()])));

        }
    }
}
