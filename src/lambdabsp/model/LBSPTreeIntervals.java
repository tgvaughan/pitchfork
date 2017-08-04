package lambdabsp.model;

import beast.core.CalculationNode;
import beast.core.Input;
import beast.evolution.tree.Node;
import beast.evolution.tree.Tree;
import beast.evolution.tree.coalescent.IntervalList;
import javafx.collections.transformation.SortedList;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.TreeSet;

/**
 * Computes tree intervals, skipping zero-length intervals.
 */
public class LBSPTreeIntervals extends CalculationNode implements IntervalList {

    public Input<Tree> treeInput = new Input<>("tree",
            "Tree from which to compute coalescent intervals.",
            Input.Validate.REQUIRED);

    boolean isDirty;

    List<Integer> lineageCounts;
    List<Double> intervalDurations;
    int nIntervals, nSamples;
    boolean isBinaryTree;

    Tree tree;

    TreeSet<Node> sortedNodeSet = new TreeSet<>((n1, n2) -> {
        if (n1.getHeight() < n2.getHeight())
            return -1;

        if (n2.getHeight() < n1.getHeight())
            return 1;

        return 0;
    });

    @Override
    public void initAndValidate() {

        tree = treeInput.get();

        lineageCounts = new ArrayList<>();
        intervalDurations = new ArrayList<>();
        isDirty = true;
    }

    void update() {
        if (!isDirty)
            return;

        // Do update

        sortedNodeSet.addAll(Arrays.asList(tree.getNodesAsArray()));
        int lineages = 0;
        int intervalIndex = 0;
        double prevHeight = 0;
        for (Node node : sortedNodeSet) {

            if (node.isLeaf())
                lineages += 1;
            else
                lineages -= 1;


        }

        isDirty = false;
    }

    @Override
    public int getIntervalCount() {
        update();



        return nIntervals;
    }

    @Override
    public int getSampleCount() {
        return nSamples;
    }

    @Override
    public double getInterval(int i) {
        return intervalDurations.get(i);
    }

    @Override
    public int getLineageCount(int i) {
        return lineageCounts.get(i);
    }

    @Override
    public int getCoalescentEvents(int i) {
        return 0;
    }

    @Override
    public IntervalType getIntervalType(int i) {
        return null;
    }

    @Override
    public double getTotalDuration() {
        return tree.getRoot().getHeight();
    }

    @Override
    public boolean isBinaryCoalescent() {
        return isBinaryTree;
    }

    @Override
    public boolean isCoalescentOnly() {
        return nSamples == 0;
    }

    @Override
    protected boolean requiresRecalculation() {
        isDirty = true;

        return true;
    }
}
