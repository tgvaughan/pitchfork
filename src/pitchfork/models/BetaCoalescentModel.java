package pitchfork.models;

import beast.core.CalculationNode;
import beast.core.Input;
import beast.core.parameter.RealParameter;
import beast.evolution.alignment.TaxonSet;
import beast.evolution.tree.Tree;
import beast.math.Binomial;
import org.apache.commons.math.special.Beta;
import org.apache.commons.math.special.Gamma;

public class BetaCoalescentModel extends CalculationNode {

    public Input<RealParameter> alphaInput = new Input<>(
            "alpha",
            "Alpha parameter for Beta-coalescent process",
            Input.Validate.REQUIRED);

    public Input<TaxonSet> taxonSetInput = new Input<>(
            "taxonSet",
            "Taxon set used to define maximum number of extant lineages.");

    public Input<Tree> treeInput = new Input<>(
            "tree",
            "Tree used to define maximum number of extant lineages.",
            Input.Validate.XOR, taxonSetInput);

    private int nLeaves;
    private RealParameter alpha;
    private double logLambdaOffset;
    private double[][] logLambdaValues;
    private double[][] cumulativeCoalRates;
    private boolean[] dirtyFlags;

    private double storedLogLambdaOffset;
    private double[][] storedLogLambdaValues;
    private double[][] storedCumulativeCoalRates;
    private boolean[] storedDirtyFlags;


    @Override
    public void initAndValidate() {
        if (taxonSetInput.get() != null)
            nLeaves = taxonSetInput.get().getTaxonCount();
        else
            nLeaves = treeInput.get().getLeafNodeCount();

        alpha = alphaInput.get();

        logLambdaValues = new double[nLeaves-1][];
        storedLogLambdaValues = new double[nLeaves-1][];
        cumulativeCoalRates = new double[nLeaves-1][];
        storedCumulativeCoalRates = new double[nLeaves-1][];
        for (int n=2; n<=nLeaves; n++) {
            logLambdaValues[n - 2] = new double[n - 1];
            storedLogLambdaValues[n - 2] = new double[n - 1];
            cumulativeCoalRates[n - 2] = new double[n - 1];
            storedCumulativeCoalRates[n - 2] = new double[n - 1];
        }

        dirtyFlags = new boolean[nLeaves-1];
        storedDirtyFlags = new boolean[nLeaves-1];
        for (int n=2; n<=nLeaves; n++) {
            dirtyFlags[n-2] = true;
            storedDirtyFlags[n-2] = true;
        }
    }

    private void makeDirty(int n) {
        dirtyFlags[n-2] = true;
    }

    private void makeClean(int n) {
        dirtyFlags[n-2] = false;
    }

    private void makeAllDirty() {
        for (int n=2; n<=nLeaves; n++)
            makeDirty(n);
    }

    private boolean isDirty(int n) {
        return n>=2 && dirtyFlags[n-2];
    }

    private void computeCoalRateDistribs(int n) {
        logLambdaOffset = -Beta.logBeta(2-alpha.getValue(), alpha.getValue());

        logLambdaValues[n-2][0] = logLambdaOffset + Beta.logBeta(2-alpha.getValue(), n-2+alpha.getValue());
        cumulativeCoalRates[n-2][0] = Math.exp(logLambdaValues[n-2][0] + Binomial.logChoose(n, 2));

        for (int k=3; k<=n; k++) {
            logLambdaValues[n-2][k-2] = logLambdaOffset + Beta.logBeta(k-alpha.getValue(), n-k+alpha.getValue());
            cumulativeCoalRates[n-2][k-2] = cumulativeCoalRates[n-2][k-3]
                    + Math.exp(logLambdaValues[n-2][k-2] + Binomial.logChoose(n, k));
        }
    }

    private void update(int n) {
        if (!isDirty(n))
            return;

        computeCoalRateDistribs(n);

        makeClean(n);
    }

    public double getLogLambda(int n, int k) {
        update(n);

        if (n<2)
            return 0;
        else
            return logLambdaValues[n-2][k-2];
    }

    public double getTotalCoalRate(int n) {
        update(n);

        if (n<2)
            return 0;
        else
            return cumulativeCoalRates[n-2][n-2];
    }

    public double[] getCumulativeCoalRateArray(int n) {
        update(n);

        return cumulativeCoalRates[n-2];
    }

    @Override
    protected void store() {

        storedLogLambdaOffset = logLambdaOffset;

        for (int n=2; n<=nLeaves; n++) {
            for (int k=2; k<=n; k++) {
                storedLogLambdaValues[n-2][k-2] = logLambdaValues[n-2][k-2];
                storedCumulativeCoalRates[n-2][k-2] = cumulativeCoalRates[n-2][k-2];
            }
            storedDirtyFlags[n-2] = dirtyFlags[n-2];
        }

        super.store();
    }

    @Override
    protected void restore() {

        logLambdaOffset = storedLogLambdaOffset;

        double [][] tmp = logLambdaValues;
        logLambdaValues = storedLogLambdaValues;
        storedLogLambdaValues = tmp;

        tmp = cumulativeCoalRates;
        cumulativeCoalRates = storedCumulativeCoalRates;
        storedCumulativeCoalRates = tmp;

        boolean [] tmpFlags = dirtyFlags;
        dirtyFlags = storedDirtyFlags;
        storedDirtyFlags = tmpFlags;

        super.restore();
    }

    @Override
    protected boolean requiresRecalculation() {
        makeAllDirty();
        return true;
    }
}
