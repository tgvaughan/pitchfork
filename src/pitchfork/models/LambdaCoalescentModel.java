package pitchfork.models;

import beast.core.CalculationNode;
import beast.core.Input;
import beast.core.parameter.RealParameter;
import beast.evolution.alignment.TaxonSet;
import beast.math.Binomial;
import org.apache.commons.math.special.Beta;

public class LambdaCoalescentModel extends CalculationNode {

    public Input<RealParameter> alphaInput = new Input<>(
            "alpha",
            "Alpha parameter for Beta-coalescent process",
            Input.Validate.REQUIRED);

    public Input<TaxonSet> taxonSetInput = new Input<>(
            "taxonSet",
            "Taxon set used to define maximum number of extant lineages.");

    public Input<Integer> maxExtantCountInput = new Input<>(
            "maxExtantLineages",
            "Maximum number of extant lineages allowed.",
            Input.Validate.XOR, taxonSetInput);

    int nLeaves;
    RealParameter alpha;
    double[][] cumulativeCoalRates;

    boolean isDirty;

    @Override
    public void initAndValidate() {
        if (taxonSetInput.get() != null)
            nLeaves = taxonSetInput.get().getTaxonCount();
        else
            nLeaves = maxExtantCountInput.get();

        alpha = alphaInput.get();

       cumulativeCoalRates = new double[nLeaves-1][];
        for (int n=2; n<=nLeaves; n++)
            cumulativeCoalRates[n-2] = new double[n-1];

        isDirty = true;
    }

    public double getLogLambda(int n, int k) {
        return Beta.logBeta(k - alpha.getValue(), n - k + alpha.getValue())
                - Beta.logBeta(2-alpha.getValue(), alpha.getValue());
    }

    private void computeCoalRateDistribs() {
        for (int n=2; n<=nLeaves; n++) {
            cumulativeCoalRates[n-2][0] = Math.exp(getLogLambda(n, 2) + Binomial.logChoose(n, 2));

            for (int k=3; k<=n; k++) {
                cumulativeCoalRates[n-2][k-2] = cumulativeCoalRates[n-2][k-3]
                        + Math.exp(getLogLambda(n, k) + Binomial.logChoose(n, k));
            }
        }
    }

    private void update() {
        if (!isDirty)
            return;

        computeCoalRateDistribs();
        isDirty = false;
    }

    public double getTotalCoalRate(int n) {
        update();

        if (n<2)
            return 0;
        else
            return cumulativeCoalRates[n-2][n-2];
    }

    public double[] getCumulativeCoalRateArray(int n) {
        update();

        return cumulativeCoalRates[n-2];
    }

    @Override
    protected void restore() {
        isDirty = true;
        super.restore();
    }

    @Override
    protected boolean requiresRecalculation() {
        isDirty = true;
        return true;
    }
}
