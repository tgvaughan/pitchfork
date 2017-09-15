package lambdabsp.operators;

import beast.evolution.operators.TreeOperator;

/**
 * Operator to traverse the space of multifurcating phylogenetic trees.
 */
public class multifurcationOperartor extends TreeOperator{

    @Override
    public void initAndValidate() { }

    @Override
    public double proposal() {
        return 0;
    }
}
