package lambdabsp.model;

import beast.core.parameter.RealParameter;
import beast.evolution.alignment.Taxon;
import beast.evolution.alignment.TaxonSet;
import beast.evolution.tree.TraitSet;
import beast.evolution.tree.Tree;
import beast.evolution.tree.coalescent.ConstantPopulation;
import beast.math.statistic.DiscreteStatistics;
import junit.framework.Assert;
import org.junit.Test;

import java.util.ArrayList;
import java.util.List;
import java.util.function.Function;

public class SimulatedLambdaCoalescentTreeTest {

    private static <T> List<T> performSimulation(double alpha, int nLeaves,
                                     Function<Tree,T> summaryFunction) {

        int nSims = 100000;

        List<T> results = new ArrayList<>();

        ConstantPopulation popFun = new ConstantPopulation();
        popFun.initByName("popSize", new RealParameter("1.0"));

        List<Taxon> taxonList = new ArrayList<>();
        StringBuilder traitValueBuilder = new StringBuilder();

        for (int i=1; i<=nLeaves; i++) {
            String id = "t" + i;
            taxonList.add(new Taxon(id));

            if (i>1)
                traitValueBuilder.append(",");

            traitValueBuilder.append(id).append("=0.0");
        }

        TraitSet dateTrait = new TraitSet();
        dateTrait.initByName( "traitname", "date-backward",
                "taxa", new TaxonSet(taxonList),
                "value", traitValueBuilder.toString());

        SimulatedLambdaCoalescentTree tree = new SimulatedLambdaCoalescentTree();
        tree.initByName(
                "alpha", new RealParameter(String.valueOf(alpha)),
                "populationFunction", popFun,
                "trait", dateTrait);

        for (int i = 0; i< nSims; i++) {
            if (i>0)
                tree.initAndValidate();

            results.add(summaryFunction.apply(tree));
        }

        return results;
    }

    private static double getListMean(List<Double> values) {
        double[] array = new double[values.size()];
        for (int i=0; i<values.size(); i++)
            array[i] = values.get(i);

        return DiscreteStatistics.mean(array);
    }

    private static double getListVariance(List<Double> values) {
        double[] array = new double[values.size()];
        for (int i=0; i<values.size(); i++)
            array[i] = values.get(i);

        return DiscreteStatistics.variance(array);
    }

    @Test
    public void testKingmanCoalescentTimes() {

        List<Double> times = performSimulation(1.5, 2,
                t -> t.getRoot().getHeight());

        double mean = getListMean(times);
        double var = getListVariance(times);


        System.out.println("Mean age: " + mean);
        Assert.assertEquals(1.0, mean, 1e-2);

        System.out.println("Variance: " + var);
        Assert.assertEquals(1.0, var, 1e-2);
    }
}
