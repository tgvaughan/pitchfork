package pitchfork.models;

import beast.core.parameter.RealParameter;
import beast.evolution.alignment.Taxon;
import beast.evolution.alignment.TaxonSet;
import beast.evolution.tree.RandomTree;
import beast.evolution.tree.TraitSet;
import beast.evolution.tree.Tree;
import beast.evolution.tree.coalescent.PopulationFunction;
import beast.math.statistic.DiscreteStatistics;
import beast.util.Randomizer;
import junit.framework.Assert;
import pitchfork.LBSPTest;
import org.junit.Test;

import java.util.ArrayList;
import java.util.List;
import java.util.function.Function;

public class SimulatedLambdaCoalescentTreeTest extends LBSPTest {

    private static SimulatedLambdaCoalescentTree getSimulatedLambdaCoalescentTree(int nLeaves, double alpha,
                                                                                  PopulationFunction populationFunction) {

        List<Taxon> taxonList = new ArrayList<>();
        StringBuilder traitValueBuilder = new StringBuilder();

        for (int i=1; i<=nLeaves; i++) {
            String id = "t" + i;
            taxonList.add(new Taxon(id));

            if (i>1)
                traitValueBuilder.append(",");

            traitValueBuilder.append(id).append("=0.0");
        }
        TaxonSet taxonSet = new TaxonSet(taxonList);

        TraitSet dateTrait = new TraitSet();
        dateTrait.initByName( "traitname", "date-backward",
                "taxa", taxonSet,
                "value", traitValueBuilder.toString());

        LambdaCoalescentModel lcModel = new LambdaCoalescentModel();
        lcModel.initByName("alpha", new RealParameter(String.valueOf(alpha)),
                "taxonSet", taxonSet);

        SimulatedLambdaCoalescentTree tree = new SimulatedLambdaCoalescentTree();
        tree.initByName(
                "model", lcModel,
                "populationFunction", populationFunction,
                "trait", dateTrait);

        return tree;
    }

    private static RandomTree getSimulatedKingmanCoalescentTree(int nLeaves, PopulationFunction populationFunction) {

        List<Taxon> taxonList = new ArrayList<>();
        StringBuilder traitValueBuilder = new StringBuilder();

        for (int i=1; i<=nLeaves; i++) {
            String id = "t" + i;
            taxonList.add(new Taxon(id));

            if (i>1)
                traitValueBuilder.append(",");

            traitValueBuilder.append(id).append("=0.0");
        }

        TaxonSet taxonSet = new TaxonSet(taxonList);
        TraitSet dateTrait = new TraitSet();
        dateTrait.initByName( "traitname", "date-backward",
                "taxa", taxonSet,
                "value", traitValueBuilder.toString());

        RandomTree tree = new RandomTree();
        tree.initByName(
                "populationModel", populationFunction,
                "trait", dateTrait,
                "taxonset", taxonSet);

        return tree;
    }

    private static <T> List<T> performSimulation(Tree simulatedTree,
                                     Function<Tree,T> summaryFunction) {

        int nSims = 100000;

        List<T> results = new ArrayList<>();

        for (int i = 0; i< nSims; i++) {
            if (i>0)
                simulatedTree.initAndValidate();

            results.add(summaryFunction.apply(simulatedTree));
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
    public void testAlpha2CoalescentTimes2Taxon() {
        Randomizer.setSeed(2);

        Tree tree = getSimulatedLambdaCoalescentTree(2, 1.5,
                getConstantPopulation(1.0));
        List<Double> times = performSimulation(tree, t -> t.getRoot().getHeight());

        double mean = getListMean(times);
        double var = getListVariance(times);


        System.out.println("Mean age: " + mean);
        Assert.assertEquals(1.0, mean, 1e-2);

        System.out.println("Variance: " + var);
        Assert.assertEquals(1.0, var, 1e-2);
    }


    @Test
    public void testAlpha2CoalescentTimes10Taxon() {
        Randomizer.setSeed(1);

        Tree tree = getSimulatedLambdaCoalescentTree(10, 1.999999,
                getExponentialPopulation(1.0, 10.0));
        List<Double> times = performSimulation(tree, t -> t.getRoot().getHeight());
        double mean = getListMean(times);
        double var = getListVariance(times);

        Tree kingmanCoalescentTree = getSimulatedKingmanCoalescentTree(10,
                getExponentialPopulation(1, 10.0));
        List<Double> kingmanTimes = performSimulation(kingmanCoalescentTree, t->t.getRoot().getHeight());
        double kingmanMean = getListMean(kingmanTimes);
        double kingmanVar = getListVariance(kingmanTimes);

        System.out.println("Mean age: " + mean);
        System.out.println("Mean age (Kingman): " + kingmanMean);
        Assert.assertEquals(kingmanMean, mean, kingmanMean*1e-2);

        System.out.println("Variance: " + var);
        System.out.println("Variance (Kingman): " + kingmanVar);
        Assert.assertEquals(kingmanVar, var, kingmanVar*1e-2);
    }
}
