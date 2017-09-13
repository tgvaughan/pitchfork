package lambdabsp.model;

import beast.core.Input;
import beast.core.Runnable;
import beast.core.parameter.RealParameter;
import beast.evolution.tree.coalescent.PopulationFunction;

public class LambdaCoalescentTreeSimulator extends Runnable {

    public Input<RealParameter> alphaInput = new Input<>(
            "alpha",
            "Parameter determining shape of Lambda distribution",
            Input.Validate.REQUIRED);

    public Input<PopulationFunction> populationFunctionInput = new Input<>(
            "populationFunction",
            "Population function used in coalescent simulation.",
            Input.Validate.REQUIRED);

    public Input<String> fileNameInput =  new Input<>(
            "fileName",
            "Name of file to save Newick representation of tree to.",
            Input.Validate.REQUIRED);

    @Override
    public void initAndValidate() { }

    @Override
    public void run() throws Exception {

        SimulatedLambdaCoalescentTree tree = new SimulatedLambdaCoalescentTree();
        tree.initByName(
                "alpha", alphaInput.get(),
                "populationFunction", populationFunctionInput.get(),
                "fileName", fileNameInput.get());


    }
}
