package pitchfork.models;

import beast.base.core.BEASTObject;
import beast.base.core.Input;
import beast.base.core.Loggable;

import java.io.PrintStream;

public class BetaSkylineLogger extends BEASTObject implements Loggable {

    public Input<AbstractBetaSkylineDistribution> skylineDistributionInput = new Input<>("skylineDistribution",
            "Skyline distribution from which the population function should be " +
                    "logged",
            Input.Validate.REQUIRED);

    public Input<Integer> gridSizeInput = new Input<>("gridSize",
            "Number of points at which population size should be logged.",
            101);

    int gridSize;
    AbstractBetaSkylineDistribution distr;

    @Override
    public void initAndValidate() {
        gridSize = gridSizeInput.get();
        distr = skylineDistributionInput.get();
    }

    @Override
    public void init(PrintStream out) {
        for (int i=0; i<gridSize; i++)
            out.print(getID() + ".popSize" + i + "\t");
    }

    @Override
    public void log(long sample, PrintStream out) {

        double [] popSizes = distr.getPopSizes(gridSize);

        for (int i=0; i<gridSize; i++)
            out.print(popSizes[i] + "\t");
    }

    @Override
    public void close(PrintStream out) {

    }
}
