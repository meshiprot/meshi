package meshi.applications.optimize;

import meshi.energy.EnergyCreator;
import meshi.energy.TotalEnergy;
import meshi.energy.simpleEnergyTerms.tether.TetherEnergy;
import meshi.molecularElements.Protein;
import meshi.molecularElements.extendedAtoms.Pro;
import meshi.scoringFunctions.GdtScore;
import meshi.scoringFunctions.Score;
import meshi.util.CommandList;

import java.util.ArrayList;

import static meshi.applications.optimize.OptimizeConstants.*;

public class OptimizeEnergies {
    public final TotalEnergy perturbationEnergy1, perturbationEnergy2, perturbationEnergy3, perturbationEnergy4;
    public final TotalEnergy minimizationEnergy;
    public final ArrayList<Score> scoreFunctions;

    public OptimizeEnergies(Protein model, CommandList commands) {
        perturbationEnergy1 = new TotalEnergy(model,
                excludeEnergy(energyCreators, excludeFromPerturbation1),
                commands,"perturbationEnergy1");
        perturbationEnergy2 = new TotalEnergy(model,
                excludeEnergy(energyCreators, excludeFromPerturbation2),
                commands,"perturbationEnergy2");
        perturbationEnergy3 = new TotalEnergy(model,
                excludeEnergy(energyCreators, excludeFromPerturbation3),
                commands,"perturbationEnergy3");
        perturbationEnergy4 = new TotalEnergy(model,
                excludeEnergy(energyCreators, excludeFromPerturbation4),
                commands,"perturbationEnergy4");
        minimizationEnergy = new TotalEnergy(model,
                excludeEnergy(energyCreators, excludeFromMinimization),
                commands,"minimization energy");
        //Instantiation TotalEnergy object may kill the other ones. Here we want all of them to be simultaneously alive (though not used).
        perturbationEnergy1.terminator.reset();
        perturbationEnergy2.terminator.reset();
        perturbationEnergy3.terminator.reset();
        perturbationEnergy4.terminator.reset();
        minimizationEnergy.terminator.reset();
        ((TetherEnergy) residueTetherCreator.term()).setResidue(commands);
        scoreFunctions = GdtScore.getScoreFunctions(commands);
        ArrayList<Score> energyScores = new ArrayList<Score>();
        if (scoreFunctions != null) {
            for (Score scoreFunction : scoreFunctions)
                if (scoreFunctions.toString().equals("energy"))
                    energyScores.add(scoreFunction);
        }
    }
    private static EnergyCreator[] excludeEnergy(EnergyCreator[] source, EnergyCreator[] exclude) {
        EnergyCreator[] out = new EnergyCreator[source.length - exclude.length];
        int j = 0;
        for (int i = 0; i < source.length; i++) {
            if (notFound(source[i], exclude)) {
                out[j] = source[i];
                j++;
            }
        }
        return out;
    }

    private static boolean notFound(EnergyCreator energyCreator, EnergyCreator[] energyCreators) {
        for (int i = 0; i < energyCreators.length; i++) {
            if (energyCreators[i].equals(energyCreator)) return false;
        }
        return true;
    }
}
