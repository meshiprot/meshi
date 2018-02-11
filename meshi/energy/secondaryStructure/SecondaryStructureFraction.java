package meshi.energy.secondaryStructure;

import meshi.energy.AbstractEnergy;
import meshi.energy.EnergyInfoElement;
import meshi.energy.EnergyType;
import meshi.energy.TotalEnergy;
import meshi.molecularElements.Protein;
import meshi.molecularElements.Residue;
import meshi.molecularElements.atoms.Atom;
import meshi.parameters.SecondaryStructure;
import meshi.scoringFunctions.Sigmoid;
import meshi.util.info.InfoType;

/**
 * Created by IntelliJ IDEA.
 * User: chen
 * Date: 20/05/12
 * Time: 11:26
 * To change this template use File | Settings | File Templates.
 */
public class SecondaryStructureFraction extends AbstractEnergy {
    private final double SS_FRACTION_SATURATION_MIDPOINT = 0.4;
    private final double SS_FRACTION_SATURATION_SLOPE = 100;
    private final double SS_SATURATION_MIDPOINT = 0.5;
    private final double SS_SATURATION_SLOPE = 10;
    Protein protein;
    public SecondaryStructureFraction(EnergyInfoElement info, Protein protein){
        super(toArray(),info, EnergyType.NON_DIFFERENTIAL);
        this.protein = protein;
    }

    public EnergyInfoElement evaluate() {
        // sigmoid parameters - maps [0 1] to [0 1] with 0.4 mapped to 0.5
        int sum = 0;
        double sumSS = 0;
        double sumHelix = 0;
        double sumSheet = 0;
        for (Residue residue : protein.residues()) {
            if (!residue.dummy()) {
                sum++;
                if (residue.getSecondaryStructure() != SecondaryStructure.COIL)
                    sumSS += 1;
                if (residue.getSecondaryStructure() == SecondaryStructure.HELIX)
                    sumHelix += 1;
                if (residue.getSecondaryStructure() == SecondaryStructure.SHEET)
                    sumSheet += 1;
            }
        }
        double ssFraction = sumSS/sum;
        info().setValue(ssFraction);
        EnergyInfoElement helixFraction = new EnergyInfoElement(InfoType.HELIX_FRACTION, "Fraction of alpha helices among SS elements");
        EnergyInfoElement sheetFraction = new EnergyInfoElement(InfoType.SHEET_FRACTION, "Fraction of beta sheet among SS elements");
        EnergyInfoElement saturatedSSfraction = new EnergyInfoElement(InfoType.SATURATED_SECONDARY_STRUCTURE_FRACTION, "A non-linear version of SECONDARY_STRUCTURE_FRACTION");
        info.getChildren().add(saturatedSSfraction);
        info.getChildren().add(helixFraction);
        info.getChildren().add(sheetFraction);
        Sigmoid sigmoid = new Sigmoid(SS_SATURATION_MIDPOINT, SS_SATURATION_SLOPE);
        double helixF;
        if (sumHelix == 0) helixF = 0.05; // for numerical stability, in case sumSheet is also zero.
        else helixF = sigmoid.sigmoid(sumHelix/(sumSheet+sumHelix))*0.9+0.05;
        helixFraction.setValue(helixF);
        double sheetF;
        if (sumSheet == 0) sheetF = 0.05;
        else sheetF = sigmoid.sigmoid(sumSheet/(sumSheet+sumHelix))*0.0+0.05;
        sheetFraction.setValue(sheetF);
        sigmoid = new Sigmoid(SS_FRACTION_SATURATION_MIDPOINT, SS_FRACTION_SATURATION_SLOPE);
        saturatedSSfraction.setValue(sigmoid.sigmoid(ssFraction));
        return info();
    }

    public void evaluateAtoms() {}
    public void test(TotalEnergy te, Atom a) {
             System.out.println("SecondaryStructureFraction tested.");
     }

}

