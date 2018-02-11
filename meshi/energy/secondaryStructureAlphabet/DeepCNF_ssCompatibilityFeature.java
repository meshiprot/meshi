package meshi.energy.secondaryStructureAlphabet;

import meshi.energy.AbstractEnergy;
import meshi.energy.EnergyInfoElement;
import meshi.energy.EnergyType;
import meshi.energy.TotalEnergy;
import meshi.molecularElements.Chain;
import meshi.molecularElements.Protein;
import meshi.molecularElements.Residue;
import meshi.molecularElements.atoms.Atom;
import meshi.parameters.DsspLocalStructureLetter;
import meshi.parameters.LocalStructureAlphabetType;
import meshi.sequences.*;
import meshi.util.MeshiAttribute;
import meshi.util.Utils;
import meshi.util.filters.Filter;
import meshi.util.info.InfoType;

import java.util.Random;

/**
 * Created by chen on 23/03/2016.
 */
public class DeepCNF_ssCompatibilityFeature extends AbstractEnergy {

    private static final Filter NON_GAP = new NonGapFilter();
    private static final Filter SS = new SsFilter();
    private SequenceAlignment ss3Alignment;
    private SequenceAlignment ss8Alignment;
    Protein model;
    MeshiSequence sequenceWithPrediction3;
    MeshiSequence sequenceWithPrediction8;
    ResidueSequence residueSequence;


    public DeepCNF_ssCompatibilityFeature(Protein model, MeshiSequence sequenceWithPrediction3,MeshiSequence sequenceWithPrediction8) {
        super(toArray(), new StaticFeaturesInfo(), EnergyType.NON_DIFFERENTIAL);

        if (model.chains().size() > 1)
            throw new RuntimeException("The current version can handle only monomers");

        this.model = model;
        this.sequenceWithPrediction3 = sequenceWithPrediction3;
        this.sequenceWithPrediction8 = sequenceWithPrediction8;
        ss3Alignment = new SequenceAlignment();
        ss8Alignment = new SequenceAlignment();
        for (int chainI = 0; chainI < model.chains().size(); chainI++){
            Chain chain = model.chains().get(chainI);
            residueSequence = chain.sequence();
            ss3Alignment.addAll(SequenceAlignment.identityAlignment(sequenceWithPrediction3, residueSequence));
            ss8Alignment.addAll(SequenceAlignment.identityAlignment(sequenceWithPrediction8, residueSequence));
        }
        evaluate();
    }

    public EnergyInfoElement evaluate() {
        info.setValue(ssCometability(ss3Alignment,LocalStructureAlphabetType.DSSP3));
        //((StaticFeaturesInfo) info).ss3.setValue(ssCometability(ss3Alignment,LocalStructureAlphabetType.DSSP3));
        ((StaticFeaturesInfo) info).ss8.setValue(ssCometability(ss8Alignment,LocalStructureAlphabetType.DSSP7));

        return info;
    }

    private static double ssCometability(SequenceAlignment alignment, LocalStructureAlphabetType localStructureAlphabetType)   {
        double sumPredicted = 0;
        double sumObserved  = 0;
	Random rnd = new Random();
	double num = rnd.nextDouble();

        for (SequenceAlignmentColumn column : alignment) {
            ResidueSsPrediction prediction = (ResidueSsPrediction) column.cell0().getAttribute(MeshiAttribute.SECONDARY_STRUCTURE_ATTRIBUTE);
            Residue residue = (Residue) column.cell1().getAttribute(MeshiAttribute.RESIDUE_ATTRIBUTE);
            num = rnd.nextDouble();

            if (residue != null && prediction != null) {
                double observed = prediction.getProbabilityOf(residue.getLocalStructure().getDSSPReduction(localStructureAlphabetType));
                double predicted = prediction.highestProbability;
                sumObserved += observed;
                sumPredicted += predicted;
            } else {
                String errorMessage = "Something is wrong about the secondary structure prediction. Apparently it has holes in it.\n" +
                        "Residue is " + residue;
                if (Utils.isStrict()) throw new RuntimeException(errorMessage);
                else Utils.println(errorMessage);
            }
        }
            return sumObserved/sumPredicted;


    }


    private static class NonGapFilter implements Filter {
        public boolean accept(Object obj) {
            AlignmentColumn column = (AlignmentColumn) obj;
            return  (!column.cell0().gap()) && (!column.cell1().gap());
        }
    }

    private static class SsFilter implements Filter {
        public boolean accept(Object obj) {
            AlignmentColumn column = (AlignmentColumn) obj;
            ResidueSsPrediction ssp = (ResidueSsPrediction) column.cell0().getAttribute(MeshiAttribute.SECONDARY_STRUCTURE_ATTRIBUTE);
            DsspLocalStructureLetter ss = ssp.secondaryStructure;
            if (ss == DsspLocalStructureLetter.HELIX) return true;
            if (ss == DsspLocalStructureLetter.SHEET) return true;
            return false;
        }
    }

    public void evaluateAtoms(){}
    public void test(TotalEnergy energy, Atom atom) {}





    private static class StaticFeaturesInfo extends EnergyInfoElement{
        public EnergyInfoElement ss3,ss8;
        public StaticFeaturesInfo() {
            super(InfoType.SS_COMPATIBILITY_DEEPCNF3, "deepCNF3Compatibility");
            //ss3 = new EnergyInfoElement(InfoType.SECONDARY_STRUCTURE_COMPATIBILITY_DEEPCNF3, "...");
            ss8 = new EnergyInfoElement(InfoType.SS_COMPATIBILITY_DEEPCNF8, "...");

            getChildren().add(ss8);
            //getChildren().add(ss3);
        }
    }

}
