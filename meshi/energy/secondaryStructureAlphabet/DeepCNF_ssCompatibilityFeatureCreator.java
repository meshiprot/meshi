package meshi.energy.secondaryStructureAlphabet;

import meshi.energy.AbstractEnergy;
import meshi.energy.EnergyCreator;
import meshi.geometry.DistanceMatrix;
import meshi.molecularElements.Protein;
import meshi.sequences.MeshiSequence;
import meshi.sequences.SequenceAlignmentCell;
import meshi.util.CommandList;
import meshi.util.dssp.DSSPFull;
import meshi.util.file.MeshiLineReader;
import meshi.util.info.InfoType;

import java.io.IOException;

/**
 * Created by chen on 23/03/2016.
 *
 * ssCompetability
 */
public class DeepCNF_ssCompatibilityFeatureCreator extends EnergyCreator{
    MeshiSequence Ss3Sequence;
    MeshiSequence Ss8Sequence;

    public DeepCNF_ssCompatibilityFeatureCreator() {
        super(InfoType.SS_COMPATIBILITY_DEEPCNF3);


    }
    // The protein model is already assigned a DSSPFull structure for each residue
    public AbstractEnergy createEnergyTerm(Protein model, DistanceMatrix dm, CommandList commands) {
        String ssPredictionFileName = getDeepCNF3FileName(commands);
        Ss3Sequence = getDeepCNFSequence(ssPredictionFileName);

        ssPredictionFileName = getDeepCNF8FileName(commands);
        Ss8Sequence = getDeepCNFSequence(ssPredictionFileName);

        term =  new DeepCNF_ssCompatibilityFeature(model, Ss3Sequence,Ss8Sequence);
        return term;

    }

    private static String getDeepCNF3FileName(CommandList commands) {
        return commands.firstWord("deepCNF3FileName").secondWord();
    }

    private static String getDeepCNF8FileName(CommandList commands) {
        return commands.firstWord("deepCNF8FileName").secondWord();
    }

    private static MeshiSequence[] getDsspMeshiSequence(DSSPFull dssp){
        if (dssp.length() == 0) throw new RuntimeException("problem wth DSSP file");
        return dssp.getResidueSequenceWithFullSs();

    }

    private static MeshiSequence getDeepCNFSequence(String inFile){
        MeshiSequence SsSequence = new MeshiSequence("Sequence of SS3 prediction (DeepCNF)");
        MeshiLineReader reader = new MeshiLineReader(inFile);
        String line;
        int i = 0;
        String line1, line2;
        try {
            while ((line = reader.readLine()) != null) {
                if (line.startsWith("#")) continue;
                if (line.length() < 10) continue;
                ResidueSsPrediction residueSsPrediction = new ResidueSsPrediction(line,"DEEPCNF");
                SequenceAlignmentCell cell = new SequenceAlignmentCell(residueSsPrediction.residueType.nameOneLetter().charAt(0), i);
                i++;
                SsSequence.add(cell);
                cell.addAttribute(residueSsPrediction);
            }

            reader.close();
        } catch (IOException ex){
            ex.printStackTrace();
            throw new RuntimeException(ex.getMessage());
        }
        return SsSequence;
    }


}
