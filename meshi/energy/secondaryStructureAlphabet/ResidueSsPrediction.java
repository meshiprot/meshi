package meshi.energy.secondaryStructureAlphabet;

import meshi.parameters.DsspLocalStructureLetter;
import meshi.parameters.ResidueType;
import meshi.util.MeshiAttribute;

/**
 * Created by chen on 23/03/2016.
 */
public class ResidueSsPrediction implements MeshiAttribute{
    ResidueType residueType;
    public DsspLocalStructureLetter secondaryStructure;
    double[] probs;
    //probs = [H E C  ,  G I B T S]
    double helixProbability;
    double extendedProbability;
    double coilProbability;


    double gProbability;
    double iProbability;
    double tProbability;
    double sProbability;
    double bProbability;

    double highestProbability;

    public ResidueSsPrediction(String line) {
        class3ProbPsiPred(line);
    }
    public ResidueSsPrediction(String line,String type) {
        String[] words;
        this.probs = new double[8];
        words = line.split("\\s+");
        int t = type.compareTo("DEEPCNF");
        if (words.length == 7 && type.compareTo("DEEPCNF")==0) class3ProbDeepCNF(line);
        else if (words.length == 7 && type.compareTo("PSIPRED")==0) class3ProbPsiPred(line);
        else if (words.length == 12) class8ProbDeepCNF(line);
        else throw new RuntimeException("Error - ResidueSsPrediction - Unknown ssp line.");
    }

    public void class3ProbDeepCNF(String line) {
        String[] words;
        words = line.split("\\s+");

        residueType         = ResidueType.type(words[2]);
        if (words[3].charAt(0) == 'C') secondaryStructure = DsspLocalStructureLetter.COIL;
        else secondaryStructure  = DsspLocalStructureLetter.dsspLocalStructureLetter(1,words[3].charAt(0));

        probs[0] = Double.valueOf(words[4]);//Helix
        probs[1] = Double.valueOf(words[5]);//Extended
        probs[2] = Double.valueOf(words[6]);//Coil


        highestProbability = -1;
        for (int i=0; i<probs.length; i++) {
            if (probs[i]>highestProbability)
                highestProbability = probs[i];
        }
    }
    public void class3ProbPsiPred(String line) {
        String[] words;
        words = line.split("\\s+");

        residueType         = ResidueType.type(words[2]);
        if (words[3].charAt(0) == 'C') secondaryStructure = DsspLocalStructureLetter.COIL;
        else secondaryStructure  = DsspLocalStructureLetter.dsspLocalStructureLetter(1,words[3].charAt(0));

        probs[2] = Double.valueOf(words[4]);//Coil
        probs[0] = Double.valueOf(words[5]);//Helix
        probs[1] = Double.valueOf(words[6]);//Extended


        highestProbability = -1;
        for (int i=0; i<probs.length; i++) {
            if (probs[i]>highestProbability)
                highestProbability = probs[i];
        }
    }

    public void class8ProbDeepCNF(String line) {
        //probabilities are in the order of H G I E B T S  L(loops), the 8 secondary structure types used in DSSP
        //                                  4 5 6 7 8 9 10 11
        //probs = [H E C  ,  G I B T S]
        String[] words;
        words = line.split("\\s+");

        residueType         = ResidueType.type(words[2]);
        if (words[3].charAt(0) == 'L') secondaryStructure = DsspLocalStructureLetter.COIL;
        else secondaryStructure  = DsspLocalStructureLetter.dsspLocalStructureLetter(1,words[3].charAt(0));


        probs[0] = Double.valueOf(words[4]);//Helix
        probs[1] = Double.valueOf(words[7]);//Extended
        probs[2] = Double.valueOf(words[11]);//Coil

        probs[3] = Double.valueOf(words[5]);//G
        probs[4] = Double.valueOf(words[6]);//I
        probs[6] = Double.valueOf(words[9]);//T
        probs[7] = Double.valueOf(words[10]);//S
        probs[5] = Double.valueOf(words[8]);//B

        highestProbability = -1;
        for (int i=0; i<probs.length; i++) {
            if (probs[i]>highestProbability)
                highestProbability = probs[i];
        }
    }


    public ResidueType getResidueType() {return residueType;}

    public double getProbabilityOf(DsspLocalStructureLetter secondaryStructure) {
        //probs = [H E C  ,  G I B T S]
        if (secondaryStructure == DsspLocalStructureLetter.HELIX) return probs[0];
        else if (secondaryStructure == DsspLocalStructureLetter.SHEET) return probs[1];
        else if (secondaryStructure == DsspLocalStructureLetter.COIL) return probs[2];

        else if (secondaryStructure == DsspLocalStructureLetter.THREE10HELIX) return probs[3];
        else if (secondaryStructure == DsspLocalStructureLetter.PIHELIX) return probs[4];
        else if (secondaryStructure == DsspLocalStructureLetter.BETABRIDGE) return probs[5];
        else if (secondaryStructure == DsspLocalStructureLetter.TURN) return probs[6];
        else if (secondaryStructure == DsspLocalStructureLetter.BEND) return probs[7];
        else if (secondaryStructure == DsspLocalStructureLetter.UNK) return 0.0;
        else throw new RuntimeException("Error - ResidueSsPrediction - getProbabilityOf - Unknown DsspLocalStructureLetter.");
    }



    public int key() {
        return MeshiAttribute.SECONDARY_STRUCTURE_ATTRIBUTE;
    }




    public String toString() {
        String out = "ResidueSsPrediction "+residueType+" "+secondaryStructure+" ";
        for (int i=0; i<probs.length; i++) {
                out += probs[i]+" ";
        }
        return out;

    }

}
