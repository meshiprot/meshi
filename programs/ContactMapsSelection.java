package programs;

import meshi.applications.contactMap.ContactMap;
import meshi.applications.contactMap.ContactMaps;
import meshi.applications.contactMap.ExponentialMeanProbabilityEstimator;
import meshi.applications.optimize.OptimizeEnergies;
import meshi.applications.optimize.OptimizeFiles;
import meshi.applications.optimize.OptimizeUtils;
import meshi.energy.TotalEnergy;
import meshi.molecularElements.Protein;
import meshi.molecularElements.atoms.AtomList;
import meshi.molecularElements.extendedAtoms.ResidueExtendedAtomsCreator;
import meshi.optimizers.LBFGS;
import meshi.scoringFunctions.Score;
import meshi.sequences.ResidueAlignment;
import meshi.sequences.ResidueAlignmentMethod;
import meshi.sequences.ResidueSequence;
import meshi.util.CommandList;
import meshi.util.ModelAnalyzer;
import meshi.util.Utils;
import meshi.util.info.ChainsInfo;
import meshi.util.info.MeshiInfoElement;
import meshi.util.info.ProteinInfoOLd;

import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.ObjectOutputStream;
import java.util.ArrayList;
import java.util.Iterator;


public class ContactMapsSelection {
    public static void main(String[] args) throws Exception {
        double bestMccToNative = -1;
        double topScore = -1;
        double topScoreMcc = -1;
        double topMCCtoPredicted = -1;
        double topMccMcc = -1;
        Protein bestModel = null;

        ContactMaps contactMaps = new ContactMaps(args);
        Iterator<OptimizeFiles> modelsIterator = contactMaps.modelsIterator();
        ContactMap predictedCM = ContactMap.load(contactMaps.predictedCMfileName+".cm");
        while (modelsIterator.hasNext()) {
            OptimizeFiles files = modelsIterator.next();
            Protein model = OptimizeUtils.getModel(contactMaps.commands,files);
            ContactMap modelCM = new ContactMap(model, 8, ContactMap.Mode.CB, ContactMap.Type.PROTEIN);
            ResidueSequence nativeSequence = new ResidueSequence(contactMaps.nativeStructure.residues(), "native structure");
            ResidueSequence modelSequence = new ResidueSequence(model.residues(), "model");
            ResidueAlignment alignment = new ResidueAlignment(nativeSequence, modelSequence);
            double mccToNative    = ContactMap.MCC(contactMaps.nativeCM, modelCM, alignment, 0.5);
            if (mccToNative > bestMccToNative) bestMccToNative = mccToNative;
            double mccToPredicted = ContactMap.MCC(predictedCM, modelCM, alignment, 0.01);
            if (mccToPredicted > topMCCtoPredicted) {
                topMCCtoPredicted = mccToPredicted;
                topMccMcc = mccToNative;
                bestModel = model;
            }
            try {
                ProteinInfoOLd proteinInfo = ContactMaps.evaluate(model, contactMaps.commands);
                double score = -1;
                for (MeshiInfoElement element : proteinInfo)
                    if (element.comment.indexOf(contactMaps.scoreName) != -1)
                        score = element.doubleValue();
                if (score == -1)
                    throw new RuntimeException("Score " + contactMaps.scoreName + " not found.");
                if (score > topScore) {
                    topScore = score;
                    topScoreMcc = mccToNative;
                }
                //predictedCM.estimateProbabilities(new SimpleMeanProbabilityEstimator());
                System.out.println((model.name()+"                    ").substring(0,20) + " " + score + " " + mccToNative + " " + mccToPredicted +" "+topMccMcc+" "+topScoreMcc +" "+bestMccToNative +" "+bestModel.name());
            } catch (Exception ex) {
                //throw ex;
                System.out.println(model.name()+" failed");
            }
        }
    }


}
