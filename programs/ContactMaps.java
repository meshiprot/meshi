package programs;

import meshi.applications.contactMap.ContactMap;
import meshi.applications.contactMap.ExponentialMeanProbabilityEstimator;
import meshi.applications.contactMap.SimpleMeanProbabilityEstimator;
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



public class ContactMaps {
    public static void main(String[] args) throws Exception{
        double bestMcc = -1;
        double topScore = -1;
        double topScoreMcc = -1;
        ArrayList<ContactMap> maps = new ArrayList<>();
        int topScoreNumberOfContacts = -1;
        CommandList commands = new CommandList(args[0]);
        if (commands.keyExists("verbose")) Utils.verboseOn();
        else Utils.verboseOff();
        String scoreName = args[4];
        Protein nativeStructure = new Protein(new AtomList(args[1]), ResidueExtendedAtomsCreator.creator);
        ContactMap nativeCM = new ContactMap(nativeStructure, 8, ContactMap.Mode.CB, ContactMap.Type.PROTEIN);
        ContactMap predictedCM = new ContactMap(nativeStructure, 8, ContactMap.Mode.CB, ContactMap.Type.PREDICTED);
        maps.add(predictedCM);
        predictedCM.reset(0.01);
        nativeCM.print("nativeCM.csv");
        File modelsDir = new File(args[2]);
        File[] modelFiles = modelsDir.listFiles();
        File dsspDir = new File(args[3]);
        File[] dsspFiles = dsspDir.listFiles();
        for (int i = 0; i < modelFiles.length; i++) {
            String fileName = modelFiles[i].getAbsolutePath();
            if (fileName.endsWith("pdb")) {
                if (fileName.indexOf(".N.") >= 0) continue;
                String dsspFileName = getDsspFile(fileName, dsspFiles);
                OptimizeFiles files = new OptimizeFiles(fileName, dsspFileName,"", "");
                Protein model = OptimizeUtils.getModel(commands,files);
                ContactMap modelCM = new ContactMap(model, 8, ContactMap.Mode.CB, ContactMap.Type.PROTEIN);
                maps.add(modelCM);
                ResidueSequence nativeSequence = new ResidueSequence(nativeStructure.residues(), "native structure");
                ResidueSequence modelSequence = new ResidueSequence(model.residues(), "model");
                ResidueAlignment alignment = new ResidueAlignment(nativeSequence, modelSequence);
                double mcc = ContactMap.MCC(nativeCM, modelCM, alignment, 0.5);
                if (mcc > bestMcc) bestMcc = mcc;
                try {
                    ProteinInfoOLd proteinInfo = evaluate(model, commands);
                    double score = -1;
                    for (MeshiInfoElement element : proteinInfo)
                        if (element.comment.indexOf(scoreName) != -1)
                            score = element.doubleValue();
                    if (score == -1)
                        throw new RuntimeException("Score " + scoreName + " not found.");
                    predictedCM.addEvidences(modelCM, alignment, score);
                    if (score > topScore) {
                        topScore = score;
                        topScoreMcc = mcc;
                        topScoreNumberOfContacts = modelCM.getNumberOfContacts();
                    }
                    //predictedCM.estimateProbabilities(new SimpleMeanProbabilityEstimator());
                    double threshold = predictedCM.estimateProbabilities(new ExponentialMeanProbabilityEstimator(5));
                    if (threshold < 0.5) threshold = 0.5;
                    ResidueAlignment alignment1 = new ResidueAlignment(nativeSequence, nativeSequence);
                    double predictionMcc = ContactMap.MCC(nativeCM, predictedCM, alignment1, threshold);
                   // if (i%10 == 0)
                        System.out.println(model.name() + " " + score + " " + mcc + " " + predictionMcc +" "+bestMcc+" "+topScoreMcc +
                                       " "+modelCM.getNumberOfContacts()+" "+predictedCM.getNumberOfContacts()+" "+topScoreNumberOfContacts+
                                       " "+nativeCM.getNumberOfContacts());
                } catch (Exception ex) {
                    //throw ex;
                    System.out.println(model.name()+" failed");
                }
            }
        }
        predictedCM.print("predictedCM.csv");
        save("contactMaps.obj", maps);
    }

    public static void save(String fileName, ArrayList<ContactMap> maps) throws IOException{
        FileOutputStream fos = new FileOutputStream(fileName);
        ObjectOutputStream oos = new ObjectOutputStream(fos);
        oos.writeObject(maps);
        oos.close();
    }


    private static ProteinInfoOLd evaluate(Protein model, CommandList commands) {
        OptimizeEnergies energies = new OptimizeEnergies(model, commands);
        TotalEnergy energy = energies.minimizationEnergy;
        LBFGS lbfgs;
        try {
            lbfgs = new LBFGS(energy, 0.01, 3000, 200);
        } catch (Exception ex) {
            throw new RuntimeException(ex);
        }
        //lbfgs.run();
        energy.setCurrent();
        ArrayList<Score> scoreFunctions = energies.scoreFunctions;
        ModelAnalyzer modelAnalyzer = new ModelAnalyzer(model, null, null, energy, scoreFunctions, ResidueAlignmentMethod.IDENTITY);
        ChainsInfo chainsInfo = new ChainsInfo(model);
        return modelAnalyzer.analyze("Evaluate "+model.name(), chainsInfo);
    }

    private static String getDsspFile(String modelFile, File[] files){
        String modelName;
        int tsIndex = modelFile.indexOf("_TS");
        int alIndex = modelFile.indexOf("_AL");
        int nIndex  = modelFile.indexOf("N.scwrl");
        if (alIndex > tsIndex) tsIndex = alIndex;
        if (nIndex  > tsIndex) tsIndex = nIndex;
        int dotIndex = tsIndex;
        while (modelFile.charAt(dotIndex) != '.' & dotIndex > -1) dotIndex--;
        if (dotIndex <= -1)
            throw new RuntimeException("This is weird.");
        modelName = modelFile.substring(dotIndex+1,tsIndex);
        for (File file : files) {
            String fileName = file.getAbsolutePath();
            if (fileName.indexOf(modelName) != -1) return fileName;
        }
        throw new RuntimeException("Dssp file not found "+modelName);
    }

}
