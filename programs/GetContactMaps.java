package programs;

import meshi.applications.contactMap.ContactMap;
import meshi.applications.contactMap.ContactMaps;
import meshi.applications.contactMap.ExponentialMeanProbabilityEstimator;
import meshi.applications.optimize.OptimizeEnergies;
import meshi.applications.optimize.OptimizeFiles;
import meshi.applications.optimize.OptimizeUtils;
import meshi.energy.TotalEnergy;
import meshi.molecularElements.Protein;
import meshi.optimizers.LBFGS;
import meshi.scoringFunctions.Score;
import meshi.sequences.ResidueAlignment;
import meshi.sequences.ResidueAlignmentMethod;
import meshi.sequences.ResidueSequence;
import meshi.util.CommandList;
import meshi.util.ModelAnalyzer;
import meshi.util.info.ChainsInfo;
import meshi.util.info.MeshiInfoElement;
import meshi.util.info.ProteinInfoOLd;

import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.ObjectOutputStream;
import java.util.ArrayList;
import java.util.Iterator;


public class GetContactMaps {
    public static void main(String[] args) throws Exception {
        double bestMcc = -1;
        double topScore = -1;
        double topScoreMcc = -1;
        int topScoreNumberOfContacts = -1;
        ContactMap predictedCM;


        ContactMaps contactMaps = new ContactMaps(args);
        predictedCM = new ContactMap(contactMaps.nativeStructure, contactMaps.threshold, contactMaps.mode, ContactMap.Type.PREDICTED);
        predictedCM.reset(0.01);
        //double[] redundancyWeight = new double[contactMaps.modelFiles.length];
        double[] scores = new double[contactMaps.modelFiles.length]; // -1 because of the native structure
        ContactMap[] maps = new ContactMap[contactMaps.modelFiles.length];
        Iterator<OptimizeFiles> modelsIteratorI = contactMaps.modelsIterator();
        int index = 0;
        double [] nContacts = new double[contactMaps.modelFiles.length];
        boolean [] good = new boolean[contactMaps.modelFiles.length];
        while (modelsIteratorI.hasNext()) {
            OptimizeFiles filesI = modelsIteratorI.next();
            Protein modelI = OptimizeUtils.getModel(contactMaps.commands, filesI);
            ContactMap cmI = new ContactMap(modelI, contactMaps.threshold, contactMaps.mode, ContactMap.Type.PROTEIN);
            ResidueSequence seqI = new ResidueSequence(modelI.residues(), "modelI");
            maps[index] = cmI;
            nContacts[index] = cmI.getNumberOfContacts();
            good[index] = true;
            double score = -1;
            try {
                ProteinInfoOLd proteinInfo = evaluate(modelI, contactMaps.commands);
                for (MeshiInfoElement element : proteinInfo)
                    if (element.comment.indexOf(contactMaps.scoreName) != -1)
                        score = element.doubleValue();
                if (score == -1)
                    throw new RuntimeException("Score " + contactMaps.scoreName + " not found.");
            } catch (Exception ex) {
                System.out.println(modelI.name()+" faild in score evaluation.");
                good[index] = false;
            }
            scores[index] = score;
            System.out.println(" number of contacts: "+modelI.name()+" "+cmI.getNumberOfContacts()+" "+score);
            index++;
        }
        for (int k = 0; k < 6; k++) {
            double contactMean = getMean(nContacts, good);
            double contactsStd = getStd(nContacts, good);
            double scoreMean = getMean(scores, good);
            double scoreStd = getStd(scores, good);
            for (int i = 0; i < nContacts.length; i++) {
                if (maps[i].getNumberOfContacts() > contactMean + 4 * contactsStd) {
                    System.out.println(contactMaps.modelFiles[i].getName() + " removed " + maps[i].getNumberOfContacts());
                    good[i] = false;
                } else if (maps[i].getNumberOfContacts() < contactMean - 2 * contactsStd) {
                    System.out.println(contactMaps.modelFiles[i].getName() + " removed " + maps[i].getNumberOfContacts());
                    good[i] = false;
                }
                if (scores[i] < scoreMean - 2 * scoreStd) {
                    System.out.println(contactMaps.modelFiles[i].getName() + " removed " + scores[i]);
                    good[i] = false;
                }
            }
        }


//        for (int i = 0; i < maps.length; i++) {
//            if (!good[i]) continue;
//            ContactMap cmI = maps[i];
//            ResidueSequence seqI = residueSequences[i];
//            for (int j = 0; j < maps.length; j++) {
//                if (!good[j]) continue;
//                ContactMap cmJ = maps[j];
//                ResidueSequence seqJ = residueSequences[j];
//                    double mcc = ContactMap.MCC(cmI, cmJ, new ResidueAlignment(seqI, seqJ), 0);
//                    redundancyWeight[i] += mcc;
//            }
//            System.out.println("sum of MCCs: "+contactMaps.modelFiles[i].getName() + " " + i + " " + redundancyWeight[i]);
//
//        }
//
//        double sum = 0;
//        for (int i = 0; i < redundancyWeight.length; i++) {
//            if (!good[i]) continue;
//            redundancyWeight[i] = 1/redundancyWeight[i];
//            sum += redundancyWeight[i];
//        }
//        for (int i = 0; i < redundancyWeight.length; i++) {
//            if (!good[i]) continue;
//            redundancyWeight[i] = (contactMaps.modelFiles.length)*redundancyWeight[i]/sum;
//        }

        Iterator<OptimizeFiles> modelsIterator = contactMaps.modelsIterator();
        index = 0;
        while (modelsIterator.hasNext()) {
            if (!good[index]) {
                System.out.println(modelsIterator.next().modelIn.getName()+" skiped.");
                index++;
                continue;
            }
            OptimizeFiles files = modelsIterator.next();
            Protein model = OptimizeUtils.getModel(contactMaps.commands, files);
            ContactMap modelCM = maps[index];
            ResidueSequence nativeSequence = new ResidueSequence(contactMaps.nativeStructure.residues(), "native structure");
            ResidueSequence modelSequence = new ResidueSequence(model.residues(), model.name());
            ResidueAlignment alignment = new ResidueAlignment(nativeSequence, modelSequence);
            double mcc = ContactMap.MCC(contactMaps.nativeCM, modelCM, alignment, 0.5);
            if (mcc > bestMcc) bestMcc = mcc;
            double score = scores[index];
                predictedCM.addEvidences(modelCM, alignment, score);
                if (score > topScore) {
                    topScore = score;
                    topScoreMcc = mcc;
                    topScoreNumberOfContacts = modelCM.getNumberOfContacts();
                }
                //predictedCM.estimateProbabilities(new SimpleMeanProbabilityEstimator());
                double threshold = predictedCM.estimateProbabilities(new ExponentialMeanProbabilityEstimator(10));
                if (threshold < 0.5) threshold = 0.5;
                ResidueAlignment alignment1 = new ResidueAlignment(nativeSequence, nativeSequence);
                double predictionMcc = ContactMap.MCC(contactMaps.nativeCM, predictedCM, alignment1, threshold);
                // if (i%10 == 0)
                System.out.println((contactMaps.modelFiles[index]+"                    ").substring(0,20) + " " + score + " "+ mcc + " " + predictionMcc +" "+bestMcc+" "+topScoreMcc +
                        " "+modelCM.getNumberOfContacts()+" "+predictedCM.getNumberOfContacts()+" "+topScoreNumberOfContacts+
                        " "+contactMaps.nativeCM.getNumberOfContacts());
            index++;
        }
        String outFileName = contactMaps.predictedCMfileName;
        predictedCM.print(outFileName+".csv");
        predictedCM.save(outFileName+".CM");
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


    private static double getMean(double[] nContacts, boolean[] good) {
        int n = 0;
        double sum = 0;
        for (int i = 0; i < nContacts.length; i++) {
            if (good[i]){
                sum += nContacts[i];
                n++;
            }
        }
        double mean = sum/n;
        System.out.println("Mean # of contacts: "+mean);
        return mean;
    }

    private static double getStd(double[] nContacts, boolean[] good) {
        int n = 0;
        double sum = 0;
        double sum2 = 0;
        for (int i = 0; i < nContacts.length; i++) {
            if (good[i]){
                sum += nContacts[i];
                sum2 += nContacts[i]*nContacts[i];
                n++;
            }
        }
        double mean = sum/n;
        double sum2mean = sum2/n;
        double variance = sum2mean - mean*mean;
        double std = Math.sqrt(variance);
        System.out.println("STD # of contacts: "+std);
        return std;
    }

}
