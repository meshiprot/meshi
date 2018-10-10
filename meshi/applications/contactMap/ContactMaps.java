package meshi.applications.contactMap;

import meshi.applications.optimize.OptimizeEnergies;
import meshi.applications.optimize.OptimizeFiles;
import meshi.energy.TotalEnergy;
import meshi.molecularElements.Protein;
import meshi.molecularElements.atoms.AtomList;
import meshi.molecularElements.extendedAtoms.ResidueExtendedAtomsCreator;
import meshi.optimizers.LBFGS;
import meshi.scoringFunctions.Score;
import meshi.sequences.ResidueAlignmentMethod;
import meshi.util.CommandList;
import meshi.util.ModelAnalyzer;
import meshi.util.Utils;
import meshi.util.info.ChainsInfo;
import meshi.util.info.ProteinInfoOLd;
import programs.Optimize;

import java.io.*;
import java.util.ArrayList;
import java.util.Iterator;

public class ContactMaps {
    public final CommandList commands;
    public final String scoreName;
    public final double threshold;
    public final ContactMap.Mode mode;
    public final File modelsDir;
    public final File[] modelFiles;
    public final File dsspDir;
    public final File[] dsspFiles;
    public final Protein nativeStructure;
    public final ContactMap nativeCM;
    public final String predictedCMfileName;

    public ContactMaps(String[] args){

        commands = new CommandList(args[0]);
        if (commands.keyExists("verbose")) Utils.verboseOn();
        else Utils.verboseOff();

        threshold = getDouble(args[1]);
        mode = getMode(args[2]);
        nativeStructure = new Protein(new AtomList(args[3]), ResidueExtendedAtomsCreator.creator);
        nativeCM = new ContactMap(nativeStructure, threshold, mode, ContactMap.Type.PROTEIN);
        nativeCM.print("nativeCM.csv");
        modelsDir = new File(args[4]);
        File[] tmp = modelsDir.listFiles();
        int count = 0;
        for (File file : tmp)
            if (file.getName().endsWith(".pdb") & file.getName().indexOf(".N.") == -1)
                count++;
        modelFiles = new File[count];
        int index = 0;
        for (File file : tmp)
            if (file.getName().endsWith(".pdb") & file.getName().indexOf(".N.") == -1) {
                modelFiles[index] = file;
                index++;
            }
        dsspDir = new File(args[5]);
        dsspFiles = dsspDir.listFiles();
        scoreName = args[6];
        predictedCMfileName = "predictedCM_"+threshold+"_"+mode+"_"+nativeStructure.name();

    }

    public Iterator<OptimizeFiles> modelsIterator() {
        return new ModelsIterator(this);
    }

    private static class ModelsIterator implements Iterator<OptimizeFiles> {
        Iterator<File> filesIterator;
        ContactMaps contactMaps;
        OptimizeFiles next;
        public ModelsIterator(ContactMaps contactMaps) {
            this.contactMaps = contactMaps;
            ArrayList<File> tmp1 = new ArrayList<>();
            for (File file : contactMaps.modelFiles)
                tmp1.add(file);
            filesIterator = tmp1.iterator();

            next = null;
            while ((next == null) & (filesIterator.hasNext()) ){
                File modelFile = filesIterator.next();
                String fileName = modelFile.getAbsolutePath();
                if ((fileName.endsWith("pdb")) & (fileName.indexOf(".N.") == -1)) {
                    String dsspFileName = getDsspFile(fileName, contactMaps.dsspFiles);
                    next = new OptimizeFiles(fileName, dsspFileName,"", "");
                }
            }
        }

        public boolean hasNext() { return next != null;}

        public OptimizeFiles next() {
            if (! hasNext())
                throw new RuntimeException("No more items.");
            OptimizeFiles out = next;
            next = null;
            while ((next == null) & (filesIterator.hasNext()) ){
                File modelFile = filesIterator.next();
                String fileName = modelFile.getAbsolutePath();
                if ((fileName.endsWith("pdb")) & (fileName.indexOf(".N.") == -1)) {
                    String dsspFileName = getDsspFile(fileName, contactMaps.dsspFiles);
                    next = new OptimizeFiles(fileName, dsspFileName,"", "");
                }
            }
            return out;
        }

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



    private double getDouble(String s) {
        try {
            return (new Double(s)).doubleValue();
        } catch (Exception ex) {
            throw new RuntimeException("Filed to create a double from " + s + "\n"+ex);
        }
    }

    private ContactMap.Mode getMode(String s) {
        for (ContactMap.Mode mode : ContactMap.Mode.values())
            if (mode.toString().equals(s))
                return mode;
        throw new RuntimeException("Filed to find mode " + mode);
    }


    public static void save(String fileName, ArrayList<ContactMap> maps) throws IOException {
        FileOutputStream fos = new FileOutputStream(fileName);
        ObjectOutputStream oos = new ObjectOutputStream(fos);
        oos.writeObject(maps);
        oos.close();
    }


    public static ProteinInfoOLd evaluate(Protein model, CommandList commands) {
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



}
