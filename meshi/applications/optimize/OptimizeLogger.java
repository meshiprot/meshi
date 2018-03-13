package meshi.applications.optimize;

import meshi.energy.EvaluationException;
import meshi.energy.TotalEnergy;
import meshi.molecularElements.Protein;
import meshi.molecularElements.extendedAtoms.ResidueExtendedAtoms;
import meshi.molecularElements.extendedAtoms.ResidueExtendedAtomsCreator;
import meshi.optimizers.MCM;
import meshi.scoringFunctions.CombinedEnergyScore;
import meshi.scoringFunctions.Score;
import meshi.sequences.AlignmentException;
import meshi.sequences.ResidueAlignmentMethod;
import meshi.util.*;
import meshi.util.file.MeshiWriter;
import meshi.util.info.ChainsInfo;
import meshi.util.info.MeshiInfo;
import meshi.util.info.MeshiInfoXMLwriter;
import meshi.util.info.ProteinInfoListOld;
import programs.Optimize;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Date;

public class OptimizeLogger extends ProteinInfoListOld implements Logger , OptimizeConstants {
    private ModelAnalyzer analyzer;
    private String parentString;
    private TotalEnergy energy;
    private Protein model, originalModel, nativeStructure;
    OptimizeFiles files;


    public OptimizeLogger(Protein model,
                          OptimizeFiles files, TotalEnergy energy, String parentString, CommandList commands) {
        super("Optimization history of " + model);
        this.model = model;
        this.originalModel = Utils.getProtein(commands, files.modelIn.getAbsolutePath(), ResidueExtendedAtomsCreator.creator, Utils.defaultExceptionHandler);;
        this.files = files;
        if (!(files.nativeStructure.getAbsolutePath().equals("NONE") || files.nativeStructure.getAbsolutePath().equals("none"))) {
            nativeStructure = Utils.getProtein(commands, files.nativeStructure.getAbsolutePath(), ResidueExtendedAtoms.creator, Utils.defaultExceptionHandler);
            if (nativeStructure.chains().size() > 1)
                throw new RuntimeException("Current version does not allow complexes in the native structure");
        } else {
            nativeStructure = null;
        }
        analyzer = new ModelAnalyzer(model, nativeStructure, originalModel, energy, ResidueAlignmentMethod.IDENTITY);
        this.parentString = parentString;
        this.energy = energy;
    }

    public void setEnergy(TotalEnergy energy) {
        this.energy = energy;
        analyzer.setEnergy(energy);
    }

    public double rms() {
        if (nativeStructure == null) return -1;
        else {
            try {
                return analyzer.rms();
            } catch (Exception ex) {
                OptimizeUtils.writeFailureXml(model, ex, files);
                throw new RuntimeException(ex.getMessage());
            }
        }
    }

    public void log(String comment, boolean printFlag) {
        log(comment, printFlag, null);
    }

    public void log(String comment, boolean printFlag, ChainsInfo chainsInfo)  {
        MeshiWriter output;
        MeshiWriter infoTable;
        MeshiInfoXMLwriter infoXML;
        add(analyzer.analyze(comment, chainsInfo));
        if (printFlag) {
            try {
                output = new MeshiWriter(files.modelOut.getAbsolutePath());
                infoXML = new MeshiInfoXMLwriter(files.infoXML.getAbsolutePath());
                infoTable = new MeshiWriter(files.infoTable.getAbsolutePath());
                output.println(parentString);
            } catch (Exception ex) {
                OptimizeUtils.writeFailureXml(model, ex, files);
                throw new RuntimeException(ex);
            }
            try {
                print(output, true);
                print(infoTable, false);
            } catch (IOException ex) {
                OptimizeUtils.writeFailureXml(model, ex, files);
                Utils.throwException(this, ex, "failed to print");
            }
            try {
                print(infoXML);
            } catch (IOException ex) {
                OptimizeUtils.writeFailureXml(model, ex, files);
                Utils.throwException(this, ex, "failed to print");
            }
            energy.on();
            try {
                energy.evaluateAtoms();
            } catch (UpdateableException ex) {
                OptimizeUtils.writeFailureXml(model, ex, files);
                System.out.println("log failed due to " + ex);
                ex.printStackTrace();
                throw new RuntimeException("quiting");
            } catch (EvaluationException ee) {
                OptimizeUtils.writeFailureXml(model, ee, files);
                System.out.println("log failed due to " + ee);
                ee.printStackTrace();
                throw new RuntimeException("quiting");
            }
            Utils.colorByEnergy(model.atoms());
            analyzer.model.atoms().print(output);
            output.close();
            infoTable.close();
            infoXML.close();
        }
    }
    public void mcm(ArrayList<Score> scoreFunctions, TotalEnergy energy, String label, ChainsInfo chainsInfo)  {
        analyzer.setEnergy(energy);
        long time = (new Date()).getTime() - energy.startTime;
        MeshiInfo infoList = energy.energyInfo();
        if (scoreFunctions != null) {
            for (Score scoreFunction : scoreFunctions) {
                Utils.println(" Calculating score " + scoreFunction);


                double rmsFromOriginal;
                try {
                    rmsFromOriginal = Rms.rms(originalModel, model, ResidueAlignmentMethod.IDENTITY);
                } catch (Exception ex) {
                    OptimizeUtils.writeFailureXml(model, ex, files);
                    ex.printStackTrace();
                    throw new RuntimeException(ex);
                }
                ((CombinedEnergyScore) scoreFunction).setChangeElement(rmsFromOriginal);
                MeshiInfo scores = scoreFunction.score(infoList);
                for (MeshiInfo score : scores.flatten())
                    label = label + "\" " + scoreFunction.toString() + "_" + score.type.tag + "=\"" + score.getValue() + " ";
            }
        }
        if (Utils.verbose() || label.startsWith("MCM_END") || label.startsWith("BEGINNING"))
            log(label+"\" time=\""+time, true, chainsInfo);
    }

    public void mcm(ArrayList<Score> scoreFunctions,
                    TotalEnergy energy,
                    int i,
                    MCM.mcmStepResult mcmStepResult)  throws AlignmentException{
        mcm(scoreFunctions, energy, i, mcmStepResult, null);
    }
    public void mcm(ArrayList<Score> scoreFunctions,
                    TotalEnergy energy,
                    int i,
                    MCM.mcmStepResult mcmStepResult,
                    ChainsInfo chainsInfo)  throws AlignmentException{
        //           mcm(optimizationScore, selectionScore, energy, "MCM\"  step=\""+i+"\" score=\""+score.score()+"\" result=\""+mcmStepResult+"\"  lastSuccess=\""+mcmStepResult.lastSuccess());
        mcm(scoreFunctions, energy, "MCM\"  step=\""+i+"\" result=\""+
                mcmStepResult+"\"  lastSuccess=\""+mcmStepResult.lastSuccess(), chainsInfo);
    }
}
