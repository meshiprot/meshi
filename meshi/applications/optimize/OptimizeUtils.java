package meshi.applications.optimize;

import meshi.energy.EvaluationException;
import meshi.energy.TotalEnergy;
import meshi.energy.pairwiseNonBondedTerms.atomicPairwisePMFSumma.AtomicPairwisePMFSumma;
import meshi.energy.simpleEnergyTerms.angle.AngleParametersList;
import meshi.energy.simpleEnergyTerms.bond.BondParametersList;
import meshi.energy.simpleEnergyTerms.compositeTorsions.ramachandran.RamachandranEnergy;
import meshi.energy.simpleEnergyTerms.plane.PlaneEnergy;
import meshi.energy.simpleEnergyTerms.plane.PlaneParametersList;
import meshi.energy.simpleEnergyTerms.tether.TetherEnergy;
import meshi.geometry.putH.PutHpos;
import meshi.geometry.putH.PutHposLog;
import meshi.molecularElements.Protein;
import meshi.molecularElements.Residue;
import meshi.molecularElements.ResidueList;
import meshi.molecularElements.atoms.Atom;
import meshi.molecularElements.extendedAtoms.Pro;
import meshi.molecularElements.extendedAtoms.ResidueExtendedAtomsCreator;
import meshi.molecularElements.loops.AtomFinding;
import meshi.optimizers.*;
import meshi.parameters.MeshiPotential;
import meshi.scoringFunctions.Score;
import meshi.sequences.AlignmentException;
import meshi.util.*;
import meshi.util.file.MeshiLineReader;
import meshi.util.file.MeshiWriter;
import meshi.util.info.MeshiInfo;
import meshi.util.info.ProteinInfo;
import programs.Optimize;

import java.io.IOException;
import java.util.ArrayList;

import static meshi.util.KeyWords.MINIMIZE;

public class OptimizeUtils implements OptimizeConstants, KeyWords{
    public static String[] getArguments(String[] keys,String[] args){
        String[] out = new String[keys.length];
        try {
            for (int i = 0; i < keys.length; i++) {
                String key = keys[i];
                out[i] = getArgument(key,args,0);
            }
        }
        catch (MeshiException ex) {
            System.err.print("Required arguments: ");
            for (String key : keys) System.err.print(" "+key);
            System.err.println();
            System.err.print("found\n");
            for (String arg : args) System.err.print(" "+arg);
            System.err.println();
        }

        return out;
    }
    private static String getArgument(String key, String[] args, int start) {
        String temp = null;
        for (int i = start; i < args.length; i++) {
            String fullKey = "-"+key+"=";
            if (args[i].startsWith(fullKey)) {
                temp = args[i].substring(fullKey.length());
                if (i < args.length-1) {
                    if (getArgument(key,args,i+1)!=null)
                        throw new MeshiException("Flag " + key +
                                " appear more then once in command line:\n");
                }
                return temp;
            }
        }
        return null;
    }

    public static void writeFailureXml(Protein model, Exception exception, OptimizeFiles files)  {
        MeshiWriter writer;
        try {
            writer = new MeshiWriter(files.infoXML.getAbsolutePath());
        }
        catch (IOException ex) {
            throw new RuntimeException(ex);
        }
        writer.println("<?xml version=\"1.0\" encoding=\"UTF-8\" ?> ");
        writer.println("<ProteinInfoList name=\"Failure report for Protein: " + model.sourceFile() + "\">");
        writer.print("<ProteinInfo  value=\"MCM_END\" step=\"0\" ");
        writer.println("time=\"0\" fileName=\""+files.modelIn+"\" >");
        writer.println("<Exception>\n"+ exception + "\n</Exception>" );
        writer.println("</ProteinInfo>\n" + "</ProteinInfoList>\n");
        writer.close();
    }

    public static String getParentString(String fileName) throws IOException {
        MeshiLineReader reader = new MeshiLineReader(fileName);
        String out;
        while ((out = reader.readLine()) != null) {
            if (out.startsWith("PARENT")) {
                reader.close();
                return out;
            }
        }
        reader.close();
        reader = new MeshiLineReader(fileName);
        out = "PARENT ";
        String line;
        while ((line = reader.readLine()) != null) {
            if (line.startsWith("REMARK   6 TEMPLATE:"))
                out = out + line.substring(21,28)+" ";
        }
        reader.close();
        if (!out.equals("PARENT "))
            return out;
        return "PARENT N/A";
    }

    public static Protein getModel(CommandList commands, OptimizeFiles optimizeFiles) throws IOException{
        Protein model = Utils.getProtein(commands, optimizeFiles.modelIn.getAbsolutePath(), ResidueExtendedAtomsCreator.creator, Utils.defaultExceptionHandler);;
        Utils.alignToX(model);
        addAtoms(model, commands);
        model.resetBonds();

        //      Secondary structure;
        if (!optimizeFiles.dssp.getAbsolutePath().equals("NONE"))
            Utils.AssignFullDSSP(model,optimizeFiles.dssp.getAbsolutePath());
        Utils.setSS(model, commands);
        Utils.println("secondary structure 1: ");
        int i = 0;
        for (Residue residue : model.residues()) {
            Utils.print(residue.getSecondaryStructure() + "  ");
            if (i++ % 10 == 0) Utils.println();
        }


        return model;
    }

    private static ResidueList getConservedResidues(Protein model, CommandList commands) {
        ResidueList out = new ResidueList();
        String line = commands.firstWord("conserved").secondWord();
        String[] words = line.split(",");
        for (String word : words) {
            int number = Integer.valueOf(word);
            for (Residue residue : model.residues()) {
                if (residue.number() == number) {
                    out.add(residue);
                    Utils.println(residue + " marked as conserved");
                }
            }
        }
        return out;
    }

    public static void addAtoms(Protein model, CommandList commands) throws IOException{
        Command command = commands.firstWord(KeyWords.PARAMETERS_DIRECTORY);
        String parametersDirectory = command.secondWord();
        BondParametersList bondParametersList  = new BondParametersList(parametersDirectory+"/" + MeshiPotential.BOND_PARAMETERS);
        AngleParametersList angleParametersList = new AngleParametersList(parametersDirectory+"/" + MeshiPotential.ANGLE_PARAMETERS);
        PlaneParametersList planeParametersList = new PlaneParametersList(parametersDirectory+"/" + MeshiPotential.PLANE_PARAMETERS);

        boolean hydrogenFailure = false;
        for (Atom atom : model.atoms()) {
            if (atom.isHydrogen() && atom.nowhere()) {
                PutHposLog puthLog = PutHpos.pos(atom, bondParametersList, angleParametersList);
                if (puthLog == null ) hydrogenFailure = true;
            }
        }
        AtomFinding.finalFindAtoms(model,bondParametersList,angleParametersList,planeParametersList);

    }

    private static Optimizer.OptimizerStatus minimize(Protein model, OptimizeFiles files,
                                                      Minimizer lbfgs, TotalEnergy energy,
                                                      TETHER_FLAG tetherFlag,
                                                      RG_FLAG rgFlag, SOLVATE_FLAG solvateFlag) throws OptimizerException {

        return minimize(model, files, lbfgs, energy, RAMACH_COOP_FLAG.ON,
                PROP_COOP_FLAG.ON,
                SUMMA_COOP_FLAG.ON,
                tetherFlag,
                rgFlag,solvateFlag);
    }
    private static Optimizer.OptimizerStatus minimize(Protein model, OptimizeFiles files,
                                                      Minimizer lbfgs,TotalEnergy energy,
                                                      RAMACH_COOP_FLAG rcFlag,
                                                      PROP_COOP_FLAG prFlag,
                                                      SUMMA_COOP_FLAG scFlag,
                                                      TETHER_FLAG tetherFlag,
                                                      RG_FLAG rgFlag,
                                                      SOLVATE_FLAG solvateFlag)  {


        if ((rcFlag == RAMACH_COOP_FLAG.OFF) &&
                (cooperativeZRamachandranCreator.term() != null)) {
            cooperativeZRamachandranCreator.term().off();
            cooperativeZStdRamachandranCreator.term().off();
        }
        if ((prFlag == PROP_COOP_FLAG.OFF) &&
                (cooperativeZPropensityCreator.term() != null)){
            cooperativeZPropensityCreator.term().off();
            cooperativeZStdPropensityCreator.term().off();
        }
        if ((scFlag == SUMMA_COOP_FLAG.OFF) &&
                (cooperativeZSummaCreator.term()!= null)) {
            cooperativeZSummaCreator.term().off();
        }
        if (tetherFlag == TETHER_FLAG.OFF) {
            tetherAllCreator.term().off();
        } else if (tetherFlag == TETHER_FLAG.RESET_OFF) {
            ((TetherEnergy) tetherAllCreator.term()).reset();
            tetherAllCreator.term().off();
        } else if (tetherFlag == TETHER_FLAG.RESET_ON) {
            ((TetherEnergy) tetherAllCreator.term()).reset();
        }
        if (solvateFlag == SOLVATE_FLAG.OFF) {
            solvationCreator.term().off();
        }
        RamachandranEnergy.resetRamachandran(energy);

        try {
            return lbfgs.run();
        } catch (Exception ex) {
            writeFailureXml(model, ex, files);
            ex.printStackTrace();
            throw new RuntimeException(ex);
        }
    }
    public static void minimizeCassette(Protein model, OptimizeFiles files, OptimizeLogger log,
                                         Minimizer lbfgs, TotalEnergy energy, TetherEnergy tether)   {
        Optimizer.OptimizerStatus os = null;
        Utils.println("\nStarting first phase of minimizationCassette \n");

        AtomicPairwisePMFSumma summaTerm = (AtomicPairwisePMFSumma) energy.getEnergyTerm(new AtomicPairwisePMFSumma());
        double summaWeight = summaTerm.weight();
        if (summaTerm.weight() > 1) summaTerm.setWeight(1);

        tether.scaleWeight(0.1);
        ramachCreator.term().scaleWeight(0.1);
        hydrogenBondsPairsCreator.term().scaleWeight(0.1);
        ((PlaneEnergy)planeCreator.term()).scaleWeight(0.01);
        Utils.print("minimizationCassete   1");
        minimize(model, files,lbfgs,energy,
                RAMACH_COOP_FLAG.OFF,
                PROP_COOP_FLAG.OFF,
                SUMMA_COOP_FLAG.OFF,
                TETHER_FLAG.RESET_ON,
                RG_FLAG.ON, SOLVATE_FLAG.OFF);

        if (Utils.verbose())
            log.log("REFINE 0 ", Utils.verbose());
        tether.scaleWeight(10);
        ramachCreator.term().scaleWeight(10);
        hydrogenBondsPairsCreator.term().scaleWeight(10);
        ((PlaneEnergy)planeCreator.term()).scaleWeight(10);
        Utils.print("minimizationCassete   2");
        minimize(model, files,lbfgs,energy,
                RAMACH_COOP_FLAG.OFF,
                PROP_COOP_FLAG.OFF,
                SUMMA_COOP_FLAG.OFF,
                TETHER_FLAG.RESET_ON,
                RG_FLAG.ON, SOLVATE_FLAG.ON);
        if (Utils.verbose())
            log.log("REFINE 0.1 ", Utils.verbose());
        ((PlaneEnergy)planeCreator.term()).scaleWeight(10);
        Utils.print("minimizationCassete   3.");
        minimize(model, files,lbfgs,energy,
                RAMACH_COOP_FLAG.OFF, PROP_COOP_FLAG.OFF,
                SUMMA_COOP_FLAG.OFF, TETHER_FLAG.RESET_ON,
                RG_FLAG.ON, SOLVATE_FLAG.ON);
        if (Utils.verbose())
            log.log("REFINE 0.2 ", Utils.verbose());

        if (cooperativeZSummaCreator.term() != null)
            cooperativeZSummaCreator.term().on();
        double eTest = energy.evaluate();
        Utils.println("Energy with cooperativeSumma = "+eTest);
        //   while (((AtomicPairwisePMFSumma) summaCreator.term()).numberOfClashes()>0) {
        while (eTest > 100000) {
            int nClashes = ((AtomicPairwisePMFSumma) summaCreator.term()).numberOfClashes();
            Utils.println("\n"+nClashes+"clashes\n");
            if (nClashes > 0) Utils.println(((AtomicPairwisePMFSumma) summaCreator.term()).clashes());

            excludedVolCreator.term().scaleWeight(2);

            minimize(model, files,lbfgs,energy,
                    RAMACH_COOP_FLAG.ON, PROP_COOP_FLAG.ON,
                    SUMMA_COOP_FLAG.OFF, TETHER_FLAG.RESET_ON,
                    RG_FLAG.ON, SOLVATE_FLAG.ON);
            if (cooperativeZSummaCreator.term() != null)
                cooperativeZSummaCreator.term().on();
            eTest = energy.evaluate();
            Utils.println("Energy with cooperativeSummaPotential = "+eTest);
        }
        if (excludedVolCreator.term() != null)
            excludedVolCreator.term().off();
        if (Utils.verbose())
            log.log("REFINE 2 ", Utils.verbose());
        Utils.println("\nStarting second phase of minimizationCassette \n");
        try {
            for (int i = 0; (i < 10) & (os != Optimizer.OptimizerStatus.CONVERGED) ; i++) {
                Utils.println(energy.reportHeader());
                Utils.println(energy.report(22200));
                Utils.print("REFINE 2."+i) ;
                os = minimize(model, files,lbfgs, energy,
                        TETHER_FLAG.RESET_ON,
                        RG_FLAG.ON, SOLVATE_FLAG.ON);
                if (Utils.verbose())
                    log.log("REFINE 2."+i, Utils.verbose());
            }
            os = minimize(model, files,lbfgs, energy,
                    TETHER_FLAG.OFF, RG_FLAG.ON, SOLVATE_FLAG.ON);
            if (Utils.verbose())
                log.log("REFINE 3 ", Utils.verbose());
        }
        catch (OptimizerException ex) {
            writeFailureXml(model, ex, files);
            Utils.print("\nException caught\n" + ex + "\n");
            if (Utils.verbose()) energy.test();
            else Utils.throwException("Static method minimizeCassette",ex," minimizeCassette failed.");
        }
        summaTerm.setWeight(summaWeight);
    }

    public static MCM getMCM(Protein model,
                              TotalEnergy minimizationEnergy,
                              ArrayList<Score> scoreFunctions,
                              TotalEnergy perturbationEnergy1,
                              TotalEnergy perturbationEnergy2,
                              TotalEnergy perturbationEnergy3,
                              TotalEnergy perturbationEnergy4,
                              CommandList commands)  {
        ResidueList conservedResidues = getConservedResidues(model, commands);
        Score optimizationScore = null;
        if (scoreFunctions != null) {
            for (Score score:scoreFunctions)
                if (score.toString().equals("optimizationScore"))  {
                    optimizationScore = score;
                    break;
                }
        }
        if (optimizationScore == null)
            throw new RuntimeException("No optimization score.");

        return getMCM(model,minimizationEnergy,scoreFunctions,optimizationScore,
                perturbationEnergy1,perturbationEnergy2,perturbationEnergy3,perturbationEnergy4,
                commands, conservedResidues, MCM.mcmMode.OPTIMIZATION);
    }
    /*   private static Relaxation getRelaxation(Protein model,
                                     TotalEnergy minimizationEnergy,
                                     ArrayList<Score> scoreFunctions,
                                     TotalEnergy perturbationEnergy1,
                                     TotalEnergy perturbationEnergy2,
                                     TotalEnergy perturbationEnergy3,
                                     TotalEnergy perturbationEnergy4,
                                     CommandList commands,
                                     ResidueList conservedResidues) throws UpdateableException,EvaluationException {
               return (Relaxation) getMCM(model, minimizationEnergy, scoreFunctions,
                                                                     perturbationEnergy1, perturbationEnergy2,perturbationEnergy3,perturbationEnergy4,
                                                                     commands, conservedResidues, MCM.mcmMode.RELAXATION);
           }
      */
    private static MCM getMCM(Protein model,
                              TotalEnergy minimizationEnergy,
                              ArrayList<Score> scoreFunctions,
                              Score optimizationScore,
                              TotalEnergy perturbationEnergy1,
                              TotalEnergy perturbationEnergy2,
                              TotalEnergy perturbationEnergy3,
                              TotalEnergy perturbationEnergy4,
                              CommandList commands,
                              ResidueList conservedResidues, MCM.mcmMode mode) {
        Perturbation perturbation1 = new PerturbationByMinimization(perturbationEnergy1, commands, model, conservedResidues, " perturbation with a complete energy function");
        Perturbation perturbation2 = new PerturbationByMinimization(perturbationEnergy2, commands, model, conservedResidues, " perturbation with a minimal energy function");
        Perturbation perturbation3 = new PerturbationByMinimization(perturbationEnergy3, commands, model, conservedResidues, " perturbation with a minimal energy function");
        Perturbation perturbation4 = new PerturbationByMinimization(perturbationEnergy4, commands, model, conservedResidues, " perturbation with a minimal energy function");
        Perturbation[] perturbations = {perturbation1, perturbation2, perturbation3};
        Perturbation perturbation = new CombinedPerturbation(perturbations);
        //new ScmodPerturbations(model, commands,conservedResidues)};
        minimizationEnergy.setCurrent();
        Minimizer minimizer;
        minimizer = Utils.getLBFGS(minimizationEnergy, commands, MINIMIZE);
        double initialTemperature = commands.firstWordFilter(MC_MINIMIZATION).secondWord(INITIAL_TEMPERATURE).thirdWordDouble();
        double finalTemperature = commands.firstWordFilter(MC_MINIMIZATION).secondWord(FINAL_TEMPERATURE).thirdWordDouble();
        int nSteps = commands.firstWordFilter(MC_MINIMIZATION).secondWord(MAX_STEPS).thirdWordInt();
        TemperatureGenerator temperatureGenerator = new TemperatureGenerator(initialTemperature, finalTemperature, nSteps);
        //AbstractEnergy[] excludedTerms = {inflateCreator.term(), tetherAllCreator.term()};
        if (mode == MCM.mcmMode.RELAXATION )
            return new Relaxation(minimizationEnergy, scoreFunctions, optimizationScore, minimizer, perturbation, temperatureGenerator, nSteps);
        return new MCM(minimizationEnergy, scoreFunctions, optimizationScore, minimizer, perturbation, temperatureGenerator, nSteps);
    }

}
