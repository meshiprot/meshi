/*
 * Copyright (c) 2010. of Chen Keasar, BGU . For free use under LGPL
 */

package programs;

import meshi.energy.EvaluationException;
import meshi.energy.secondaryStructureAlphabet.ResidueSsPrediction;
import meshi.energy.simpleEnergyTerms.angle.AngleParametersList;
import meshi.energy.simpleEnergyTerms.bond.BondParametersList;
import meshi.energy.simpleEnergyTerms.inflate.inflate.InflateCreator;
import meshi.energy.simpleEnergyTerms.inflate.inflate.InflateType;
import meshi.energy.simpleEnergyTerms.plane.PlaneParametersList;
import meshi.geometry.ArcCos;
import meshi.geometry.putH.PutHpos;
import meshi.geometry.putH.PutHposLog;
import meshi.molecularElements.Protein;
import meshi.molecularElements.atoms.Atom;
import meshi.molecularElements.extendedAtoms.ResidueExtendedAtomsCreator;
import meshi.molecularElements.loops.AtomFinding;
import meshi.optimizers.OptimizerException;
import meshi.parameters.DsspLocalStructureLetter;
import meshi.parameters.LocalStructure;
import meshi.parameters.MeshiPotential;
import meshi.sequences.MeshiSequence;
import meshi.sequences.ResidueAlignmentMethod;
import meshi.sequences.SequenceAlignment;
import meshi.sequences.SequenceAlignmentColumn;
import meshi.util.*;
import meshi.util.dssp.DSSPFull;
import meshi.util.file.MeshiLineReader;
import meshi.util.file.MeshiWriter;
import meshi.util.filters.Filter;

import java.io.File;
import java.io.IOException;


public class LocalStructurePsiPredErrorAnalyzer extends MeshiProgram implements KeyWords {


    public static final String NAME = "Optimize";
    private static String inFileName, nativeFileName, outFileName,dsspFile;
    private static Protein model;
    private static CommandList commands;
    private static int seed;
    private static InflateCreator inflateByOtherModelCreator = new InflateCreator(InflateType.BY_OTHER_MODEL, new File("."), new MyFileFilter());
    private static String parentString;


    public static void main(String[] argv) throws IOException, OptimizerException, UpdateableException, EvaluationException {
        init(argv);
        ArcCos.useFastArcCos();
        model = null;

        parentString = getParentString(inFileName);
        model = Utils.getProtein(commands, inFileName, ResidueExtendedAtomsCreator.creator, Utils.defaultExceptionHandler);
        Utils.alignToX(model);
        addAtoms(model, commands);
        model.resetBonds();


        Utils.println("nAtoms " + model.atoms().size());
        Utils.println("nResidues " + model.residues().size());
        ((MyFileFilter) inflateByOtherModelCreator.filter).reset(model);

        if (!dsspFile.equals("NONE"))
            Utils.AssignFullDSSP(model, dsspFile);

        /****************************** PsiPred Error Analysis - Dssp specific string **************************/
        MeshiWriter result = new MeshiWriter(outFileName);
        MeshiSequence ssPrediction = Utils.getSSPrediction(commands.firstWord("predictedSsFileName").secondWord());
        String ssNativeDsspFileName = commands.firstWord("nativeDsspFileName").secondWord();
        DSSPFull nativeDssp = new DSSPFull(ssNativeDsspFileName);
        MeshiSequence nativeSsSequence = nativeDssp.getResidueSequenceWithFullSs()[0];

        SequenceAlignment ssAlignment = new SequenceAlignment();
        ssAlignment.addAll(SequenceAlignment.identityAlignment(ssPrediction, nativeSsSequence));

	LocalStructure dls = new LocalStructure("E____-A_");
        //LocalStructure dls = new LocalStructure("_____-__");
        //LocalStructure dls = new LocalStructure("T3__S+__");
        LocalStructureHCEprob probs = predictionSuccessDSSP3(ssAlignment,dls);

        
        //System.out.println(probs.sequence);
        //System.out.println(probs.sequenceWithDLS);
        //System.out.println("LocalStructreDsspLetter="+dls);
        //System.out.println("numberOfLetterApearancesInNative="+probs.numberOfLetterApearancesInNative);
        //System.out.println("Helix="+probs.helixProb/probs.numberOfLetterApearancesInNative+", isLetter="+probs.isLetter(DsspLocalStructureLetter.HELIX));
        //System.out.println("Sheet="+probs.sheetProb/probs.numberOfLetterApearancesInNative+", isLetter="+probs.isLetter(DsspLocalStructureLetter.SHEET));
        //System.out.println("Coil="+probs.coilProb/probs.numberOfLetterApearancesInNative+", isLetter="+probs.isLetter(DsspLocalStructureLetter.COIL));
        

        result.println(probs.sequence);
        result.println(probs.sequenceWithDLS);
        result.println("LocalStructreDsspLetter="+dls);
        /*
        result.println("numberOfLetterApearancesInNative="+probs.numberOfLetterApearancesInNative);
        result.println("Helix="+probs.helixProb/probs.numberOfLetterApearancesInNative+", isLetter="+probs.isLetter(DsspLocalStructureLetter.HELIX));
        result.println("Sheet="+probs.sheetProb/probs.numberOfLetterApearancesInNative+", isLetter="+probs.isLetter(DsspLocalStructureLetter.SHEET));
        result.println("Coil="+probs.coilProb/probs.numberOfLetterApearancesInNative+", isLetter="+probs.isLetter(DsspLocalStructureLetter.COIL));
        */
        result.println("numberOfLetterApearancesInNative="+probs.numberOfLetterApearancesInNative);
        result.println("Helix="+probs.helixProb);
        result.println("Sheet="+probs.sheetProb);
        result.println("Coil="+probs.coilProb);
        result.close();
    }
    private static class LocalStructureHCEprob{
        double helixProb = 0;
        double coilProb = 0;
        double sheetProb = 0;
        String sequence = "";
        String sequenceWithDLS = "";
        LocalStructure dls = null;
        double numberOfLetterApearancesInNative = 0;
        private static LocalStructureHCEprob probs = null;

        private LocalStructureHCEprob(LocalStructure dls){
            this.dls = dls;
        }
        protected static LocalStructureHCEprob getProbs(LocalStructure dls){
            if (probs != null){
                probs.helixProb = 0;
                probs.coilProb = 0;
                probs.sheetProb = 0;
                probs.sequence = "";
                probs.sequenceWithDLS = "";
                probs.dls = dls;
                probs.numberOfLetterApearancesInNative = 0;
                return probs;
            }
            probs = new LocalStructureHCEprob(dls);
            return probs;
        }

        protected boolean isLetter(DsspLocalStructureLetter other){
            return other.equals(dls.getDSSP3Reduction());
        }
    }

    private static LocalStructureHCEprob predictionSuccessDSSP3(SequenceAlignment alignment, LocalStructure dls){


        LocalStructureHCEprob probs = LocalStructureHCEprob.getProbs(dls);

        for (SequenceAlignmentColumn column : alignment){
            ResidueSsPrediction predictionResidue            = (ResidueSsPrediction)             column.cell0().getAttribute(MeshiAttribute.SECONDARY_STRUCTURE_ATTRIBUTE);
            LocalStructure nativeSsResidue         = (LocalStructure) column.cell1().getAttribute(MeshiAttribute.SS_DSSP);

            probs.sequence+=predictionResidue.getResidueType().nameOneLetter();
            if (predictionResidue != null && nativeSsResidue != null && dls.compareTo(nativeSsResidue)==0) {
                probs.numberOfLetterApearancesInNative+=1;

                probs.sequenceWithDLS+=predictionResidue.getResidueType().nameOneLetter();
                /*
                if (predictionResidue.secondaryStructure.equals(DsspLocalStructureLetter.COIL))
                    probs.coilProb += 1;
                else if (predictionResidue.secondaryStructure.equals(DsspLocalStructureLetter.HELIX))
                    probs.helixProb += 1;
                else if (predictionResidue.secondaryStructure.equals(DsspLocalStructureLetter.SHEET))
                    probs.sheetProb += 1;
                */
                probs.coilProb+= predictionResidue.getProbabilityOf(DsspLocalStructureLetter.COIL);
                probs.helixProb+= predictionResidue.getProbabilityOf(DsspLocalStructureLetter.HELIX);
                probs.sheetProb+= predictionResidue.getProbabilityOf(DsspLocalStructureLetter.SHEET);


            } else {
                probs.sequenceWithDLS+="_";
            }
        }

        return probs;
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

    private static void init(String[] argv) {
        //                 0            1            2            3                4              5
        String[] keys = {"commands", "inFileName", "dsspFile", "nativeFileName", "outFileName", "seed"};
        String[] arguments = getArguments(keys, argv);

        seed = Integer.parseInt(arguments[5]);
        System.out.println("seed " + seed);
        initRandom(seed);

        commands = new CommandList(arguments[0]);
        inFileName = arguments[1];
        dsspFile   = arguments[2];
        nativeFileName = arguments[3];
        outFileName = arguments[4];
        if (commands.keyExists("verbose")) Utils.verboseOn();
        else Utils.verboseOff();
    }

    private static String getParentString(String fileName) throws IOException {
        MeshiLineReader reader = new MeshiLineReader(fileName);
        String out;
        while ((out = reader.readLine()) != null) {
            if (out.startsWith("PARENT"))
                return out;
        }
        return "PARENT N/A";
    }

    private static class MyFileFilter implements Filter {
        private String prefix = null;
        Protein thisModel;

        public void reset(Protein thisModel) {
            this.thisModel = thisModel;
            int index = thisModel.name().indexOf('.');
            if (index != -1) prefix = thisModel.name().substring(0, index);
            else prefix = thisModel.name();

        }

        public boolean accept(Object obj) {
            if (prefix == null) throw new RuntimeException("prefix is " + prefix);
            File file = (File) obj;
            String path = file.getAbsolutePath();
            if (!path.endsWith("pdb")) return false;
            if (path.indexOf("out") == -1) return false;
            if (file.getName().startsWith(prefix)) {
                try {
                    double rms = Rms.rms(thisModel, Protein.getCAproteinFromApdbFile(file), ResidueAlignmentMethod.IDENTITY);
                    if (rms < 1) return false;
                } catch (Exception ex ) {
                    writeFailureXml(model, ex);
                    throw new RuntimeException(ex.getMessage());}
                return true;
            }
            return false;
        }
    }

    private static void writeFailureXml(Protein model, Exception exception)  {
        MeshiWriter writer;
        try {
            writer = new MeshiWriter(outFileName + ".xml");
        } catch (IOException ex) {throw  new RuntimeException("Cannot write failure XML file after exception:\n"+exception+"\n"+"Due to "+ex);}
        writer.println("<?xml version=\"1.0\" encoding=\"UTF-8\" ?> ");
        writer.println("<ProteinInfoList name=\"Failure report for Protein: " + model.sourceFile() + "\">");
        writer.print("<ProteinInfo  value=\"MCM_END\" step=\"0\" ");
//        for (Score scoreFunction : scoreFunctions) {
//            writer.print(scoreFunction.toString()+"_weightedMedianScore=\"0.05\" "+scoreFunction.toString()+"_interdecile=\"0\" ");
//        }
        writer.println("time=\"0\" fileName=\""+inFileName+"\" >");
        writer.println("<Exception>\n"+ exception + "\n</Exception>" );
        writer.println("</ProteinInfo>\n" + "</ProteinInfoList>\n");
        writer.close();
    }
}

