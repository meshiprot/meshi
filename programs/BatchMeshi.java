package programs;

import meshi.PDB.pdbFilters.PdbLineMultipleModelsFilter;
import meshi.energy.EvaluationException;
import meshi.energy.simpleEnergyTerms.angle.AngleParametersList;
import meshi.energy.simpleEnergyTerms.bond.BondParametersList;
import meshi.energy.simpleEnergyTerms.plane.PlaneParametersList;
import meshi.geometry.putH.PutHpos;
import meshi.geometry.putH.PutHposLog;
import meshi.molecularElements.Protein;
import meshi.molecularElements.atoms.Atom;
import meshi.molecularElements.extendedAtoms.ResidueExtendedAtomsCreator;
import meshi.molecularElements.loops.AtomFinding;
import meshi.optimizers.OptimizerException;
import meshi.parameters.MeshiPotential;
import meshi.sequences.AlignmentException;
import meshi.util.*;
import meshi.util.externalProgExec.DsspExec;
import meshi.util.externalProgExec.ExternalProgExecutioner;
import meshi.util.externalProgExec.ScwrlExec;
import meshi.util.file.MeshiLineReader;
import meshi.util.file.MeshiWriter;

import java.nio.file.Paths;

import java.io.File;
import java.io.IOException;

/**
 * Created by siditom on 08/03/2018.
 * This program is designed to get a Batch PDB file - a pdb file woth multiple models, and an indexing file - with the start position of each model.
 * For each model it generates the model's DSSP,SCWRL files.
 * Then it runs the meshi optimization program.
 * @Return XML Meshi Feature files, Meshi PDB file.
 */
public class BatchMeshi extends MeshiProgram implements KeyWords {

    public static final String NAME = "BatchMeshi";
    private static int seed;
    private static String tmpPath="/tmp/";
    private static String inFileName, nativeFileName, outFileName,dsspFile;
    private static int nModels,iMGroup;
    private static long loc;
    private static Protein model, originalModel;
    private static CommandList commands;
    private static String parentString;

    public static void main(String[] argv) throws IOException, OptimizerException, UpdateableException, EvaluationException, AlignmentException {
        init(argv);
        model = null;

        //TODO - add options: (1-save dssp and scwrl files)
        //TODO - add options: (2-ouput only meshi-score, without any files)
        //TODO - add options: (3-output a single line of features, seperated by ',')
        //TODO - add options: (4-incorparate rosetta and voronoi score in the final xml).
        //TODO - add options: (5-incorparate rosetta and voronoi score in the meshi optimization process).
        ExternalProgExecutioner dsspExec = new DsspExec(commands);
        ExternalProgExecutioner scwrlExec = new ScwrlExec(commands);
        for (int i=0;i<nModels; i++) {

            //Step1 - Create tmp folder
            String tmp = tmpPath+tmpPath+"Meshi." + iMGroup + "." + Paths.get(inFileName).getFileName().toString()+File.separator;
            String modelFilePath=tmp+"MODEL."+iMGroup+"." + i;
	        new File(tmp).mkdirs();

            //Step1.a - generate a single pdb file for the current model
            ProteinGenerator pg = new ProteinGenerator(inFileName, loc, nModels, new PdbLineMultipleModelsFilter());
            model = pg.getProtein(commands, ResidueExtendedAtomsCreator.creator, Utils.defaultExceptionHandler);
            Utils.alignToX(model);
            addAtoms(model, commands);
            model.resetBonds();

            MeshiWriter mw = new MeshiWriter(modelFilePath+".pdb");
            mw.print(pg.getPdbAsString());
            mw.close();
            //model.printAtomsToFile(modelFilePath+".pdb");

            //Step1.b - generate scwrl4 files - pdb and log (scwrl score).
            //Step1.b - generate the model's dssp file.
            try {
                scwrlExec.exec(new String[]{modelFilePath+".pdb",modelFilePath+".scwrl.pdb"});
                dsspExec.exec(new String[]{modelFilePath+".scwrl.pdb",modelFilePath+".scwrl.pdb.dssp"});
            }catch(Exception e){
                System.err.println(e);
                System.err.println(e.getStackTrace());
            }

            //Step2 - activate Meshi Optimize with the generated scwrl4 and dssp files.
            //Step3 - copy the meshi result files - pdb and xml - to the out directory.
            //Step4 - Delete the tmp folder (dssp, and scwrl file will be lost).
        }

    }

    private static void init(String[] argv) {
        //                 0            1            2            3                4              5
        String[] keys = {"commands", "inFileName","iGroup","seed"};
        String[] arguments = getArguments(keys, argv);

        commands = new CommandList(arguments[0]);
        inFileName = arguments[1];
        iMGroup = Integer.parseInt(arguments[2]);
        loc = 0; //TODO
        nModels = 1;//TODO
        seed = Integer.parseInt(arguments[3]);
        System.out.println("seed " + seed);
        initRandom(seed);


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
}
