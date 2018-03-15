package programs;

import meshi.energy.EvaluationException;
import meshi.energy.simpleEnergyTerms.angle.AngleParametersList;
import meshi.energy.simpleEnergyTerms.bond.BondParametersList;
import meshi.energy.simpleEnergyTerms.plane.PlaneParametersList;
import meshi.geometry.putH.PutHpos;
import meshi.geometry.putH.PutHposLog;
import meshi.molecularElements.Protein;
import meshi.molecularElements.atoms.Atom;
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

import java.io.File;
import java.io.IOException;
import java.io.RandomAccessFile;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.nio.file.StandardCopyOption;

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
    private static String outPath;
    private static String modelFilePath,dsspFilePath,scwrlPdbFilePath;

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
            modelFilePath=tmp+"MODEL."+iMGroup+"." + i;
            dsspFilePath= modelFilePath+".scwrl.pdb.dssp";
            scwrlPdbFilePath=modelFilePath+".scwrl.pdb";
	        new File(tmp).mkdirs();

            //Step1.a - generate a single pdb file for the current model

            seperateModel(inFileName,loc,modelFilePath+".pdb");

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
            //String[] keys = {"commands", "inFileName", "dsspFile", "nativeFileName", "outFileName", "seed"};
            outFileName=modelFilePath+".out.pdb";
            Optimize.main(new String[]{argv[0],"-inFileName="+modelFilePath+".pdb","-dsspFile="+dsspFilePath,"-nativeFileName=NONE","-outFileName="+outFileName,"-seed="+seed});
            //Step3 - copy the meshi result files - pdb and xml - to the out directory.
            new File(outPath).mkdirs();
            Files.move(Paths.get(outFileName),Paths.get(outPath+Paths.get(outFileName).getFileName().toString()), StandardCopyOption.REPLACE_EXISTING);
            Files.move(Paths.get(outFileName+".xml"),Paths.get(outPath+Paths.get(outFileName+".xml").getFileName().toString()), StandardCopyOption.REPLACE_EXISTING);

            //Step4 - Delete the tmp folder (dssp, and scwrl file will be lost).
            Files.delete(Paths.get(tmpPath));
        }

    }

    private static void init(String[] argv) {
        //                 0            1            2            3                4              5
        String[] keys = {"commands", "inFileName","iGroup","seed","out"};
        String[] arguments = getArguments(keys, argv);

        commands = new CommandList(arguments[0]);
        inFileName = arguments[1];
        iMGroup = Integer.parseInt(arguments[2]);
        loc = 0; //TODO
        nModels = 1;//TODO
        seed = Integer.parseInt(arguments[3]);
        outPath=arguments[4];

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

    private static void seperateModel(String filepath,Long fLoc,String outPdb) {
        try {
            MeshiWriter mw = new MeshiWriter(outPdb);
            RandomAccessFile raf = new RandomAccessFile(filepath, "r");//Open our file with read/write access

            raf.seek(fLoc);
            String line = "";
            line = raf.readLine();
            while (line != null && !line.contains("END MODEL")) {
                mw.println(line);
                line = raf.readLine();
            }
            if (line !=null) mw.println(line);

            raf.close();
            mw.close();
        } catch (Exception e) {
            e.printStackTrace();
        }
    }
}
