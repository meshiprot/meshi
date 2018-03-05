package programs;

import meshi.PDB.pdbFilters.PdbLineMultipleModelsFilter;
import meshi.energy.EvaluationException;
import meshi.energy.simpleEnergyTerms.angle.AngleParametersList;
import meshi.energy.simpleEnergyTerms.bond.BondParametersList;
import meshi.energy.simpleEnergyTerms.plane.PlaneParametersList;
import meshi.geometry.ArcCos;
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
import meshi.util.externalProgExec.ExternalFeatureExtractor;
import meshi.util.file.MeshiLineReader;

import java.io.IOException;

/**
 * Created by user on 04/03/2018.
 */
public class ProteinGeneratorTest extends MeshiProgram implements KeyWords {

    public static final String NAME = "ProteinGeneratorTest";
    private static int seed;
    private static String inFileName, nativeFileName, outFileName,dsspFile;
    private static int nModels;
    private static long loc;
    private static Protein model, originalModel;
    private static CommandList commands;
    private static String parentString;

    public static void main(String[] argv) throws IOException, OptimizerException, UpdateableException, EvaluationException, AlignmentException {
        init(argv);
        ArcCos.useFastArcCos();
        boolean stepwiseFlag = commands.keyExists("stepwise");
        // Chemistry
        model = null;

        parentString = getParentString(inFileName);
        ProteinGenerator pg = new ProteinGenerator(inFileName, loc, nModels, new PdbLineMultipleModelsFilter());
        for (int i=0;i<nModels; i++) {
            model = pg.getProtein(commands, ResidueExtendedAtomsCreator.creator, Utils.defaultExceptionHandler);
            Utils.alignToX(model);
            addAtoms(model, commands);
            model.resetBonds();
            model.printAtomsToFile("test"+i+".pdb");
            model = null;
            try {
                ExternalFeatureExtractor.getExternals(commands, pg.getPdbAsString());
            }catch(Exception e){
                System.err.println(e);
                System.err.println(e.getStackTrace());
            }
        }





    }

    private static void init(String[] argv) {
        //                 0            1            2            3                4              5
        String[] keys = {"commands", "inFileName","fLoc","nModels","seed"};
        String[] arguments = getArguments(keys, argv);

        commands = new CommandList(arguments[0]);
        inFileName = arguments[1];

        loc = Long.parseLong(arguments[2]);
        nModels = Integer.parseInt(arguments[3]);
        seed = Integer.parseInt(arguments[4]);
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
