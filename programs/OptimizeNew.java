package programs;

import meshi.applications.optimize.*;
import meshi.energy.EvaluationException;
import meshi.energy.TotalEnergy;
import meshi.energy.simpleEnergyTerms.tether.TetherEnergy;
import meshi.geometry.ArcCos;
import meshi.molecularElements.Protein;
import meshi.molecularElements.atoms.Atom;
import meshi.molecularElements.extendedAtoms.ResidueExtendedAtomsCreator;
import meshi.optimizers.LBFGS;
import meshi.optimizers.MCM;
import meshi.optimizers.OptimizerException;
import meshi.optimizers.SteepestDecent;
import meshi.sequences.AlignmentException;
import meshi.util.CommandList;
import meshi.util.MeshiProgram;
import meshi.util.UpdateableException;
import meshi.util.Utils;
import meshi.util.info.ChainsInfo;

import java.io.IOException;

import static meshi.util.KeyWords.RELAX;


public class OptimizeNew extends MeshiProgram implements OptimizeConstants{
    private static CommandList commands = null;
    private static OptimizeFiles files;
    private static String parentString = "N/A";

//    private static MCM mcm = null;
//    private static Relaxation relaxation = null;
//    private static LBFGS lbfgs = null;
//    private static String inFileName = null, nativeFileName = null, outFileName = null,dsspFile = null;
//    private static Protein model = null, originalModel = null;
//    private static Boolean OK = null;
//        private static TotalEnergy minimizationEnergy = null;
//    private static ArrayList<Score> scoreFunctions = null;
//    private static DistanceMatrix distanceMatrix = null;
//    private static OptimizeLogger log = null;

    public static void main(String[] argv) throws IOException, OptimizerException, UpdateableException, EvaluationException, AlignmentException {
        init(argv);
        Protein model = OptimizeUtils.getModel(commands,files);
        Protein originalModel = Utils.getProtein(commands, files.modelIn.getAbsolutePath(), ResidueExtendedAtomsCreator.creator, Utils.defaultExceptionHandler);
        run(model, originalModel, commands, files);
    }

    private static void init(String[] argv) throws IOException{
        //                 0            1            2            3                4              5
        String[] keys = {"commands", "inFileName", "dsspFile", "nativeFileName", "outFileName", "seed"};
        String[] arguments = OptimizeUtils.getArguments(keys, argv);

        int seed = Integer.parseInt(arguments[5]);
        try {
            System.out.println("seed " + seed);
            initRandom(seed);
        } catch (RuntimeException e){
            Utils.println(e.getMessage() + "\n seed wasn't updated.\n Continue...");
        }

        commands = new CommandList(arguments[0]);
        files = new OptimizeFiles(arguments[1], arguments[2], arguments[3],arguments[4]);//inFileName dsspFile  nativeFileName outFileName
        if (commands.keyExists("verbose")) Utils.verboseOn();
        else Utils.verboseOff();
        ArcCos.useFastArcCos();
        parentString = OptimizeUtils.getParentString(files.modelIn.getAbsolutePath());
    }

    public static void run(Protein model, Protein originalModel, CommandList commands, OptimizeFiles files) {
        firstRelax(model);
        OptimizeEnergies energies = new OptimizeEnergies(model, commands);
        energies.minimizationEnergy.setCurrent();
        LBFGS lbfgs = Utils.getLBFGS(energies.minimizationEnergy, commands, RELAX);
        OptimizeLogger log = new OptimizeLogger(model, files, energies.minimizationEnergy, parentString, commands);
        ChainsInfo chainsInfo = new ChainsInfo(model);
        try {
            log.mcm(energies.scoreFunctions, energies.minimizationEnergy, "BEGINNING", null);  //the null is for residue info list
            OptimizeUtils.minimizeCassette(model, files, log, lbfgs, energies.minimizationEnergy, (TetherEnergy) tetherAllCreator.term());

            log.setEnergy(energies.minimizationEnergy);
            MCM mcm = OptimizeUtils.getMCM(model, energies.minimizationEnergy, energies.scoreFunctions, energies.perturbationEnergy1,
                    energies.perturbationEnergy2, energies.perturbationEnergy3, energies.perturbationEnergy3, commands);
            mcm.run(log);
            log.mcm(energies.scoreFunctions, energies.minimizationEnergy, "MCM_END\" step=\"" + mcm.maxSteps, chainsInfo);
        } catch (RuntimeException ex) {
            OptimizeUtils.writeFailureXml(model,ex, files);
            ex.printStackTrace();
            throw ex;
        }
    }

    private static void firstRelax(Protein model) {
        TotalEnergy firstRelaxEnergy = new TotalEnergy(model, relax1Creators, commands,"relaxEnergy");
        try {
            SteepestDecent steepestDecent = new SteepestDecent(firstRelaxEnergy,1,1000,200,0.00000001,0.5,1.5);
            steepestDecent.run();
        } catch (Exception ex) {
            throw new RuntimeException(ex);
        }
        Utils.println("nAtoms " + model.atoms().size());
        Utils.println("nResidues " + model.residues().size());
        ((MyFileFilter) inflateByOtherModelCreator.filter).reset(model);
        model.atoms().defrost();
        for (Atom atom : model.atoms()) {
            if (!atom.nowhere()) {
                for (Atom neighbor : atom.bonded()) {
                    if (neighbor.nowhere()) atom.freeze();
                }
            }
        }
    }


}
