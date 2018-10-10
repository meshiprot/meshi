/*
 * Copyright (c) 2010. of Chen Keasar, BGU . For free use under LGPL
 */

package programs;

import meshi.applications.optimize.OptimizeEnergies;
import meshi.applications.optimize.OptimizeFiles;
import meshi.applications.optimize.OptimizeUtils;
import meshi.energy.EvaluationException;
import meshi.energy.TotalEnergy;
import meshi.geometry.ArcCos;
import meshi.molecularElements.Protein;
import meshi.molecularElements.extendedAtoms.ResidueExtendedAtoms;
import meshi.molecularElements.extendedAtoms.ResidueExtendedAtomsCreator;
import meshi.optimizers.*;
import meshi.scoringFunctions.Score;
import meshi.sequences.AlignmentException;
import meshi.sequences.ResidueAlignmentMethod;
import meshi.util.*;
import meshi.util.info.*;

import java.io.IOException;
import java.util.ArrayList;

//import meshi.energy.one.OneCreator;

public class EvaluateModel extends MeshiProgram implements KeyWords {
    private static CommandList commands = null;
    private static OptimizeFiles files;


    public static void main(String[] argv) throws IOException, OptimizerException, UpdateableException, EvaluationException, AlignmentException {
        init(argv);
        Protein model = OptimizeUtils.getModel(commands,files);
        Protein originalModel = Utils.getProtein(commands, files.modelIn.getAbsolutePath(), ResidueExtendedAtomsCreator.creator, Utils.defaultExceptionHandler);
        OptimizeEnergies energies = new OptimizeEnergies(model, commands);
        TotalEnergy energy = energies.minimizationEnergy;
        LBFGS lbfgs = new LBFGS(energy, 0.01, 3000, 200);
        //lbfgs.run();
        energy.setCurrent();
        ArrayList<Score> scoreFunctions = energies.scoreFunctions;
        Protein nativeStructure;
            if (!(files.nativeStructure.getName().equals("NONE") || files.nativeStructure.getName().equals("none"))) {
                nativeStructure = Utils.getProtein(commands, files.nativeStructure.getAbsolutePath(), ResidueExtendedAtoms.creator, Utils.defaultExceptionHandler);
                if (nativeStructure.chains().size() > 1)  throw new RuntimeException("Current version does not allow complexes in the native structure");
            } else {
                nativeStructure = null;
            }

            ModelAnalyzer modelAnalyzer = new ModelAnalyzer(model, nativeStructure, originalModel, energy, scoreFunctions, ResidueAlignmentMethod.IDENTITY);
            ChainsInfo chainsInfo = new ChainsInfo(model);
            ProteinInfoOLd proteinInfo = modelAnalyzer.analyze("Evaluate "+model.name(), chainsInfo);
            System.out.println(proteinInfo.toXML(4));
    }

    private static void init(String[] argv) throws IOException{
        //                 0            1            2            3                4              5
        String[] keys = {"commands", "inFileName", "dsspFile", "nativeFileName", "outFileName", "seed"};
        String[] arguments = OptimizeUtils.getArguments(keys, argv);

        int seed = Integer.parseInt(arguments[5]);
        try {
            initRandom(seed);
            System.out.println("seed " + seed);
        } catch (RuntimeException e){
            Utils.println(e.getMessage() + "\n seed wasn't updated.\n Continue...");
        }

        commands = new CommandList(arguments[0]);
        files = new OptimizeFiles(arguments[1], arguments[2], arguments[3],arguments[4]);//inFileName dsspFile  nativeFileName outFileName
        if (commands.keyExists("verbose")) Utils.verboseOn();
        else Utils.verboseOff();
        ArcCos.useFastArcCos();
    }

}
