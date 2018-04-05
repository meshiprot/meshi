/*
 * Copyright (c) 2010. of Chen Keasar, BGU . For free use under LGPL
 */

package programs;

import meshi.applications.optimize.OptimizeEnergies;
import meshi.applications.optimize.OptimizeFiles;
import meshi.applications.optimize.OptimizeUtils;
import meshi.energy.EnergyCreator;
import meshi.energy.EvaluationException;
import meshi.energy.TotalEnergy;
import meshi.energy.compatebility.StaticFeaturesCreator;
import meshi.energy.conservationContacts.ConservationContactsCreator;
import meshi.energy.conservationContacts.ConservationContactsHrCreator;
import meshi.energy.contacts.ContactsCreator;
import meshi.energy.goap.GoapCreator;
import meshi.energy.hydrogenBond.HydrogenBondsCreator;
import meshi.energy.hydrogenBondsAngle.HBondsPunishOHNAngleCreator;
import meshi.energy.hydrogenBondsAngle.HbondsPunishHOCAngleCreator;
import meshi.energy.hydrogenBondsPairs.HydrogenBondsPairsCreator;
import meshi.energy.pairwiseNonBondedTerms.CooperativeSumma.CooperativeZ5typesSummaCreator;
import meshi.energy.pairwiseNonBondedTerms.CooperativeSumma.CooperativeZStd5typesSummaCreator;
import meshi.energy.pairwiseNonBondedTerms.ExcludedVol.ExcludedVolCreator;
import meshi.energy.pairwiseNonBondedTerms.atomicPairwisePMFSumma.AtomicPairwisePMFSumma;
import meshi.energy.pairwiseNonBondedTerms.atomicPairwisePMFSumma.AtomicPairwisePMFSummaCreator;
import meshi.energy.rapdf.RapdfCreator;
import meshi.energy.rg.RgCreator;
import meshi.energy.secondaryStructureAlphabet.DeepCNF_ssCompatibilityFeatureCreator;
import meshi.energy.simpleEnergyTerms.DistanceConstraints.DistanceConstraintsCreator;
import meshi.energy.simpleEnergyTerms.angle.AngleCreator;
import meshi.energy.simpleEnergyTerms.angle.AngleParametersList;
import meshi.energy.simpleEnergyTerms.bond.BondCreator;
import meshi.energy.simpleEnergyTerms.bond.BondParametersList;
import meshi.energy.simpleEnergyTerms.compositeTorsions.ResidueTorsions;
import meshi.energy.simpleEnergyTerms.compositeTorsions.ResidueTorsionsList;
import meshi.energy.simpleEnergyTerms.compositeTorsions.compositePropensity.CompositePropensityCreator;
import meshi.energy.simpleEnergyTerms.compositeTorsions.compositePropensity.CooperativeZPropensityCreator;
import meshi.energy.simpleEnergyTerms.compositeTorsions.compositePropensity.CooperativeZStdPropensityCreator;
import meshi.energy.simpleEnergyTerms.compositeTorsions.ramachandran.RamachandranCreator;
import meshi.energy.simpleEnergyTerms.compositeTorsions.ramachandran.RamachandranEnergy;
import meshi.energy.simpleEnergyTerms.compositeTorsions.ramachandran.RamachandranEnergyElement;
import meshi.energy.simpleEnergyTerms.compositeTorsions.ramachandranCore.RamachandranCoreCreator;
import meshi.energy.simpleEnergyTerms.compositeTorsions.ramachandranSidechain.CooperativeZRamachandranCreator;
import meshi.energy.simpleEnergyTerms.compositeTorsions.ramachandranSidechain.CooperativeZStdRamachandranCreator;
import meshi.energy.simpleEnergyTerms.compositeTorsions.ramachandranSidechain.RamachandranSidechainEnergyCreator;
import meshi.energy.simpleEnergyTerms.inflate.inflate.InflateCreator;
import meshi.energy.simpleEnergyTerms.inflate.inflate.InflateType;
import meshi.energy.simpleEnergyTerms.outOfPlane.OutOfPlaneCreator;
import meshi.energy.simpleEnergyTerms.plane.PlaneCreator;
import meshi.energy.simpleEnergyTerms.plane.PlaneEnergy;
import meshi.energy.simpleEnergyTerms.plane.PlaneParametersList;
import meshi.energy.simpleEnergyTerms.tether.TetherCreator;
import meshi.energy.simpleEnergyTerms.tether.TetherEnergy;
import meshi.energy.solvate.SolvateCreatorHBforMinimization;
import meshi.energy.solvation.AtomEnvironmentCreator;
import meshi.energy.solvation.AtomEnvironmentCreatorSS;
import meshi.energy.solvation.SolvationCreator;
import meshi.energy.twoTorsions.FlatRamachCreator;
import meshi.geometry.ArcCos;
import meshi.geometry.DistanceMatrix;
import meshi.geometry.putH.PutHpos;
import meshi.geometry.putH.PutHposLog;
import meshi.molecularElements.Protein;
import meshi.molecularElements.Residue;
import meshi.molecularElements.ResidueList;
import meshi.molecularElements.atoms.Atom;
import meshi.molecularElements.extendedAtoms.ResidueExtendedAtoms;
import meshi.molecularElements.extendedAtoms.ResidueExtendedAtomsCreator;
import meshi.molecularElements.loops.AtomFinding;
import meshi.optimizers.*;
import meshi.parameters.MeshiPotential;
import meshi.parameters.SecondaryStructure;
import meshi.scoringFunctions.CombinedEnergyScore;
import meshi.scoringFunctions.GdtScore;
import meshi.scoringFunctions.Score;
import meshi.sequences.AlignmentException;
import meshi.sequences.ResidueAlignment;
import meshi.sequences.ResidueAlignmentMethod;
import meshi.util.*;
import meshi.util.file.MeshiLineReader;
import meshi.util.file.MeshiWriter;
import meshi.util.filters.Filter;
import meshi.util.info.*;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Date;

//import meshi.energy.one.OneCreator;

public class EvaluateModel extends MeshiProgram implements KeyWords {
    private static CommandList commands = null;
    private static OptimizeFiles files;


    public static void main(String[] argv) throws IOException, OptimizerException, UpdateableException, EvaluationException, AlignmentException {
        init(argv);
        Protein model = OptimizeUtils.getModel(commands,files);
        OptimizeEnergies energies = new OptimizeEnergies(model, commands);
        TotalEnergy energy = energies.minimizationEnergy;
        energy.setCurrent();
        ArrayList<Score> scoreFunctions = energies.scoreFunctions;
        Protein nativeStructure;
            if (!(files.nativeStructure.getName().equals("NONE") || files.nativeStructure.getName().equals("none"))) {
                nativeStructure = Utils.getProtein(commands, files.nativeStructure.getAbsolutePath(), ResidueExtendedAtoms.creator, Utils.defaultExceptionHandler);
                if (nativeStructure.chains().size() > 1)  throw new RuntimeException("Current version does not allow complexes in the native structure");
            } else {
                nativeStructure = null;
            }

            ModelAnalyzer modelAnalyzer = new ModelAnalyzer(model, nativeStructure, null, energy, scoreFunctions, ResidueAlignmentMethod.IDENTITY);
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
