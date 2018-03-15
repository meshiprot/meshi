package meshi.applications.optimize;

import meshi.energy.EnergyCreator;
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
import meshi.energy.pairwiseNonBondedTerms.atomicPairwisePMFSumma.AtomicPairwisePMFSummaCreator;
import meshi.energy.rapdf.RapdfCreator;
import meshi.energy.rg.RgCreator;
import meshi.energy.secondaryStructureAlphabet.DeepCNF_ssCompatibilityFeatureCreator;
import meshi.energy.simpleEnergyTerms.DistanceConstraints.DistanceConstraintsCreator;
import meshi.energy.simpleEnergyTerms.angle.AngleCreator;
import meshi.energy.simpleEnergyTerms.bond.BondCreator;
import meshi.energy.simpleEnergyTerms.compositeTorsions.compositePropensity.CompositePropensityCreator;
import meshi.energy.simpleEnergyTerms.compositeTorsions.compositePropensity.CooperativeZPropensityCreator;
import meshi.energy.simpleEnergyTerms.compositeTorsions.compositePropensity.CooperativeZStdPropensityCreator;
import meshi.energy.simpleEnergyTerms.compositeTorsions.ramachandran.RamachandranCreator;
import meshi.energy.simpleEnergyTerms.compositeTorsions.ramachandranCore.RamachandranCoreCreator;
import meshi.energy.simpleEnergyTerms.compositeTorsions.ramachandranSidechain.CooperativeZRamachandranCreator;
import meshi.energy.simpleEnergyTerms.compositeTorsions.ramachandranSidechain.CooperativeZStdRamachandranCreator;
import meshi.energy.simpleEnergyTerms.compositeTorsions.ramachandranSidechain.RamachandranSidechainEnergyCreator;
import meshi.energy.simpleEnergyTerms.inflate.inflate.InflateCreator;
import meshi.energy.simpleEnergyTerms.inflate.inflate.InflateType;
import meshi.energy.simpleEnergyTerms.outOfPlane.OutOfPlaneCreator;
import meshi.energy.simpleEnergyTerms.plane.PlaneCreator;
import meshi.energy.simpleEnergyTerms.tether.TetherCreator;
import meshi.energy.solvate.SolvateCreatorHBforMinimization;
import meshi.energy.solvation.AtomEnvironmentCreator;
import meshi.energy.solvation.AtomEnvironmentCreatorSS;
import meshi.energy.solvation.SolvationCreator;
import meshi.energy.twoTorsions.FlatRamachCreator;
import meshi.geometry.DistanceMatrix;
import meshi.molecularElements.Protein;
import meshi.optimizers.LBFGS;
import meshi.optimizers.MCM;
import meshi.optimizers.Relaxation;
import meshi.scoringFunctions.Score;
import meshi.util.CommandList;
import meshi.util.info.InfoType;
import programs.Optimize;

import java.io.File;
import java.util.ArrayList;

public interface OptimizeConstants {
    enum TETHER_FLAG {
        ON, OFF, RESET_ON, RESET_OFF
    }

    enum SUMMA_COOP_FLAG {
        ON, OFF
    }



    enum RAMACH_COOP_FLAG {
        ON, OFF
    }



    enum PROP_COOP_FLAG {
        ON, OFF
    }



    enum RG_FLAG {
        ON, OFF
    }

    enum SOLVATE_FLAG{
        ON, OFF
    }

    TetherCreator tetherAllCreator = new TetherCreator(InfoType.TETHER_ALL);
    TetherCreator residueTetherCreator = new TetherCreator(InfoType.TETHER);
    HydrogenBondsCreator hydrogenBondsCreator = new HydrogenBondsCreator();
    SolvateCreatorHBforMinimization solvateCreator = new SolvateCreatorHBforMinimization();
    SolvationCreator solvationCreator = new SolvationCreator();
    AtomEnvironmentCreator atomEnvironmentCreator = new AtomEnvironmentCreator(solvationCreator);
    AtomEnvironmentCreatorSS atomEnvironmentCreatorSS = new AtomEnvironmentCreatorSS(solvationCreator);
    RgCreator rgCreator = new RgCreator();
    //ContactsCreator contactsCreator = new ContactsCreator(InfoType.CONTACTS_ETC);
    RapdfCreator rapdfCreator = new RapdfCreator();

    CompositePropensityCreator propensityCreator = new CompositePropensityCreator();
    CooperativeZPropensityCreator cooperativeZPropensityCreator = new CooperativeZPropensityCreator(propensityCreator);
    CooperativeZStdPropensityCreator cooperativeZStdPropensityCreator = new CooperativeZStdPropensityCreator(cooperativeZPropensityCreator);

    AtomicPairwisePMFSummaCreator summaCreator = new AtomicPairwisePMFSummaCreator();
    CooperativeZ5typesSummaCreator cooperativeZSummaCreator = new CooperativeZ5typesSummaCreator(summaCreator);
    CooperativeZStd5typesSummaCreator cooperativeZStdSummaCreator = new CooperativeZStd5typesSummaCreator(cooperativeZSummaCreator);
    ExcludedVolCreator excludedVolCreator = new ExcludedVolCreator();
    RamachandranSidechainEnergyCreator ramachCreator = new RamachandranSidechainEnergyCreator();
    RamachandranCoreCreator ramachandranCoreCreator = new RamachandranCoreCreator(ramachCreator);
    CooperativeZRamachandranCreator cooperativeZRamachandranCreator = new CooperativeZRamachandranCreator(ramachCreator);
    CooperativeZStdRamachandranCreator cooperativeZStdRamachandranCreator = new CooperativeZStdRamachandranCreator(cooperativeZRamachandranCreator);
    RamachandranCreator ramachandranCreator = new RamachandranCreator();
    InflateCreator inflateCreator = new InflateCreator(InflateType.SIMPLE);
    InflateCreator inflatePerSegmentCreator = new InflateCreator(InflateType.PER_SEGMENT);
    InflateCreator inflateBySegmentCreator = new InflateCreator(InflateType.BY_SEGMENT);
    InflateCreator inflateByOtherModelCreator = new InflateCreator(InflateType.BY_OTHER_MODEL, new File("."), new MyFileFilter());
    ConservationContactsCreator conservationContactsCreator8 = new ConservationContactsCreator(InfoType.CONSERVATION_CONTACTS8);
    ConservationContactsCreator conservationContactsCreator11 = new ConservationContactsCreator(InfoType.CONSERVATION_CONTACTS11);
    ConservationContactsCreator conservationContactsCreator15 = new ConservationContactsCreator(InfoType.CONSERVATION_CONTACTS15);
    ConservationContactsHrCreator conservationContactsHrCreator = new ConservationContactsHrCreator();
    HydrogenBondsPairsCreator hydrogenBondsPairsCreator = new HydrogenBondsPairsCreator(hydrogenBondsCreator);
    FlatRamachCreator flatRamachCreator = new FlatRamachCreator();    //ExcludedVolCreator excludedVolCreator = new ExcludedVolCreator();
    //OneCreator oneCreator = new OneCreator();
    PlaneCreator planeCreator = new PlaneCreator();
    GoapCreator goapCreator = new GoapCreator();
    final ContactsCreator contactsCreator  = new ContactsCreator();

    EnergyCreator[] energyCreators = {
            new DistanceConstraintsCreator(),
            new BondCreator(),
            new AngleCreator(),
            planeCreator,
            new OutOfPlaneCreator(),
            //excludedVolCreator,
            ramachandranCreator,
            ramachCreator,
            cooperativeZRamachandranCreator,
            cooperativeZStdRamachandranCreator,
            ramachandranCoreCreator,
            propensityCreator,
            cooperativeZPropensityCreator,
            cooperativeZStdPropensityCreator,
            summaCreator,
            excludedVolCreator,
            cooperativeZSummaCreator,
            cooperativeZStdSummaCreator,
            solvationCreator,
            atomEnvironmentCreator,
            atomEnvironmentCreatorSS,
            solvateCreator,
            hydrogenBondsCreator,
            hydrogenBondsPairsCreator,
            new HbondsPunishHOCAngleCreator(hydrogenBondsCreator),
            new HBondsPunishOHNAngleCreator(hydrogenBondsCreator),
            contactsCreator,
            rgCreator,
            goapCreator,
            conservationContactsCreator8,
            conservationContactsCreator11,
            conservationContactsCreator15,
            conservationContactsHrCreator,
            //oneCreator,
            inflateCreator,
            inflatePerSegmentCreator,
            inflateBySegmentCreator,
            inflateByOtherModelCreator,
            residueTetherCreator,
            flatRamachCreator,
            tetherAllCreator,
            rapdfCreator,
            new StaticFeaturesCreator(),
            //new SecondaryStructureAlphabetFeatures_plusNoiseCreator(),
            new DeepCNF_ssCompatibilityFeatureCreator()};
    EnergyCreator[] excludeFromMinimization = {inflateCreator, inflatePerSegmentCreator,
            inflateBySegmentCreator, inflateByOtherModelCreator, ramachandranCreator};       // inflate is noise  and the other two have almost arbitrary values.
    EnergyCreator[] excludeFromPerturbation1 = {tetherAllCreator, residueTetherCreator, inflatePerSegmentCreator, inflateBySegmentCreator,
            inflateByOtherModelCreator, rapdfCreator,ramachandranCoreCreator};
    EnergyCreator[] excludeFromPerturbation2 = {tetherAllCreator, residueTetherCreator, inflateCreator, inflateBySegmentCreator,
            inflateByOtherModelCreator , rapdfCreator, ramachandranCoreCreator};
    EnergyCreator[] excludeFromPerturbation3 = {tetherAllCreator, residueTetherCreator, inflateCreator, inflatePerSegmentCreator,
            inflateByOtherModelCreator , rapdfCreator, ramachandranCoreCreator};
    EnergyCreator[] excludeFromPerturbation4 = {tetherAllCreator, residueTetherCreator, inflateCreator, inflatePerSegmentCreator,
            inflateBySegmentCreator, rapdfCreator, ramachandranCoreCreator};
    EnergyCreator[] relax1Creators = {new BondCreator(), new AngleCreator(), tetherAllCreator, residueTetherCreator};

    
}
