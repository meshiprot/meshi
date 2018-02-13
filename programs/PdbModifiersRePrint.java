package programs;

import meshi.PDB.PdbHeader;
import meshi.PDB.pdbFilters.PdbLineATOM;
import meshi.molecularElements.Chain;
import meshi.molecularElements.Protein;
import meshi.molecularElements.atoms.Atom;
import meshi.molecularElements.atoms.AtomList;
import meshi.molecularElements.extendedAtoms.ResidueExtendedAtomsCreator;
import meshi.sequences.ResidueAlignmentMethod;
import meshi.sequences.ResidueSequence;
import meshi.sequences.SequenceAlignment;
import meshi.sequences.aligner.IdentityMatrix;
import meshi.util.CommandList;
import meshi.util.MeshiProgram;
import meshi.util.file.MeshiLineReader;
import meshi.util.file.MeshiWriter;

import java.io.IOException;

/**
 * Created by siditom on 02/03/2017.
 */
public class PdbModifiersRePrint extends MeshiProgram {
    private static CommandList commands;
    private static String in3DFileName,inSeqFileName,outFileName;

    public static void main(String[] args) throws IOException {
        init(args);

        PdbHeader pdbHeader = new PdbHeader(in3DFileName);
        pdbHeader.updateSequenceWithModifiedResidues();
        AtomList al = new AtomList(in3DFileName, pdbHeader, new PdbLineATOM());

        //System.out.println(getSequenceFromFasta());
        ResidueSequence rsOriginalCASPSequence = new ResidueSequence(getSequenceFromFasta(),"Native sequence - CASP designated chain.");

        Protein protein = new Protein(al, ResidueExtendedAtomsCreator.creator);
        Chain maxScoreChain = null;
        double currentScore = Double.MIN_VALUE;
        for (Chain chain : protein.chains()) {
            ResidueSequence rsModifiedPDBbyChain = new ResidueSequence(chain.getNonDummyResidues(),"Native sequence - By chains of the native PDB. Chain "+chain.name());
            //ResidueSequence rsModifiedPDBbyChain = chain.sequence();
            SequenceAlignment s = SequenceAlignment.substitutionAlignment(rsOriginalCASPSequence,rsModifiedPDBbyChain, new IdentityMatrix(), ResidueAlignmentMethod.IDENTITY);
            if (s.score() > currentScore){
                currentScore = s.score();
                maxScoreChain = chain;
            }
        }

        MeshiWriter mw = new MeshiWriter(outFileName);

        //for (AtomList alChain :al.splitToChains()) {
            AtomList atoms = maxScoreChain.getNonDummyResidues().atoms();
            for (Atom a :atoms) {
                if (a.normal() && a.residue().number()>=0) mw.println(a.toString());
                //mw.println(a.toString());
            }
            mw.println("TER");
        //}

        mw.close();
    }
    private static void init(String[] argv) {
        //                 0            1            2            3                4              5
        String[] keys = {"commands", "in3DFileName","inSeqFileName", "outFileName"};
        String[] arguments = getArguments(keys, argv);



        commands = new CommandList(arguments[0]);
        in3DFileName = arguments[1];
        inSeqFileName = arguments[2];
        outFileName = arguments[3];
        //if (commands.keyExists("verbose")) Utils.verboseOn();
        //else Utils.verboseOff();
    }

    public static String getSequenceFromFasta(){
        String seq = null;
        MeshiLineReader mr = new MeshiLineReader(inSeqFileName);
        try {
            seq = "";
            String line = mr.readLine();
            if (line.lastIndexOf('>')>=0) line = mr.readLine();
            while (line != null){
		if (line.lastIndexOf('>')>=0) throw new RuntimeException("Error reading fasta file - multiple sequences in one file.");
                seq += line;
                line = mr.readLine();
            }
            mr.close();
        }catch(Exception e){
            System.out.println("ERROR - failed to read sequence file - "+inSeqFileName);
        }
        return seq;
    }

}
