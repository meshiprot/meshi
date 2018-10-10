package meshi.util;


import meshi.PDB.pdbLines.PdbATOMLine;
import meshi.PDB.pdbLines.PdbLine;
import meshi.energy.compatebility.ResidueSsPrediction;
import meshi.molecularElements.Protein;
import meshi.molecularElements.Residue;
import meshi.molecularElements.atoms.Atom;
import meshi.sequences.SequenceAlignment;
import meshi.sequences.SequenceAlignmentColumn;
import meshi.util.file.MeshiWriter;

import javax.swing.*;

public class PdbAlignmentWriter extends PdbWriter {
    SequenceAlignment sequenceAlignment;

    public PdbAlignmentWriter(MeshiWriter writer, Protein protein, SequenceAlignment sequenceAlignment){
        super(writer, protein);
        this.sequenceAlignment = sequenceAlignment;
    }

    public void print() {
        for (SequenceAlignmentColumn column : sequenceAlignment) {
            ResidueSsPrediction prediction         = (ResidueSsPrediction) column.cell0().getAttribute(MeshiAttribute.SECONDARY_STRUCTURE_ATTRIBUTE);
            Residue residue                        = (Residue)             column.cell1().getAttribute(MeshiAttribute.RESIDUE_ATTRIBUTE);
            if (residue != null) {
                for (Atom atom : residue.getAtoms()) {
                    String line = atom.toString();
                    String numberS = "     "+prediction.getNumber();
                    numberS = numberS.substring(numberS.length() - 4);
                    String newLine = line.substring(0,22)+numberS+line.substring(26);
                    writer.println(newLine);
                }
            }

        }
    }
}
