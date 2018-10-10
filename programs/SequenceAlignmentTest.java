package programs;

import meshi.sequences.MeshiSequence;
import meshi.sequences.ResidueAlignmentMethod;
import meshi.sequences.SequenceAlignment;
import meshi.sequences.aligner.IdentityMatrix;

import javax.sound.midi.Sequence;

public class SequenceAlignmentTest {
    public static void main(String[] args) {
        MeshiSequence sequence1 = new MeshiSequence(args[0],"s1");
        MeshiSequence sequence2 = new MeshiSequence(args[1], "s2");
        SequenceAlignment alignment = SequenceAlignment.substitutionAlignment(sequence1, sequence2, new IdentityMatrix(), ResidueAlignmentMethod.IDENTITY);
        System.out.println(alignment);

    }
}
