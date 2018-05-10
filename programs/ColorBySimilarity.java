package programs;

import meshi.molecularElements.Protein;
import meshi.molecularElements.Residue;
import meshi.molecularElements.atoms.Atom;
import meshi.molecularElements.extendedAtoms.ResidueExtendedAtoms;
import meshi.sequences.ResidueAlignment;
import meshi.sequences.ResidueAlignmentColumn;
import meshi.util.*;
import meshi.util.file.MeshiLineReader;
import meshi.util.file.MeshiWriter;
import meshi.util.info.ChainInfo;
import meshi.util.info.DoubleInfoElement;
import meshi.util.info.InfoType;
import meshi.util.info.ResidueInfo;

import java.io.File;
import java.util.ArrayList;

public class ColorBySimilarity {
    public static void main(String[] args) {
        CommandList commands = new CommandList(args[0]);
        File inDir = new File(args[1]);
        File outDir = new File(args[2]);
        ArrayList<DecoyAndMeta> decoys = getProteins(inDir, commands);
        for (DecoyAndMeta decoyAndMeta : decoys) {
            color(decoyAndMeta.decoy, decoys);
            writeDecoy(decoyAndMeta, outDir);
        }
    }

    public static void writeDecoy(DecoyAndMeta decoyAndMeta, File outDir) {
        if (!outDir.exists())
            outDir.mkdir();
        try {
            MeshiWriter writer = new MeshiWriter(outDir.getAbsolutePath() + "/" + decoyAndMeta.decoy.name());
            for (String s : decoyAndMeta.meta)
                if ((!s.startsWith("BEGIN")) & (!s.startsWith("MCM")) & (!s.startsWith("name")) )
                    writer.println(s);
            decoyAndMeta.decoy.atoms().print(writer);
            writer.println("TER\nEND");
            writer.close();
        } catch (Exception ex) {throw new RuntimeException(ex);}
    }
    public static void color(Protein decoy, ArrayList<DecoyAndMeta> decoys) {
        double[] displacements = new double[decoy.residues().get(decoy.residues().size()-1).number()+1];
        int sum = 0;
        for (DecoyAndMeta other : decoys)  {
            //Utils.printDebug(" colorBySimilarity ", decoy.name() +" "+other.decoy.name());
            if (decoy == other.decoy)
                continue;
            ResidueAlignment residueAlignment = ModelAnalyzer.getResidueAlignment(decoy, other.decoy);
            double[][] decoyCoor = new double[3][residueAlignment.size()];
            double[][] otherCoor  = new double[3][residueAlignment.size()];
            double[][] newCoor    = new double[3][residueAlignment.size()]; // Model CA coordinates after overlap transformation
            ModelAnalyzer.overlap(residueAlignment, decoyCoor, otherCoor, newCoor, 0.8);
            for (int i = 0; i < residueAlignment.size(); i++) {
                ResidueAlignmentColumn column = residueAlignment.get(i);
                Residue decoyResidue = column.residue0();
                displacements[decoyResidue.number()] += ModelAnalyzer.dis(decoyCoor, newCoor, i);
            }
            sum++;
        }
        for (int i = 0; i < displacements.length; i++)
            displacements[i] /= sum;

         for (Residue residue : decoy.residues()) {
                for (Atom atom : residue.getAtoms()) {
                    atom.setTemperatureFactor(displacements[residue.number()]);
                }
         }
    }

    public static ArrayList<DecoyAndMeta> getProteins(File inDir, CommandList commands) {
        ArrayList<DecoyAndMeta> decoys = new ArrayList<>();
        for (File file : inDir.listFiles()){
            if (file.getAbsolutePath().endsWith(".pdb")) {
                System.out.println("loading "+file);
                DecoyAndMeta decoy = new DecoyAndMeta(commands, file.getAbsolutePath());
                decoys.add(decoy);
            }
        }
        return decoys;
    }

    private static class DecoyAndMeta {
        Protein decoy;
        ArrayList<String> meta;
        public DecoyAndMeta(CommandList commands, String fileName) {
            meta = new ArrayList<>();
            MeshiLineReader reader = new MeshiLineReader((fileName));
            String line;
            try {
                while ((line = reader.readLine()) != null){
                    if (line.startsWith("ATOM"))
                        break;
                    meta.add(line);
                }
                reader.close();
            } catch (Exception ex) {throw new RuntimeException(ex);}
            decoy = Utils.getProtein(commands, fileName, ResidueExtendedAtoms.creator, Utils.defaultExceptionHandler);
        }
    }
}
