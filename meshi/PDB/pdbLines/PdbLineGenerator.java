package meshi.PDB.pdbLines;

/**
 * Created by user on 02/11/2017.
 */
public class PdbLineGenerator {

    public static PdbLine getPdbLine(String line){
        PdbLineType pdbLineType = PdbLineGenerator.getPdbLineType(line);
        PdbLine pdbLine = new PdbLine(line);;
        switch (pdbLineType){
            case ATOM:
                pdbLine = new PdbATOMLine(line);
                break;
            case MODRES:
                pdbLine = new PdbMODRESLine(line);
                break;
            case SEQRES:
                pdbLine = new PdbSEQRESLine(line);
                break;
            case HETATOM:
                pdbLine = new PdbATOMLine(line);
                break;
            case BASIC:
                pdbLine = new PdbLine(line);
                break;
        }
        return pdbLine;

    }

    private static PdbLineType getPdbLineType(String line){
        if (PdbLineGenerator.isAnAtom(line)) return PdbLineType.ATOM;
        if (PdbLineGenerator.isMODRES(line)) return PdbLineType.MODRES;
        if (PdbLineGenerator.isSEQRES(line)) return PdbLineType.SEQRES;
        if (PdbLineGenerator.isAHeteroAtom(line)) return PdbLineType.HETATOM;
        return PdbLineType.BASIC;
    }

    /*******  Static methods - PDB line generator *******/


    public static boolean isSEQRES(String line) {
        return line.startsWith("SEQRES");
    }

    public static boolean isMODRES(String line) {
        return line.startsWith("MODRES");
    }

    public static boolean isAnAtom(String line) {
        return ((line.length() >= 54) && line.startsWith("ATOM"));
    }

    public static boolean isAHeteroAtom(String line) {
        return ((line.length() >= 54) && line.startsWith("HETATM"));
    }

    public static boolean isAnAtomOrHeteroAtom(String line) {
        return (PdbLineGenerator.isAnAtom(line) || PdbLineGenerator.isAHeteroAtom(line));
    }

    public static boolean isAComment(String line) {
        return (!PdbLineGenerator.isAnAtomOrHeteroAtom(line));
    }

}
