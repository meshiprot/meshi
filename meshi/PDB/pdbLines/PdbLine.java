package meshi.PDB.pdbLines;

/**
 * Created by user on 26/10/2017.
 */
public class PdbLine {

    protected String line;
    private final PdbLineType pdbLineType;

    protected PdbLine(PdbLineType plt){
        this.pdbLineType = plt;
    }

    protected PdbLine(){
        this.pdbLineType = PdbLineType.BASIC;
        throw new RuntimeException("PdbLine: PDB line must have a pdbLineType");
    }
    public PdbLine(String line){
        this.pdbLineType = PdbLineType.BASIC;
        this.line = line;
    }

    public PdbLineType getPdbLineType(){
        return pdbLineType;
    }

    public boolean isSEQRES() {
        return line.startsWith("SEQRES");
    }

    public boolean isMODRES() {
        return line.startsWith("MODRES");
    }

    public boolean isAnAtom() {
        return ((line.length() >= 54) && line.startsWith("ATOM"));
    }

    public boolean isAHeteroAtom() {
        return ((line.length() >= 54) && line.startsWith("HETATM"));
    }

    public boolean isAnAtomOrHeteroAtom() {
        return (isAnAtom() || isAHeteroAtom());
    }

    public boolean isAComment() {
        return (!isAnAtomOrHeteroAtom());
    }

    /*
     *Check if this is a MODEL line.
     */
    public boolean isAModel() {
        return line.startsWith("MODEL");
    }

    public int length() {
        return line.length();
    }

    public boolean isTer() {
        return line.startsWith("TER");
    }

    public boolean isEnd() {
        return line.startsWith("END");
    }

    public String toString(){
        return line;
    }



}
