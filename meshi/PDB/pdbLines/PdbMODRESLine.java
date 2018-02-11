package meshi.PDB.pdbLines;

/**
 * Created by user on 26/10/2017.
 */
public class PdbMODRESLine extends PdbLine {

    public String proteinName;
    public String residueExperimentName;
    public String residueGenericName;
    public String chain;
    public int residueInx;

    public PdbMODRESLine(){
        super(PdbLineType.MODRES);
        this.proteinName = "XX";
        this.residueExperimentName="X";
        this.residueGenericName="X";
        this.chain="x";
        this.residueInx =-1;
        throw new RuntimeException("PdbMODRESLine: ERROR - No line entered.");
    }

    public PdbMODRESLine(String line){
        super(PdbLineType.MODRES);
        String[] args = line.split("\\s+");
        if (args.length < 7) throw new RuntimeException("PdbMODRESLine: ERROR - MODRES line configuration error, " + args.length +"/7 items in line. "+line);
        this.proteinName = args[1];
        this.residueExperimentName = args[2];
        this.chain = args[3];
        try { this.residueInx = Integer.parseInt(args[4].trim());} catch (Exception e) {throw new RuntimeException("PdbMODRESLine: ERROR - MODRES line configuration error, line number is not a number. " + line);}
        this.residueGenericName = args[5];
        this.line = line;
    }


}
