package meshi.PDB.pdbLines;

/**
 * Created by user on 26/10/2017.
 */
public class PdbSEQRESLine extends PdbLine {

    public int lineID;
    public String chain;
    public PdbSEQRES[] residues;

    private final static int residuesNumInLine = 13;


    public PdbSEQRESLine(){
        super(PdbLineType.SEQRES);
        this.residues = new PdbSEQRES[0];
        this.chain="x";
        throw new RuntimeException("PdbSEQRESLine: ERROR - No line entered.");
    }

    public PdbSEQRESLine(String line){
        super(PdbLineType.SEQRES);
        this.line = line;

        if (line.charAt(11) == ' ') line = line.substring(0,11) + "A" + line.substring(13,line.length());
        String[] args = line.split("\\s+");
        if (args.length < 5) throw new RuntimeException("PdbSEQRESLine: ERROR - SEQRES line configuration error, " + args.length +"/5 items in line. "+line);
        try { this.lineID = Integer.parseInt(args[1].trim());} catch (Exception e) {throw new RuntimeException("PdbSEQRESLine: ERROR - SEQRES line configuration error, lineID number is not a number. " + line);}
        this.chain = args[2];
        this.residues = new PdbSEQRES[args.length - 4];
        for (int i=4; i<args.length; i++){
            int resNum = (i-4) + PdbSEQRESLine.residuesNumInLine * (this.lineID - 1) +1;
            this.residues[i-4] = new PdbSEQRES(resNum,lineID,chain,args[i]);
        }
    }


}
