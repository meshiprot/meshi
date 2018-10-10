package meshi.PDB.pdbLines;

import meshi.parameters.ResidueType;

/**
 * Created by user on 01/11/2017.
 */
public class PdbSEQRES {
    public int residueInx;
    public int lineID;
    public String chain;
    public String residueExperimentName;
    public ResidueType residueGenericType;
    public String residueGenericName;

    public PdbSEQRES(){
        throw new RuntimeException("PdbSEQRESLine: ERROR - No sequence resdidue entered - unable to use the empty constructor.");
    }

    public PdbSEQRES(int residueNum, int lineID, String chain, String residueName){
        this.residueExperimentName = residueName;
        this.lineID = lineID;
        this.residueInx = residueNum;
        this.residueGenericName = this.residueExperimentName;
        this.residueGenericType = ResidueType.type(this.residueGenericName);
        this.chain = chain;
    }
    public boolean check(){
        return (this.residueInx>=0);
    }

}
