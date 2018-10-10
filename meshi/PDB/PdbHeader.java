package meshi.PDB;

import meshi.PDB.pdbLines.*;

import java.util.ArrayList;

import static meshi.PDB.pdbLines.PdbLineType.SEQRES;

/**
 * Created by user on 26/10/2017.
 * Currently holds only MODRES lines. TODO - REMARK and SEQRES
 */
public class PdbHeader {
    private String pdbFilePath;
    private ArrayList<PdbMODRESLine> residueNameModifiers;
    private ArrayList<PdbSEQRES> residueSequence;

    public PdbHeader(){
        this.residueNameModifiers = new ArrayList<PdbMODRESLine>();
        this.residueSequence = new ArrayList<PdbSEQRES>();
    }
    public PdbHeader(String pdbFilePath){
        init(pdbFilePath);
    }

    public void init(String pdbFilePath){
        this.pdbFilePath = pdbFilePath;
        this.residueNameModifiers = new ArrayList<PdbMODRESLine>();
        this.residueSequence = new ArrayList<PdbSEQRES>();
        PdbLineIterator pdbLineIterator = new PdbLineIterator(new PdbReader(this.pdbFilePath));
        while (pdbLineIterator.hasNext()){
            PdbLine pdbLine = pdbLineIterator.next();
            if (pdbLine.isMODRES()) this.addModifier(((PdbMODRESLine) pdbLine));
            if (pdbLine.isSEQRES()){
                this.addSequenceSnipet((PdbSEQRESLine) pdbLine);
            }
        }
    }

    public void addModifier(PdbMODRESLine modres){
        if (modres == null) throw new RuntimeException("PdbHeader:addModifier: Error - the added modres line is null.");
        if (modres.getPdbLineType() != PdbLineType.MODRES) throw new RuntimeException("PdbHeader:addModifier: Error - the added line is not a modifier (MODRES).");
        residueNameModifiers.add(modres);
    }

    public void addSequenceSnipet(PdbSEQRESLine seqres){
        if (seqres == null) throw new RuntimeException("PdbHeader:addSequenceSnipet: Error - the added SEQRES line is null.");
        if (seqres.getPdbLineType() != SEQRES) throw new RuntimeException("PdbHeader:addSequenceSnipet: Error - the added line is not a sequence identifier (SEQRES).");
        for (PdbSEQRES sr: seqres.residues) {
            if (sr.check()){
                int i = this.residueSequence.size();
                this.residueSequence.add(i,sr);
            }
        }
    }

    public void updateSequenceWithModifiedResidues(){

        for (PdbMODRESLine modifier:residueNameModifiers){
            PdbSEQRES seqres = residueSequence.get(modifier.residueInx -1 );
            if (seqres.residueInx != modifier.residueInx) throw new RuntimeException("PdbHeader:updateSequenceWithModifiedResidues: Error - SEQRES index and MODRES index are different. "+seqres.lineID+" "+modifier.residueInx);
            seqres.residueGenericName = modifier.residueGenericName;
        }
    }

    public String getTrueResidueName(PdbATOMLine pdbLine) {
        if (!pdbLine.isAnAtomOrHeteroAtom()) throw new RuntimeException("PdbHeader:getTrueResidueName: Error - pdbLine is not an atom line.");
        for (PdbMODRESLine modifier:residueNameModifiers){
            if (modifier.residueInx == pdbLine.residueNumber()
                    && modifier.chain.equals(pdbLine.chain())){
                return modifier.residueGenericName;
            }
        }
        return null;
    }

    public ArrayList<String> getChainsIDs(){
        ArrayList<String> chains = new ArrayList<String>();
        PdbSEQRES seqres = this.residueSequence.get(1);

        String currChain = seqres.chain;
        String lastChain = currChain;
        chains.add(currChain);
        for (int i=1; i < this.residueSequence.size(); i++ ){
            seqres = this.residueSequence.get(i);
            currChain = seqres.chain;

            if (!currChain.equals(lastChain)){
                chains.add(currChain);
                lastChain=currChain;
            }
        }
        return chains;
    }

    public String toStringSequence(){
        String str = "";
        PdbSEQRES seqres;
        for (int i=1; i < this.residueSequence.size(); i++ ){
            seqres = this.residueSequence.get(i);
            String chain = seqres.chain;
            str += ">Chain"+chain+"\n";
            str+=seqres.residueGenericType.nameOneLetter();
            i++;
            while (i < this.residueSequence.size() && chain.equals(seqres.chain)){
                seqres = this.residueSequence.get(i);
                str+=seqres.residueGenericType.nameOneLetter();
                i++;
            }
            str+="\n";
        }
        return str;
    }
    public String toStringSequence(String c){
        String str = "";
        PdbSEQRES seqres;
        int i=0;
        seqres = this.residueSequence.get(i);
        while (i < this.residueSequence.size() && !this.residueSequence.get(i).chain.equals(c)){
            i++;
        }
        seqres = this.residueSequence.get(i);
        str += ">Chain"+seqres.chain+"\n";
        while (i < this.residueSequence.size() && seqres.chain.equals(c)){
            //System.out.println("inx="+i+ " chain="+seqres.chain+" resID="+seqres.residueGenericType);
            str += seqres.residueGenericType.nameOneLetter();
            i++;
            if (i < this.residueSequence.size()) seqres = this.residueSequence.get(i);





        }
        str+="\n";

        return str;
    }
}
