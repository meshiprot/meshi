package meshi.PDB.pdbIterator;

import meshi.PDB.PdbLineIterator;
import meshi.PDB.PdbReader;
import meshi.PDB.pdbLines.PdbLine;
import meshi.util.filters.Filter;

import java.io.IOException;

/**
 * Created by user on 01/03/2018.
 */
public class PdbProteinIterator extends PdbLineIterator {
    private int numOfModels = 0;
    private int inxModel;
    private String lastProtein;
    private String currProtein;

    public PdbProteinIterator(String path, long fLoc, int numOfModels) throws IOException {
        super(new PdbReader(path,fLoc));
        this.numOfModels = numOfModels;
        this.inxModel = 0;
        this.currProtein = "";
    }
    public PdbProteinIterator(String path, long fLoc, int numOfModels, Filter pdbFilter) throws IOException {
        super(new PdbReader(path,fLoc),pdbFilter);
        this.numOfModels = numOfModels;
        this.inxModel = 0;
        this.currProtein = "";
    }
    public boolean nextProtein(){
        lastProtein = currProtein;
        currProtein = "";

        if (inxModel >= numOfModels) return false;
        while (this.hasNext()){
            super.next();
        }
        if (this.getLine() == null) return false;
        super.next();
        if (this.getLine() == null) return false;
        inxModel++;
        return true;
    }
    @Override
    public PdbLine next(){
        PdbLine line = super.next();
        currProtein += line +"\n";
        return line;
    }

    public String getLastProtein(){
        return this.lastProtein;
    }

}
