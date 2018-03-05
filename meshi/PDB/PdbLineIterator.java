package meshi.PDB;

import meshi.PDB.pdbLines.PdbLine;
import meshi.util.file.MeshiLineReader;
import meshi.util.filters.Filter;

import java.util.Iterator;

/**
 * Created by chen on 15/07/2017.
 */
public class PdbLineIterator implements Iterator<PdbLine> {
    private Filter filter;
    private PdbReader reader;
    private PdbLine nextLine;

    protected PdbLineIterator(){

    }

    public PdbLineIterator(PdbReader reader){
        this.reader = reader;
        this.filter = null;
        nextLine = reader.readPdbLine();

    }

    public PdbLineIterator(PdbReader reader, Filter pdbFilter){
        this.reader = reader;
        this.filter = pdbFilter;
        nextLine = reader.readPdbLine(pdbFilter);
        while (nextLine!=null && !this.filter.accept(nextLine)) nextLine = reader.readPdbLine();
    }

    protected PdbLine getLine(){
        return this.nextLine;
    }

    public PdbLine next() {
        PdbLine out = nextLine;
        nextLine = reader.readPdbLine();
        while (this.filter != null && nextLine!=null && !this.filter.accept(nextLine)) {
            nextLine = reader.readPdbLine();
        }
        return out;
    }
    public boolean hasNext() {
        if (nextLine == null) return false;
        return (!nextLine.isEnd());
    }

    public MeshiLineReader sourceFile() {
        return reader;
    }

    public void remove() {throw new RuntimeException("Unsupported operation");}

}
