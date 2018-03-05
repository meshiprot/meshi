package meshi.util;

import meshi.PDB.pdbFilters.PdbLineMultipleModelsFilter;
import meshi.PDB.pdbIterator.PdbProteinIterator;
import meshi.molecularElements.Protein;
import meshi.molecularElements.ResidueCreator;
import meshi.molecularElements.atoms.AtomList;
import meshi.util.filters.Filter;

/**
 * Created by user on 01/03/2018.
 */
public class ProteinGenerator {

    public PdbProteinIterator pit;

    public ProteinGenerator(String inPDBFile, long loc, int numOfModelsInGroup){
        try {
            pit = new PdbProteinIterator(inPDBFile, loc, numOfModelsInGroup,new PdbLineMultipleModelsFilter());
        } catch (Exception ex) {
            System.err.println(ex.getStackTrace());
        }
    }
    public ProteinGenerator(String inPDBFile, long loc, int numOfModelsInGroup, Filter pdbFilter){
        try {
            pit = new PdbProteinIterator(inPDBFile, loc, numOfModelsInGroup, pdbFilter);
        } catch (Exception ex) {
            System.err.println(ex.getStackTrace());
        }
    }
    public Protein getProtein(CommandList commands, ResidueCreator creator, ExceptionHandler exceptionHandler) {
        try {
            Protein out = new Protein(new AtomList(pit), creator, commands);
            Utils.addHydrogens(out, commands);
            pit.nextProtein();

            return out;
        } catch (Exception ex) {
            exceptionHandler.handle(ex);
        }
        return null;
    }
    public String getPdbAsString(){
        if (pit != null) return pit.getLastProtein();
        else return null;
    }

}
