package meshi.applications.optimize;

import meshi.molecularElements.Protein;
import meshi.sequences.ResidueAlignmentMethod;
import meshi.util.Rms;
import meshi.util.filters.Filter;

import java.io.File;

public class MyFileFilter implements Filter {
    private String prefix = null;
    Protein thisModel;

    public void reset(Protein thisModel) {
        this.thisModel = thisModel;
        int index = thisModel.name().indexOf('.');
        if (index != -1) prefix = thisModel.name().substring(0, index);
        else prefix = thisModel.name();

    }

    public boolean accept(Object obj) {
        if (prefix == null) throw new RuntimeException("prefix is " + prefix);
        File file = (File) obj;
        String path = file.getAbsolutePath();
        if (!path.endsWith("pdb")) return false;
        if (path.indexOf("out") == -1) return false;
        if (file.getName().startsWith(prefix)) {
            try {
                double rms = Rms.rms(thisModel, Protein.getCAproteinFromApdbFile(file), ResidueAlignmentMethod.IDENTITY);
                if (rms < 1) return false;
            } catch (Exception ex ) {
                throw new MeshiFailureException(ex.getMessage());}
            return true;
        }
        return false;
    }
}
