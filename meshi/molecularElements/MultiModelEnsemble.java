package meshi.molecularElements;

import meshi.molecularElements.atoms.AtomList;
import meshi.molecularElements.extendedAtoms.ResidueExtendedAtomsCreator;
import meshi.util.file.MeshiLineReader;
import meshi.util.file.MeshiWriter;

import java.io.IOException;
import java.util.ArrayList;

/**
 * Created by chen on 12/09/2016.
 */
public class MultiModelEnsemble extends ArrayList<Protein>{
    private enum Status {IN_MODEL, OUT}
    public MultiModelEnsemble(String fileName) throws IOException{
        super();
        if (isMultiModelEnsemble(fileName)) {
            MeshiLineReader reader = new MeshiLineReader(fileName);
            String line;
            Status status = Status.OUT;
            MeshiWriter tmpWriter = null;
            int modelNumber = 1;
            while ((line = reader.readLine()) != null) {
                if (status == Status.OUT) {
                    if (line.startsWith("MODEL")) {
                        status = Status.IN_MODEL;
                        tmpWriter = new MeshiWriter("temp.pdb");
                    }
                } else {
                    if (line.startsWith("END")) {
                        status = Status.OUT;
                        if (tmpWriter != null) tmpWriter.close();
                        Protein tmpProtein = new Protein(new AtomList("temp.pdb"), new ResidueExtendedAtomsCreator());
                        tmpProtein.setName(fileName+"_MODEL"+modelNumber);
                        modelNumber++;
                        add(tmpProtein);
                    } else tmpWriter.println(line);
                }
            }
            reader.close();
        }
        else add(new Protein(new AtomList(fileName), new ResidueExtendedAtomsCreator()));
    }

    public static boolean isMultiModelEnsemble(String fileName) throws IOException{
        boolean out = false;
        MeshiLineReader reader = new MeshiLineReader(fileName);
        String line;
        Status status = Status.OUT;
        while ((line = reader.readLine()) != null) {
            if (status == Status.OUT) {
                if (line.startsWith("MODEL")) {
                    status = Status.IN_MODEL;
                }
            }
            else {
                if (line.startsWith("END")){
                    status = Status.OUT;
                    out = true;
                }
            }
        }
        reader.close();
        return out;
    }
}
