package meshi.util;

import meshi.molecularElements.Protein;
import meshi.util.file.MeshiWriter;

public class PdbWriter {
    protected MeshiWriter writer;
    protected Protein protein;

    public PdbWriter(MeshiWriter writer, Protein protein) {
        this.writer = writer;
        this.protein = protein;
    }

    public void print() {
        protein.atoms().print(writer);
    }
}
