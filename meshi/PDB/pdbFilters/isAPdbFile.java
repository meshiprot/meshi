/*
 * Copyright (c) 2010. of Chen Keasar, BGU . For free use under LGPL
 */

package meshi.PDB.pdbFilters;

import meshi.util.filters.Filter;

import java.io.File;

public class isAPdbFile implements Filter {
    public boolean accept(Object obj) {
        if (!(obj instanceof File))
            throw new RuntimeException("weird input to isAPdbFile.accept\n" +
                    obj);
        File file = (File) obj;
        String name = file.getName();
        if (name.endsWith(".pdb")) return true;
        if (name.endsWith(".pdb.gz")) return true;
        if (name.endsWith(".ent.gz")) return true;
        if (name.endsWith(".ent")) return true;
        return false;
    }
}
