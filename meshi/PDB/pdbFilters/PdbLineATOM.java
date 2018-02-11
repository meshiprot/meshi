/*
 * Copyright (c) 2010. of Chen Keasar, BGU . For free use under LGPL
 */

package meshi.PDB.pdbFilters;

import meshi.PDB.pdbLines.PdbATOMLine;
import meshi.PDB.pdbLines.PdbLine;
import meshi.util.filters.Filter;

import java.util.ArrayList;
import java.util.List;

public class PdbLineATOM extends PdbLineFilter {
    public final List<String> chainNames;
    public static final Filter filter = new PdbLineATOM();

    public PdbLineATOM() {
        this(new ArrayList<String>());
    }

    public PdbLineATOM(List<String> chainNames) {
        this.chainNames = chainNames;
    }

    public boolean acceptPdbLine(PdbLine line) {
        if (chainNames.size() == 0) return line.isAnAtom();
        for (String name : chainNames) {
            if (line.isAnAtom() && ((PdbATOMLine)line).chain().equals(name)) return true;
        }
        return false;
    }
}
