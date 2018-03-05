/*
 * Copyright (c) 2010. of Chen Keasar, BGU . For free use under LGPL
 */

package meshi.PDB.pdbFilters;

import meshi.PDB.pdbLines.PdbATOMLine;
import meshi.PDB.pdbLines.PdbLine;
import meshi.util.filters.Filter;

import java.util.ArrayList;
import java.util.List;

public class PdbLineMultipleModelsFilter extends PdbLineFilter {
    public final List<String> chainNames;
    public static final Filter filter = new PdbLineMultipleModelsFilter();

    public PdbLineMultipleModelsFilter() {
        this(new ArrayList<String>());
    }

    public PdbLineMultipleModelsFilter(List<String> chainNames) {
        this.chainNames = chainNames;
    }

    public boolean acceptPdbLine(PdbLine line) {
        if (chainNames.size() == 0) return line.isAnAtom() || line.isEnd();
        for (String name : chainNames) {
            if (line.isAnAtom() && ((PdbATOMLine)line).chain().equals(name)) return true;
            if (line.isEnd())
                return true;
        }
        return false;
    }
}
