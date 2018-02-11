/*
 * Copyright (c) 2010. of Chen Keasar, BGU . For free use under LGPL
 */

package meshi.PDB.pdbFilters;

import meshi.PDB.pdbLines.PdbLine;
import meshi.util.filters.Filter;

public abstract class PdbLineFilter implements Filter {
    public boolean accept(Object obj) {
        return acceptPdbLine((PdbLine) obj);
    }

    public abstract boolean acceptPdbLine(PdbLine line);
}
