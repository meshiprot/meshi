/*
 * Copyright (c) 2010. of Chen Keasar, BGU . For free use under LGPL
 */

package meshi.PDB.pdbFilters;

import meshi.PDB.pdbLines.PdbLine;

public class PdbLineATOMorHETATOM extends PdbLineFilter {
    public PdbLineATOMorHETATOM() {
    }

    public boolean acceptPdbLine(PdbLine line) {
        return line.isAnAtomOrHeteroAtom();
    }
}
