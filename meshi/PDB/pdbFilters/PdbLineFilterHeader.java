/*
 * Copyright (c) 2010. of Chen Keasar, BGU . For free use under LGPL
 */

package meshi.PDB.pdbFilters;

import meshi.PDB.pdbLines.PdbLine;

public class PdbLineFilterHeader extends PdbLineFilter {
    public boolean acceptPdbLine(PdbLine line) {
        return line.isSEQRES() || line.isMODRES();
    }
}
