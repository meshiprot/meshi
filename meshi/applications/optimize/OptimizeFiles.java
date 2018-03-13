package meshi.applications.optimize;


import meshi.molecularElements.extendedAtoms.ResidueExtendedAtoms;
import meshi.util.Utils;

import java.io.File;


public class OptimizeFiles {
    public final File modelIn;
    public final File dssp;
    public final File nativeStructure;
    public final File modelOut;
    public final File infoXML;
    public final File infoTable;

    public OptimizeFiles(String inFileName, String dsspFileName, String nativeFileName, String outFileName) {
        modelIn = new File(inFileName);
        dssp = new File(dsspFileName);
        nativeStructure = new File(nativeFileName);
        modelOut = new File(outFileName);
        infoXML = new File(outFileName+".info.xml");
        infoTable = new File(outFileName+".info.csv");


    }
}
