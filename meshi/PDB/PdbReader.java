/*
 * Copyright (c) 2010. of Chen Keasar, BGU . For free use under LGPL
 */

package meshi.PDB;

import meshi.PDB.pdbLines.PdbATOMLine;
import meshi.PDB.pdbLines.PdbLine;
import meshi.PDB.pdbLines.PdbLineGenerator;
import meshi.util.MeshiException;
import meshi.util.file.MeshiLineReader;
import meshi.util.filters.Filter;

import java.io.File;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.StringTokenizer;

public class PdbReader extends MeshiLineReader {

    PdbHeader pdbHeader;

    public PdbReader(InputStreamReader inputStreamReader, String path) {
        super(inputStreamReader, path);
    }


    public PdbReader(String name) {
        super(getPdbFile(name));
    }
    public PdbReader(String name,long loc) throws IOException {
        super(getPdbFile(name));
        this.skip(loc);
    }
    public PdbReader(String dir, String name) {
        super(getPdbFile(dir, name));
    }
    public PdbReader(String name, PdbHeader pdbHeader) {
        super(getPdbFile("", name));
        this.pdbHeader = pdbHeader;
    }

    private static File getPdbFile(String name) {
        return getPdbFile("", name);
    }


    private static File getPdbFile(String dir, String name) {
        File temp;
        String alternativeName;
        if (name.length() <= 3)
            throw new MeshiException("PDBReader.getPdbFile(String name) error.\n" +
                    "name too short: " + name + "\n");
        String treeDir = name.substring(1, 3);
        String treeDir2 = name.substring(2, 4);
        if (!dir.equals("")) dir = dir + "/";
        String alternativeNamesString =
                dir + name + " \n" +
                        dir + name + ".gz \n" +
                        dir + name + ".pdb \n" +
                        dir + name + ".pdb.gz \n" +
                        dir + name + ".ent \n" +
                        dir + name + ".ent.gz \n" +
                        dir + "pdb" + name + ".ent \n" +
                        dir + "pdb" + name + ".ent.gz \n" +
                        dir + treeDir + "/" + name + " \n" +
                        dir + treeDir + "/" + name + ".pdb \n" +
                        dir + treeDir + "/" + name + ".pdb.gz \n" +
                        dir + treeDir + "/" + "pdb" + name + ".ent \n" +
                        dir + treeDir + "/" + "pdb" + name + ".ent.gz \n" +
                        dir + treeDir2 + "/" + name + " \n" +
                        dir + treeDir2 + "/" + name + ".pdb \n" +
                        dir + treeDir2 + "/" + name + ".pdb.gz \n" +
                        dir + treeDir2 + "/" + name + ".ent \n" +
                        dir + treeDir2 + "/" + name + ".ent.gz \n";
        StringTokenizer alternativeNames = new StringTokenizer(alternativeNamesString);
        while (alternativeNames.hasMoreTokens()) {
            alternativeName = alternativeNames.nextToken();
            temp = new File(alternativeName);
            if (temp.exists()) return temp;
        }
        throw new MeshiException("PDBReader.getPdbFile(" + dir + "," + name + ") error.\n" +
                "file not found in any of these alternatives:\n" +
                alternativeNamesString + "\n");
    }

    public PdbLine readPdbLine() {
        String line;
        try {
            line = readLine();
        }
        catch (Exception e) {
            throw new MeshiException("PDBReader Error1 while trying to read file:" +
                    "File name: " + path() + "\n" +
                    e.getMessage() + "\n");
        }
        if (line == null) return null;

        // PdbLine pdbLine = PdbATOMLine.getPdbLine(line);
        PdbLine pdbLine = PdbLineGenerator.getPdbLine(line);
        this.updateLine(pdbLine);
        return pdbLine;
    }

    private void updateLine(PdbLine pdbLine) {
        if (pdbHeader != null && pdbLine.isAHeteroAtom()){
            ((PdbATOMLine) pdbLine).modifyByHeader(pdbHeader);
        }
    }

    public PdbLine readPdbLine(Filter filter) {
        String line;
        PdbLine pdbLine;

        line = "";
        while (line != null) {
            try {
                line = readLine();
            }
            catch (Exception e) {
                throw new MeshiException("PDBReader Error2 while trying to read file:" +
                        "File name: " + path() + "\n" +
                        e.getMessage() + "\n");
            }
            if (line == null) return null;
            pdbLine = PdbLineGenerator.getPdbLine(line);
            this.updateLine(pdbLine);
            if (filter.accept(pdbLine)) return pdbLine;
        }
        return null;
    }
    public String toStringSequenceHeader(){
        if (pdbHeader != null) return pdbHeader.toStringSequence();
        return "No Sequence Header";
    }
}
