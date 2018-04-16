package meshi.applications.contactMap;

import meshi.molecularElements.Residue;

import java.io.Serializable;
import java.util.ArrayList;

public class ContactMapCell implements Comparable<ContactMapCell>, Serializable{
    public final int row, column;
    public final ContactMapResidue rowResidue, columnResidue;



    private double probabilty;
    protected ArrayList<Double> evidences;
    protected ArrayList<Double> weights;
    public ContactMapCell(int row, int column, ContactMapResidue rowResidue, ContactMapResidue columnResidue) {
     this.column = column;
     this.row = row;
     this.rowResidue = rowResidue;
     this.columnResidue = columnResidue;
     probabilty = 0;
     evidences = new ArrayList<>();
     weights   = new ArrayList<>();
    }

    public double getProbabilty() {
        return probabilty;
    }

    public void setProbabilty(double probabilty) {
        this.probabilty = probabilty;
    }

    public void addEvidence(double evidence, double weight) {
        evidences.add(evidence);
        weights.add(weight);
    }

    public int compareTo(ContactMapCell cell) {
        if (cell == null)
            throw new RuntimeException("This is weird");
        if (probabilty  > cell.probabilty) return 1;
        if (probabilty  < cell.probabilty) return -1;
        return 0;
    }
}
