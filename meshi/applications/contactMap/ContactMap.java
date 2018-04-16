package meshi.applications.contactMap;

import meshi.dataStructures.Array;
import meshi.molecularElements.Protein;
import meshi.molecularElements.Residue;
import meshi.molecularElements.ResidueList;
import meshi.molecularElements.atoms.Atom;
import meshi.parameters.ResidueType;
import meshi.sequences.ResidueAlignment;
import meshi.sequences.ResidueAlignmentColumn;
import meshi.util.ModelAnalyzer;
import meshi.util.file.MeshiWriter;

import java.io.*;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.TreeMap;

public class ContactMap implements Serializable{
    private final TreeMap<Integer, ContactMapCell>[] map;
    private final boolean[][] visited;
    private int numberOfContacts;
    private ArrayList<Integer> numbersOfContacts = new ArrayList<>();
    private ArrayList<Double>  weights = new ArrayList<>();


    public enum Mode {CA, CB}
    public enum Type {PROTEIN, PREDICTED}

    public final double threshold;
    public final Mode mode;
    public final Type type;
    public final ArrayList<ContactMapResidue> residues;
    public ContactMap(Protein protein, double threshold, Mode mode, Type type) {
        this.type = type;
        ResidueList residueList;
        numberOfContacts = 0;
        this.mode = mode;
        this.threshold = threshold;
        residueList = protein.residues();
        residues = new ArrayList<>();
        for (Residue residue : residueList)
            residues.add(new ContactMapResidue(protein.name(),residue));
        map = new TreeMap[residues.size()];
        for (int i = 0; i < residues.size(); i++)
            map[i] = new TreeMap<>();
        for (int row = 0; row < residues.size(); row++) {
            Residue rowResidue = residueList.get(row);
            for (int column = row+1; column < residueList.size(); column++) {
                Residue columnResidue = residueList.get(column);
                if (isContact(rowResidue, columnResidue)) {
                    set(row, column, new ContactMapResidue(protein.name(), rowResidue), new ContactMapResidue(protein.name(), columnResidue), 1);
                    numberOfContacts++;
                }
                else set(row, column, new ContactMapResidue(protein.name(), rowResidue), new ContactMapResidue(protein.name(), columnResidue), 0);
            }
        }
        visited = new boolean[map.length][map.length];
        resetVisited();
    }

    public void setNumbersOfContacts(ArrayList<Integer> numbersOfContacts) {
        this.numbersOfContacts = numbersOfContacts;
    }

    public ArrayList<Double> getWeights() {
        return weights;
    }

    private void resetVisited() {
        for (int i = 0; i < visited.length; i++)
            for (int j = i + 1; j < visited.length; j++)
                visited[i][j] = false;
    }

    public boolean isContact(Residue residue1, Residue residue2) {
        Atom atom1, atom2;
        if (mode == Mode.CA) {
            atom1 = residue1.ca();
            atom2 = residue2.ca();
        } else if (mode == Mode.CB) {
            if (residue1.type != ResidueType.GLY)
                atom1 = residue1.cb();
            else atom1 = residue1.ca();
            if (residue2.type != ResidueType.GLY)
                atom2 = residue2.cb();
            else atom2 = residue2.ca();
        } else throw new RuntimeException("This is weird.");
        return atom1.distanceFrom(atom2) <= threshold;
    }

    public void set(int row, int column, ContactMapResidue rowResidue, ContactMapResidue columnResidue, double probability) {
        ContactMapCell cell = getCell(row, column, rowResidue, columnResidue);
        cell.setProbabilty(probability);
    }

    private ContactMapCell getCell(ContactMapResidue rowResidue, ContactMapResidue columnResidue) {
        int row    = getIndex(rowResidue);
        int column = getIndex(columnResidue);
        return getCell(row, column);
    }

    private int getIndex(ContactMapResidue residue) {
        int index = -1;
        for (int i = 0; i < residues.size(); i++) {
            if (residues.get(i).equals(residue))
                index = i;
        }
        if (index == -1)
            throw new RuntimeException("This is weird " + residue);
        return index;
    }

    protected ContactMapCell getCell(int row, int column) {
        if (row >= column)
            throw new RuntimeException("Row mus be smaller than collumn.");
        ContactMapCell cell = map[row].get(column);
        if (cell == null)
            throw new RuntimeException("Cell does not exists.");
        return cell;
    }

    protected ContactMapCell getCell(int row, int column, ContactMapResidue rowResidue, ContactMapResidue columnResidue) {
        if (row > column)
            throw new RuntimeException("Row mus be smaller than collumn.");
        ContactMapCell cell = map[row].get(column);
        if (cell == null) {
            cell = new ContactMapCell(row, column, rowResidue, columnResidue);
            map[row].put(column, cell);
        }
        else {
            if ((cell.rowResidue != rowResidue) | (cell.columnResidue != columnResidue))
                throw new RuntimeException("THis is weird");
        }
        return cell;
    }


    public static double MCC(ContactMap nativeMap, ContactMap modelMap, ResidueAlignment alignment, double threshold) {
       double tp = 0, tn = 0, fp = 0, fn = 0;
       boolean[] visited = new boolean[nativeMap.map.length];
       for (int i = 0; i < visited.length; i++) visited[i] = false;
        for (int iResidue = 0; iResidue < alignment.size(); iResidue++) {
            ResidueAlignmentColumn columnI = alignment.get(iResidue);
            ContactMapResidue nativeResidueI = new ContactMapResidue(null, columnI.residue0());
            ContactMapResidue modelResidueI = new ContactMapResidue(null, columnI.residue1());
            visited[nativeMap.getIndex(nativeResidueI)] = true;
            for (int jResidue = iResidue + 1; jResidue < alignment.size(); jResidue++) {
                ResidueAlignmentColumn columnJ = alignment.get(jResidue);
                ContactMapResidue nativeResidueJ = new ContactMapResidue(null, columnJ.residue0());
                ContactMapResidue modelResidueJ = new ContactMapResidue(null, columnJ.residue1());
                visited[nativeMap.getIndex(nativeResidueJ)] = true;
                double nativeContact = nativeMap.getCell(nativeResidueI, nativeResidueJ).getProbabilty();
                double modelContact = modelMap.getCell(modelResidueI, modelResidueJ).getProbabilty();
                if (modelContact > threshold)
                    modelContact = 1;
                else modelContact = 0;
                tp += nativeContact * modelContact;
                tn += (1 - nativeContact) * (1 - modelContact);
                fp += (1 - nativeContact) * modelContact;
                fn += nativeContact * (1 - modelContact);
          }
        }
        for (int i = 0; i < visited.length; i++)
            if (!visited[i]) {
                for (ContactMapCell cell : nativeMap.map[i].values()){
                    tn += 1 - cell.getProbabilty();
                    fn += cell.getProbabilty();
                }
            }
       return ModelAnalyzer.matthewsCorrelationCoefficient(tp, fn, fp, tn);
    }

    public int getNumberOfContacts() {
        return numberOfContacts;
    }

    public ArrayList<Integer> getNumbersOfContacts() {
        return numbersOfContacts;
    }

    public double estimateProbabilities(ProbabilityEstimator estimator) {
        for (int iRow = 0; iRow < map.length; iRow++) {
            for (int jColumn = iRow+1; jColumn < map.length; jColumn++) {
                ContactMapCell cell = getCell(iRow, jColumn);
                cell.setProbabilty(estimator.estimateProbability(this, iRow, jColumn));
            }
        }
        numberOfContacts = (int) Math.round(estimator.estimateNumberOfContacts(this)*1.1);
        double[] probabilities = new double[map.length*(map.length+1)];
        for (int i = 0; i < map.length; i++)
            for (int j = i + 1; j < map.length; j++)
                probabilities[map.length*i+j] = getCell(i,j).getProbabilty();
        Arrays.sort(probabilities);
        return probabilities[probabilities.length - numberOfContacts];
    }

    public void print(String fileName) {
        MeshiWriter[] writers = new MeshiWriter[8];
        try {
            for (int i = 2; i < 10; i++)
                writers[i-2] = new MeshiWriter(""+i+"."+fileName);
        } catch (Exception ex) {
            throw new RuntimeException(ex);
        }
        for (int i = 0; i < map.length; i++)
            for (int j = i+1; j < map.length; j++) {
                ContactMapCell cell = getCell(i, j);
                double th = 0.2;
                for (int ith = 2; ith < 10; ith++) {
                    double probability = cell.getProbabilty();
                    if ((probability > th) & (probability <= th + 0.100000001))
                        writers[ith - 2].println(cell.rowResidue.number + " , " + cell.columnResidue.number);
                    th = th + 0.1;
                }
            }
        for (MeshiWriter mw : writers)
            mw.close();
    }

    public void reset(double weight) {
        for (int i = 0; i< map.length; i++)
            for (int j = i + 1; j < map.length; j++) {
                ContactMapCell cell = getCell(i,j);
                cell.setProbabilty(0);
                cell.addEvidence(0, weight);
            }
        numbersOfContacts.clear();
    }
    public void addEvidences(ContactMap modelMap, ResidueAlignment alignment, double score) {
        for (int iResidue = 0; iResidue < alignment.size(); iResidue++) {
            ResidueAlignmentColumn columnI = alignment.get(iResidue);
            ContactMapResidue nativeResidueI = new ContactMapResidue(null, columnI.residue0());
            ContactMapResidue modelResidueI =  new ContactMapResidue(null, columnI.residue1());
            for (int jResidue = iResidue + 1; jResidue < alignment.size(); jResidue++) {
                ResidueAlignmentColumn columnJ = alignment.get(jResidue);
                ContactMapResidue nativeResidueJ =  new ContactMapResidue(null, columnJ.residue0());
                ContactMapResidue modelResidueJ =  new ContactMapResidue(null, columnJ.residue1());
                ContactMapCell cell = getCell(nativeResidueI, nativeResidueJ);
                double modelContact = modelMap.getCell(modelResidueI, modelResidueJ).getProbabilty();
                cell.addEvidence(modelContact, score);
                visited[getIndex(nativeResidueI)][getIndex(nativeResidueJ)] = true;
            }
        }
        for (int i = 0; i < visited.length; i++) {
            for (int j = i + 1; j < visited.length; j++)
                if (visited[i][j])
                    visited[i][j] = false;
                else getCell(i, j).addEvidence(0, score);
        }
        numbersOfContacts.add(modelMap.numberOfContacts);
        weights.add(score);
    }

    public void save(String fileName) throws IOException{
        FileOutputStream fos = new FileOutputStream(fileName);
        ObjectOutputStream oos = new ObjectOutputStream(fos);
        oos.writeObject(this);
        oos.close();
    }

    public static ContactMap load(String fileName) throws IOException, ClassNotFoundException{
        FileInputStream fis = new FileInputStream(fileName);
        ObjectInputStream ois = new ObjectInputStream(fis);
        ContactMap contactMap = (ContactMap) ois.readObject();
        ois.close();
        return contactMap;
    }
}
