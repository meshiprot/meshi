package meshi.applications.contactMap;

import meshi.molecularElements.extendedAtoms.Pro;

import java.util.ArrayList;

public class SimpleMeanProbabilityEstimator implements ProbabilityEstimator {
    public double estimateProbability(ContactMap map, int row, int column) {
        double sum = 0;
        double sumWeights = 0;
        ContactMapCell cell = map.getCell(row, column);
        for (int i = 0; i < cell.evidences.size(); i++) {
            double weight = cell.weights.get(i);
            sum += cell.evidences.get(i) * weight;
            sumWeights += weight;
        }
        return sum / sumWeights;
    }

    public int estimateNumberOfContacts(ContactMap map){
        int numberOfContacts = 0;
        ArrayList<Integer> numbersOfContacts = map.getNumbersOfContacts();
        ArrayList<Double> weights = map.getWeights();
        double sumWeights = 0;
        for (int i = 0; i < weights.size(); i++) {
            double weight = weights.get(i);
            numberOfContacts += numbersOfContacts.get(i) * weight;
            sumWeights += weight;
        }
        return (int) Math.round(numberOfContacts/sumWeights);
    }
}
