package meshi.applications.contactMap;

import meshi.util.ValueWeight;

import java.util.ArrayList;

public class ExponentialMeanProbabilityEstimator implements ProbabilityEstimator {
    private double factor;
    public ExponentialMeanProbabilityEstimator(double factor){
        this.factor = factor;
    }
    public double estimateProbability(ContactMap map, int row, int column) {
        double sum = 0;
        double sumWeights = 0;
        ContactMapCell cell = map.getCell(row, column);
        ValueWeight[] valueWeights = new ValueWeight[cell.evidences.size()];
        for (int i = 0; i < cell.evidences.size(); i++) {
            valueWeights[i] = new ValueWeight(cell.evidences.get(i), Math.exp(factor*cell.weights.get(i)));
        }
        return ValueWeight.weightedMean(valueWeights);
    }

    public int estimateNumberOfContacts(ContactMap map){
        int numberOfContacts = 0;
        ArrayList<Integer> numbersOfContacts = map.getNumbersOfContacts();
        ArrayList<Double> weights = map.getWeights();
        double sumWeights = 0;
        for (int i = 0; i < weights.size(); i++) {
            double weight = weights.get(i);
            numberOfContacts += numbersOfContacts.get(i) * Math.exp(factor*weight);
            sumWeights += Math.exp(factor*weight);
        }
        return (int) Math.round(numberOfContacts/sumWeights);
    }
}
