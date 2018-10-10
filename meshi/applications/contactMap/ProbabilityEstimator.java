package meshi.applications.contactMap;

public interface ProbabilityEstimator {
    public double estimateProbability(ContactMap map, int row, int column);
    public int estimateNumberOfContacts(ContactMap map);
}
