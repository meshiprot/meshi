package meshi.util;

import java.util.Arrays;
import java.util.List;

public class ValueWeight implements Comparable<ValueWeight>{
    double value, weight;

    public ValueWeight(double value, double weight) {
        this.value = value;
        this.weight = weight;
    }

    @Override
    public int compareTo(ValueWeight o) {
        if (o == null) throw new RuntimeException("This is weird.");
        if (weight > o.weight) return 1;
        if (weight < o.weight) return -1;
        return 0;
    }

    public static double weightedMean(ValueWeight[] values) {
        double sum =0, sumWeights = 0;
        for (int i = 0; i < values.length; i++){
            sum += values[i].value*values[i].weight;
            sumWeights += values[i].weight;
        }
        return sum/sumWeights;
    }
    public static double weightedMedian(ValueWeight[] values) {
        Arrays.sort(values);
        double[] sumWeights = new double[values.length];
        sumWeights[0] = values[0].weight;
        for (int i = 1; i < sumWeights.length; i++)
            sumWeights[i] = sumWeights[i-1]+values[i].weight;
        double medianValue = sumWeights[sumWeights.length-1]/2;
        int medianIndex = linearSearch(sumWeights, medianValue);
        return values[medianIndex].value;
    }

    private static int linearSearch(double[] array, double key) {
        for (int i = 0; i < array.length - 1; i++) {
            if (array[i] == key) return i;
            if ((array[i]< key) & (array[i+1] > key))
                return i;
        }
        return -1;
    }

    private static int binarySearch0(double[] a, int fromIndex, int toIndex,
                                     double key) {
        int low = fromIndex;
        int high = toIndex - 1;

        while (low <= high) {
            int mid = (low + high) >>> 1;
            double midVal = a[mid];

            if (midVal < key)
                low = mid + 1;  // Neither val is NaN, thisVal is smaller
            else if (midVal > key)
                high = mid - 1; // Neither val is NaN, thisVal is larger
            else {
                long midBits = Double.doubleToLongBits(midVal);
                long keyBits = Double.doubleToLongBits(key);
                if (midBits == keyBits)     // Values are equal
                    return mid;             // Key found
                else if (midBits < keyBits) // (-0.0, 0.0) or (!NaN, NaN)
                    low = mid + 1;
                else                        // (0.0, -0.0) or (NaN, !NaN)
                    high = mid - 1;
            }
        }
        return -(low + 1);  // key not found.
    }

}
