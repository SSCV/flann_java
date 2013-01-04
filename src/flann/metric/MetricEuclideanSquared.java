package flann.metric;



/**
 * Euclidean squared distance.
 */
public class MetricEuclideanSquared implements Metric {
	public double distance (double[] a, double[] b) {
		double result = 0.0, diff;
		for (int i = 0; i < a.length; i++) {
			diff = a[i] - b[i];
			result += diff * diff;
		}
		return result;
	}

	public double distance (double a, double b) {
		double diff = a - b;
		return diff * diff;
	}
}