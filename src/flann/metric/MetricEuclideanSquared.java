package flann.metric;

import flann.exception.ExceptionFLANN;

/**
 * Euclidean squared distance.
 */
public class MetricEuclideanSquared implements Metric {
	@Override
	public double distance(double[] a, double[] b) {
		double result = 0.0;
		for (int i = 0; i < a.length; i++) {
			result += distance(a[i], b[i]);
		}
		return result;
	}

	@Override
	public double distance(double a, double b) {
		double diff = a - b;
		return diff * diff;
	}

	@Override
	public int distance(int[] a, int[] b) {
		throw new ExceptionFLANN("Unsupported types");
	}

	@Override
	public int distance(int a, int b) {
		throw new ExceptionFLANN("Unsupported types");
	}
}