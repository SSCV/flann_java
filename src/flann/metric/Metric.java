package flann.metric;



public interface Metric {
	public double distance (double[] a, double[] b);
	public double distance (double a, double b);
}