package flann.result_set;



public interface ResultSet {
	public boolean full();
	public void addPoint (double distance, int index);
	public double worstDistance();
}