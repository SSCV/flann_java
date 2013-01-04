package flann.index;

import flann.metric.*;



/**
 * Base class for all index structures.
 */
public class IndexBase {
	protected Metric metric;
	protected double[][] data;

	public int numberOfObjects;
	public int numberOfDimensions;

	public IndexBase (Metric metric, double[][] data) {
		this.metric = metric;
		setDataset (data);
	}

	protected void setDataset (double[][] data) {
		if (data.length == 0)
			return;

		numberOfObjects = data.length;
		numberOfDimensions = data[0].length;

		// Copy data.
		this.data = new double[numberOfObjects][numberOfDimensions];
		for (int i = 0; i < numberOfObjects; i++) {
			for (int j = 0; j < numberOfDimensions; j++) {
				this.data[i][j] = data[i][j];
			}
		}
	}
}
