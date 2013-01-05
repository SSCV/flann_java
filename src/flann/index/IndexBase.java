package flann.index;

import java.util.ArrayList;
import java.util.Collections;

import flann.metric.*;



/**
 * Base class for all index structures.
 */
public class IndexBase {
	protected Metric metric;
	protected double[][] data;
	ArrayList<Integer> objectsIndices;

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

	/**
	 * Hyperplane partitioning. Subdivide the list of points by a plane
	 * perpendicular to the axis corresponding to 'cutDimension' at value 'cutValue'.
	 * On return:
	 * data[objectsIndices.get(start..start+lim1-1)][cutDimension] < cutValue
	 * data[objectsIndices.get(start+lim1..start+lim2-1)][cutDimension] == cutValue
	 * data[objectsIndices.get(start+lim2..start + count - 1)][cutDimension] > cutValue
	 */
	protected void planeSplit (int start, int count, int cutDimension, double cutValue, int[] lim1Andlim2Wrapper) {
		int lim1, lim2;
		int left = start;
		int right = start + count - 1;
		for (;;) {
			while (left <= right && data[objectsIndices.get(left)][cutDimension] < cutValue)
				left++;
			while (left <= right && data[objectsIndices.get(right)][cutDimension] >= cutValue)
				right--; 
			if (left > right)
				break;
			Collections.swap(objectsIndices, left, right);
			left++;
			right--;
		}
		
		lim1 = left;
		right = start + count - 1;
		for (;;) {
			while (left <= right && data[objectsIndices.get(left)][cutDimension] <= cutValue)
				left++;
			while (left <= right && data[objectsIndices.get(right)][cutDimension] > cutValue)
				right--;
			if (left > right)
				break;
			Collections.swap(objectsIndices, left, right);
			left++;
			right--;
		}
		lim2 = left;
		
		// Convert 'lim1' and 'lim2' from absolute index values in 'objectsIndices',
		// to relative to 'start'.
		lim1 -= start;
		lim2 -= start;
		
		lim1Andlim2Wrapper[0] = lim1;
		lim1Andlim2Wrapper[1] = lim2;
	}
}
