package flann.index;

import flann.metric.Metric;
import flann.result_set.KNNSimpleResultSet;
import flann.util.BoundingBox;

import java.util.ArrayList;
import java.util.Collections;



public class IndexKDTreeSingle extends IndexBase {
	private int maxPointsInOneLeafNode;
	private Node root;
	private BoundingBox rootBBox;
	ArrayList<Integer> objectsIndices;
	
	private class Node {
		public Node child1, child2;
		public int cutDimension;
		public double cutDimensionLow, cutDimensionHigh;

		// Indices of points contained in this node (if it is a leaf node).
		public int leftObjectIndex, rightObjectIndex;
	}

    public IndexKDTreeSingle (Metric metric, double[][] data, int maxPointsInOneLeafNode) {
		super(metric, data);

		this.maxPointsInOneLeafNode = maxPointsInOneLeafNode;
		this.root = null;
		this.rootBBox = new BoundingBox();

		objectsIndices = new ArrayList<Integer>();
		for (int i = 0; i < numberOfObjects; i++) {
			objectsIndices.add(i);
		}
	}

	public void buildIndex() {
		rootBBox.fitToData (data);
		root = divideTree (0, numberOfObjects, rootBBox);
	}

	public void knnSearch (double[][] queries, int[][] indices, double[][] distances, int k, float eps) {
		KNNSimpleResultSet resultSet = new KNNSimpleResultSet (k);
		for (int i = 0; i < queries.length; i++) {
			resultSet.clear();
			findNeighbors (resultSet, queries[i], k, eps);
			int n = Math.min (resultSet.size(), k);
			resultSet.copy (distances[i], indices[i], n);
		}
    }

	private void findNeighbors (KNNSimpleResultSet resultSet, double[] query, int k, float eps) {
		float epsError = 1 + eps;
		ArrayList<Double> distances = new ArrayList<Double> ();
		for (int i = 0; i < numberOfDimensions; i++) {
			distances.add (0.0);
		}
		double distsq = rootBBox.getBoxPointDistancesPerDimension(query, distances, metric);
		searchLevel (root, query, resultSet, k, epsError, distsq, distances);
	}

	private void searchLevel(Node node, double[] query, KNNSimpleResultSet resultSet,
							 int k, float eps, Double mindistsq, ArrayList<Double> dists) {
		// If this is a leaf node.
		if (node.child1 == null && node.child2 == null) {
			double worstDistance = resultSet.worstDistance();
			for (int i = node.leftObjectIndex; i < node.rightObjectIndex; i++) {
				double dist = metric.distance (query, data[objectsIndices.get(i)]);
                if (dist < worstDistance) {
                    resultSet.addPoint (dist, objectsIndices.get(i));
                }
			}
			return;
		}
		
		// Which child branch should be taken first?
		int cutDimension = node.cutDimension;
		double queryValueInCutDimension = query[cutDimension];
		double diff1 = queryValueInCutDimension - node.cutDimensionLow;
		double diff2 = queryValueInCutDimension - node.cutDimensionHigh;
		Node bestChild, otherChild;
		double cutDistance;
		if (diff1 + diff2 < 0) {
			bestChild = node.child1;
            otherChild = node.child2;
            cutDistance = metric.distance (queryValueInCutDimension, node.cutDimensionHigh);
		} else {
			bestChild = node.child2;
            otherChild = node.child1;
            cutDistance = metric.distance (queryValueInCutDimension, node.cutDimensionLow);
		}
		
		// Call recursively to search next level down.
		searchLevel(bestChild, query, resultSet, k, eps, mindistsq, dists);
		double dst = dists.get (cutDimension);
		mindistsq = mindistsq - dst + cutDistance;
		dists.set (cutDimension, cutDistance);
        if (mindistsq*eps <= resultSet.worstDistance()) {
        	searchLevel(otherChild, query, resultSet, k, eps, mindistsq, dists);
        }
        dists.set (cutDimension, dst);
	}

	/**
	 * Construct a tree node that subdivides the list of objects/points from
	 * objectsIndices.get(left) to objectsIndices.get(right). The routine is called
	 * recursively on each sublist.
	 */
	private Node divideTree (int left, int right, BoundingBox bbox) {
		Node node = new Node ();
		
		// If too few objects remain, then make this a leaf node.
		if (right-left <= maxPointsInOneLeafNode) {
			node.child1 = node.child2 = null;
			node.leftObjectIndex = left;
			node.rightObjectIndex = right;

			// Compute bounding box for the objects in this leaf node.
			for (int i = 0; i < numberOfDimensions; i++) {
				double val = data[objectsIndices.get(left)][i];
				bbox.add (val, val);
			}			
			for (int k = left+1; k < right; k++) {
				for (int i = 0; i < numberOfDimensions; i++) {
					double objectKInDimensionI = data[objectsIndices.get(k)][i];
					if (objectKInDimensionI < bbox.getMin (i))
						bbox.setMin (i, objectKInDimensionI);
					if (objectKInDimensionI > bbox.getMax (i))
						bbox.setMax (i, objectKInDimensionI);
				}
			}
		} else {
			middleSplitResult out = new middleSplitResult();
            middleSplit (left, right-left, bbox, out);
            int cutObjectIndex = out.cutObjectIndex;
            int cutDimension = out.cutDimension;
    		double cutValue = out.cutValue;

            node.cutDimension = cutDimension;

            BoundingBox leftBBox = new BoundingBox (bbox);
            leftBBox.setMax (cutDimension, cutValue);
            node.child1 = divideTree (left, left+cutObjectIndex, leftBBox);

            BoundingBox rightBBox = new BoundingBox (bbox);
            rightBBox.setMin (cutDimension, cutValue);
            node.child2 = divideTree (left+cutObjectIndex, right, rightBBox);

            node.cutDimensionLow = leftBBox.getMax (cutDimension);
            node.cutDimensionHigh = rightBBox.getMin (cutDimension);

            for (int i=0; i < numberOfDimensions; ++i) {
            	bbox.setMin (i, Math.min (leftBBox.getMin (i), rightBBox.getMin (i)));
            	bbox.setMax (i, Math.max (leftBBox.getMax (i), rightBBox.getMax (i)));
            }
		}

		return node;
	}

	private class middleSplitResult {
		public int cutObjectIndex;
		public int cutDimension;
		public double cutValue;
	}

	private void middleSplit (int start, int count, BoundingBox bbox, middleSplitResult out) {
		// Compute 'cutDimension' and the appropriate 'cutValue'.
		out.cutDimension = bbox.getMaxSpanDimension();
		double min = bbox.getMin (out.cutDimension);
		double max = bbox.getMax (out.cutDimension);
		out.cutValue = (min + max) / 2;

		// Hyperplane partitioning.
		int[] lim1Andlim2Wrapper = new int[2];
		planeSplit (start, count, out.cutDimension, out.cutValue, lim1Andlim2Wrapper);
		int lim1 = lim1Andlim2Wrapper[0];
		int lim2 = lim1Andlim2Wrapper[1];

		// Choose the object through which the hyperplane goes through.
		int countHalf = count / 2;
		if (lim1 > countHalf)
			out.cutObjectIndex = lim1;
        else if (lim2 < countHalf)
        	out.cutObjectIndex = lim2;
        else
        	out.cutObjectIndex = countHalf;
	}

	/**
	 * Hyperplane partitioning. Subdivide the list of points by a plane
	 * perpendicular to the axis corresponding to 'cutDimension' at value 'cutValue'.
	 * On return:
	 * data[objectsIndices.get(start..start+lim1-1)][cutDimension] < cutValue
	 * data[objectsIndices.get(start+lim1..start+lim2-1)][cutDimension] == cutValue
	 * data[objectsIndices.get(start+lim2..start + count - 1)][cutDimension] > cutValue
	 */
	private void planeSplit (int start, int count, int cutDimension, double cutValue, int[] lim1Andlim2Wrapper) {
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