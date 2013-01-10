package flann.index;

import java.util.ArrayList;
import java.util.Collections;
import java.util.PriorityQueue;

import flann.exception.ExceptionFLANN;
import flann.metric.Metric;
import flann.result_set.ResultSet;
import flann.util.Utils;

public class IndexKMeans extends IndexBase {
	// Branching factor.
	int branching;

	// Max iterations to perform.
	int iterations;

	// Algorithm used for picking the initial cluster centers.
	CentersChooserAlgorithm centersInit;

	float cbIndex;
	Node root;

	public enum CentersChooserAlgorithm {
		FLANN_CENTERS_RANDOM, FLANN_CENTERS_GONZALES, FLANN_CENTERS_KMEANSPP
	}

	public static class BuildParams {
		public int branching;
		public int iterations;
		public CentersChooserAlgorithm centersInit;
		public float cbIndex;

		public BuildParams() {
			this.branching = 32;
			this.iterations = 11;
			this.centersInit = CentersChooserAlgorithm.FLANN_CENTERS_RANDOM;
			this.cbIndex = 0.2f;
		}

		public BuildParams(int branching, int iterations,
				CentersChooserAlgorithm centersInit, float cbIndex) {
			this.branching = branching;
			this.iterations = iterations;
			this.centersInit = centersInit;
			this.cbIndex = cbIndex;
		}
	}

	public static class SearchParams extends SearchParamsBase {
	}

	private class Node {
		// The cluster center.
		public double[] pivot;
		// The cluster radius.
		public double radius;
		// The cluster variance.
		public double variance;
		// Number of points in the cluster.
		public int size;
		// Children nodes (only for non-terminal nodes).
		public ArrayList<Node> children;
		// Node points (only for terminal nodes).
		public ArrayList<PointInfo> points;

		public Node() {
			this.children = new ArrayList<Node>();
			this.points = new ArrayList<PointInfo>();
		}
	}

	private class PointInfo {
		int index;
		double[] point;
	}

	public IndexKMeans(Metric metric, double[][] data, BuildParams buildParams) {
		super(metric, data);
		root = null;

		// Get parameters.
		branching = buildParams.branching;
		iterations = buildParams.iterations;
		centersInit = buildParams.centersInit;
		cbIndex = buildParams.cbIndex;
		if (iterations < 0) {
			iterations = Integer.MAX_VALUE;
		}
	}

	@Override
	protected void buildIndexImpl() throws ExceptionFLANN {
		if (branching < 2) {
			throw new ExceptionFLANN("Branching factor must be at least 2");
		}

		// Prepare objectsIndices.
		objectsIndices = new ArrayList<Integer>();
		for (int i = 0; i < numberOfObjects; i++) {
			objectsIndices.add(i);
		}

		root = new Node();
		computeNodeStatistics(root, objectsIndices);
		computeClustering(root, 0, numberOfObjects, branching);
	}

	/**
	 * Clustering function that takes a cut in the hierarchical k-means tree and
	 * returns the clusters centers of that clustering.
	 */
	public int getClusterCenters(double[][] centers) throws ExceptionFLANN {
		int numClusters = centers.length;
		if (numClusters < 1) {
			throw new ExceptionFLANN("Number of clusters must be at least 1.");
		}
		double[] variance = new double[1];
		Node[] clusters = new Node[numClusters];
		int clusterCount = getMinVarianceClusters(root, clusters, numClusters,
				variance);
		for (int i = 0; i < clusterCount; i++) {
			double[] center = clusters[i].pivot;
			for (int j = 0; j < numberOfDimensions; j++) {
				centers[i][j] = center[j];
			}
		}
		return clusterCount;
	}

	@Override
	protected void findNeighbors(ResultSet resultSet, double[] query,
			SearchParamsBase searchParams) {
		int maxChecks = searchParams.checks;
		if (maxChecks == -1) {
			findExactNN(root, resultSet, query);
		} else {
			// Priority queue storing intermediate branches in the
			// best-bin-first search.
			PriorityQueue<Branch<Node>> heap = new PriorityQueue<Branch<Node>>(
					numberOfObjects);

			int checks[] = new int[1];
			checks[0] = 0;
			findNN(root, resultSet, query, checks, maxChecks, heap);

			Branch<Node> branch;
			while ((branch = heap.poll()) != null
					&& (checks[0] < maxChecks || !resultSet.full())) {
				findNN(branch.node, resultSet, query, checks, maxChecks, heap);
			}
		}
	}

	/**
	 * Helper function the descends in the hierarchical k-means tree by spliting
	 * those clusters that minimize the overall variance of the clustering.
	 */
	private int getMinVarianceClusters(Node root, Node[] clusters,
			int clustersLength, double[] varianceValue) {
		int clusterCount = 1;
		clusters[0] = root;
		double meanVariance = root.variance * root.size;

		while (clusterCount < clustersLength) {
			double minVariance = Double.MAX_VALUE;
			int splitIndex = -1;

			for (int i = 0; i < clusterCount; i++) {
				if (!clusters[i].children.isEmpty()) {
					double variance = meanVariance - clusters[i].variance
							* clusters[i].size;
					for (int j = 0; j < branching; j++) {
						variance += clusters[i].children.get(j).variance
								* clusters[i].children.get(j).size;
					}
					if (variance < minVariance) {
						minVariance = variance;
						splitIndex = i;
					}
				}
			}

			if (splitIndex == -1)
				break;

			if (branching + clusterCount - 1 > clustersLength)
				break;

			meanVariance = minVariance;

			// Split node.
			Node toSplit = clusters[splitIndex];
			clusters[splitIndex] = toSplit.children.get(0);
			for (int i = 1; i < branching; i++) {
				clusters[clusterCount++] = toSplit.children.get(i);
			}
		}

		varianceValue[0] = meanVariance / root.size;
		return clusterCount;
	}

	/**
	 * Method that computes the squared distance from the query point q from
	 * inside region with center c to the border between this region and the
	 * region with center p.
	 */
	public double getDistanceToBorder(double[] p, double[] c, double[] q) {
		double sum = 0;
		double sum2 = 0;

		for (int i = 0; i < numberOfDimensions; i++) {
			double t = c[i] - p[i];
			sum += t * (q[i] - (c[i] + p[i]) / 2);
			sum2 += t * t;
		}

		return sum * sum / sum2;
	}

	/**
	 * It computes the order in which to traverse the child nodes of a
	 * particular node.
	 */
	private void getCenterOrdering(Node node, double[] q, int[] sortIndices) {
		double[] domainDistances = new double[branching];
		for (int i = 0; i < branching; ++i) {
			double dist = metric.distance(q, node.children.get(i).pivot);

			int j = 0;
			while (domainDistances[j] < dist && j < i)
				j++;
			for (int k = i; k > j; k--) {
				domainDistances[k] = domainDistances[k - 1];
				sortIndices[k] = sortIndices[k - 1];
			}

			domainDistances[j] = dist;
			sortIndices[j] = i;
		}
	}

	// Function the performs exact nearest neighbor search by traversing the
	// entire tree.
	private void findExactNN(Node node, ResultSet resultSet, double[] query) {
		// Ignore those clusters that are too far away.
		{
			double bsq = metric.distance(query, node.pivot);
			double rsq = node.radius;
			double wsq = resultSet.worstDistance();

			double val = bsq - rsq - wsq;
			double val2 = val * val - 4 * rsq * wsq;

			if (val > 0 && val2 > 0) {
				return;
			}
		}

		if (node.children.isEmpty()) {
			for (int i = 0; i < node.size; i++) {
				PointInfo pointInfo = node.points.get(i);
				int index = pointInfo.index;
				double dist = metric.distance(pointInfo.point, query);
				resultSet.addPoint(dist, index);
			}
		} else {
			int[] sortIndices = new int[branching];
			getCenterOrdering(node, query, sortIndices);
			for (int i = 0; i < branching; i++) {
				findExactNN(node.children.get(sortIndices[i]), resultSet, query);
			}
		}
	}

	/**
	 * Helper function that computes the nearest children of a node to a given
	 * query point.
	 */
	private int exploreNodeBranches(Node node, double[] q,
			PriorityQueue<Branch<Node>> heap) {
		double[] domainDistances = new double[branching];
		int bestIndex = 0;

		double d = metric.distance(q, node.children.get(bestIndex).pivot);
		domainDistances[bestIndex] = d;
		for (int i = 1; i < branching; i++) {
			domainDistances[i] = metric.distance(q, node.children.get(i).pivot);
			if (domainDistances[i] < domainDistances[bestIndex]) {
				bestIndex = i;
			}
		}

		for (int i = 0; i < branching; i++) {
			if (i != bestIndex) {
				domainDistances[i] -= cbIndex * node.children.get(i).variance;
				heap.add(new Branch<Node>(node.children.get(i),
						domainDistances[i]));
			}
		}

		return bestIndex;
	}

	private void findNN(Node node, ResultSet resultSet, double[] query,
			int[] checks, int maxChecks, PriorityQueue<Branch<Node>> heap) {
		// Ignore those clusters that are too far away
		{
			double bsq = metric.distance(query, node.pivot);
			double rsq = node.radius;
			double wsq = resultSet.worstDistance();

			double val = bsq - rsq - wsq;
			double val2 = val * val - 4 * rsq * wsq;

			if (val > 0 && val2 > 0) {
				return;
			}
		}

		if (node.children.isEmpty()) {
			if (checks[0] >= maxChecks) {
				if (resultSet.full())
					return;
			}

			for (int i = 0; i < node.size; i++) {
				PointInfo pointInfo = node.points.get(i);
				int index = pointInfo.index;
				double dist = metric.distance(pointInfo.point, query);
				resultSet.addPoint(dist, index);
				checks[0]++;
			}
		} else {
			int closestCenter = exploreNodeBranches(node, query, heap);
			findNN(node.children.get(closestCenter), resultSet, query, checks,
					maxChecks, heap);
		}
	}

	private void computeNodeStatistics(Node node, ArrayList<Integer> indices) {
		int size = indices.size();

		double[] mean = new double[numberOfDimensions];
		for (int i = 0; i < numberOfDimensions; i++) {
			mean[i] = 0;
		}

		for (int i = 0; i < size; i++) {
			double[] query = data[indices.get(i)];
			for (int j = 0; j < numberOfDimensions; j++) {
				mean[j] += query[j];
			}
		}

		double divFactor = 1.0 / size;
		for (int j = 0; j < numberOfDimensions; j++) {
			mean[j] *= divFactor;
		}

		double radius = 0;
		double variance = 0;
		for (int i = 0; i < size; i++) {
			double dist = metric.distance(mean, data[indices.get(i)]);
			if (dist > radius) {
				radius = dist;
			}
			variance += dist;
		}
		variance /= size;

		node.variance = variance;
		node.radius = radius;
		node.pivot = mean;
	}

	private void computeClustering(Node node, int start, int count,
			int branching) {
		node.size = count;

		if (count < branching) {
			for (int i = 0; i < count; i++) {
				node.points.add(new PointInfo());
				node.points.get(i).index = objectsIndices.get(start + i);
				node.points.get(i).point = data[objectsIndices.get(start + i)];
			}
			node.children.clear();
			return;
		}

		ArrayList<Integer> centersIdx = new ArrayList<Integer>(branching);

		switch (centersInit) {
		case FLANN_CENTERS_RANDOM:
			CenterChooser.Random(metric, data, branching, objectsIndices,
					start, count, centersIdx);
			break;
		case FLANN_CENTERS_GONZALES:
			CenterChooser.Gonzales(metric, data, branching, objectsIndices,
					start, count, centersIdx);
			break;
		case FLANN_CENTERS_KMEANSPP:
			CenterChooser.KMeansPP(metric, data, branching, objectsIndices,
					start, count, centersIdx);
			break;
		}

		int centersLength = centersIdx.size();
		if (centersLength < branching) {
			for (int i = 0; i < count; i++) {
				node.points.get(start + i).index = objectsIndices
						.get(start + i);
				node.points.get(start + i).point = data[objectsIndices
						.get(start + i)];
			}
			node.children.clear();
			return;
		}

		double[][] dcenters = new double[branching][numberOfDimensions];
		for (int i = 0; i < centersLength; i++) {
			double[] vec = data[centersIdx.get(i)];
			for (int k = 0; k < numberOfDimensions; k++) {
				dcenters[i][k] = vec[k];
			}
		}

		double[] radiuses = new double[branching];
		int[] count2 = new int[branching];
		for (int i = 0; i < branching; i++) {
			radiuses[i] = 0;
			count2[i] = 0;
		}
		int[] belongsTo = new int[count];

		for (int i = 0; i < count; i++) {
			double sqDist = metric.distance(
					data[objectsIndices.get(start + i)], dcenters[0]);
			belongsTo[i] = 0;
			for (int j = 1; j < branching; j++) {
				double newSqDist = metric.distance(
						data[objectsIndices.get(start + i)], dcenters[j]);
				if (sqDist > newSqDist) {
					belongsTo[start + i] = j;
					sqDist = newSqDist;
				}
			}
			if (sqDist > radiuses[belongsTo[start + i]]) {
				radiuses[belongsTo[start + i]] = sqDist;
			}
			count2[belongsTo[start + i]] = count2[start + i];
		}

		boolean converged = false;
		int iteration = 0;
		while (!converged && iteration < iterations) {
			converged = true;
			iteration++;

			// Compute the new cluster centers.
			for (int i = 0; i < branching; i++) {
				for (int j = 0; j < numberOfDimensions; j++) {
					dcenters[i][j] = 0;
				}
				radiuses[i] = 0.0;
			}

			for (int i = 0; i < count; i++) {
				double[] vec = data[objectsIndices.get(start + i)];
				double[] center = dcenters[belongsTo[start + i]];
				for (int k = 0; k < numberOfDimensions; k++) {
					center[k] += vec[k];
				}
			}

			for (int i = 0; i < branching; i++) {
				int cnt = count2[i];
				double divFactor = 1.0 / cnt;
				for (int k = 0; k < numberOfDimensions; k++) {
					dcenters[i][k] *= divFactor;
				}
			}

			for (int i = 0; i < count; i++) {
				double sqDist = metric.distance(
						data[objectsIndices.get(start + i)], dcenters[0]);
				int newCentroid = 0;
				for (int j = 1; j < branching; j++) {
					double newSqDist = metric.distance(
							data[objectsIndices.get(start + i)], dcenters[j]);
					if (sqDist > newSqDist) {
						newCentroid = j;
						sqDist = newSqDist;
					}
				}

				if (sqDist > radiuses[newCentroid]) {
					radiuses[newCentroid] = sqDist;
				}

				if (newCentroid != belongsTo[start + i]) {
					count2[belongsTo[start + i]]--;
					count2[newCentroid]++;
					belongsTo[start + i] = newCentroid;
					converged = false;
				}
			}

			for (int i = 0; i < branching; i++) {
				// If one cluster converges to an empty cluster,
				// move an element into that cluster.
				if (count2[i] == 0) {
					int j = (i + 1) % branching;
					while (count2[j] <= 1) {
						j = (j + 1) % branching;
					}

					for (int k = 0; k < count; k++) {
						if (belongsTo[start + k] == j) {
							belongsTo[start + k] = i;
							count2[j]--;
							count2[i]++;
							break;
						}
					}

					converged = false;
				}
			}
		}

		double[][] centers = new double[branching][numberOfDimensions];
		for (int i = 0; i < branching; i++) {
			for (int j = 0; j < numberOfDimensions; j++) {
				centers[i][j] = dcenters[i][j];
			}
		}

		// Compute kmeans clustering for each of the resulting clusters.
		int start2 = 0;
		int end = start2;
		for (int c = 0; c < branching; c++) {
			int s = count2[c];

			double variance = 0;
			for (int i = 0; i < count; i++) {
				if (belongsTo[start + i] == c) {
					variance += metric.distance(centers[c],
							data[objectsIndices.get(start + i)]);
					Collections.swap(objectsIndices, i, end);
					Utils.swapArray(belongsTo, i, end);
					end++;
				}
			}
			variance /= s;

			Node newNode = new Node();
			newNode.radius = radiuses[c];
			newNode.pivot = centers[c];
			newNode.variance = variance;
			node.children.add(newNode);
			computeClustering(node.children.get(c), start + start2, end
					- start2, branching);
			start2 = end;
		}
	}
}