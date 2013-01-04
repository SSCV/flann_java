import flann.metric.*;
import flann.index.*;



public class Main {
	public static void main(String[] args) {
		// Each row is a 2D point.
		double[][] data = {
				{1,1},
				{3,3},
				{3,4},
				{7,7},
				{7,6}
		};

		// Each row is a query.
		double[][] queries = {
				{3,1},
				{5,3}
		};

		// kNN search parameter.
		int k = 3;

		// The output is stored in two matrices 'indices' and 'distances',
		// such that each row corresponds to a single query. Both matrices
		// have 'k' columns for the k indices and distances to the k nearest
		// neighbors.
		int q = queries.length;
		int[][] indices = new int[q][k];
		double[][] distances = new double[q][k];

		// Specify metric that will be used to compute distances between objects.
		Metric metric = new MetricEuclideanSquared();

		// Construct/initialize the index, and then build.
		IndexKDTreeSingle.BuildParams buildParams = new IndexKDTreeSingle.BuildParams (3, false);
		IndexKDTreeSingle index = new IndexKDTreeSingle (metric, data, buildParams);
		index.buildIndex();

		// Perform kNN search.
		IndexKDTreeSingle.SearchParams searchParams = new IndexKDTreeSingle.SearchParams ();
		searchParams.eps = 0.0f;
		searchParams.numberOfNeighbors = k;
		index.knnSearch (queries, indices, distances, searchParams);

		// Print the result contained in matrices 'indices' and 'distances'.
		for (int i = 0; i < q; i++) {
			for (int j = 0; j < k; j++) {
				System.out.print ("(" + indices[i][j] + ", " + distances[i][j] + ")  ");
			}
			System.out.println ();
		}
	}
}