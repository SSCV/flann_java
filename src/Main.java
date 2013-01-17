import flann.index.IndexBase;
import flann.index.IndexKDTree;
import flann.index.IndexKDTreeSingle;
import flann.index.IndexKMeans;
import flann.index.IndexLSH;
import flann.metric.Metric;
import flann.metric.MetricEuclideanSquared;
import flann.metric.MetricHamming;

public class Main {
	public static void main(String[] args) {
		// Each row is a 2D point.
		double[][] data = { { 1, 1 }, { 3, 3 }, { 3, 4 }, { 7, 7 }, { 7, 6 } };
		int[][] dataBinary = { { 1 }, { 3 }, { 4 }, { 7 }, { 8 } };

		// Each row is a query.
		double[][] queries = { { 3, 1 }, { 5, 3 }, { 3, 3.6 } };
		int[][] queriesBinary = { { 2 }, { 5 }, { 6 } };

		// kNN search parameter.
		int k = 3;

		// The output is stored in two matrices 'indices' and 'distances',
		// such that each row corresponds to a single query. Both matrices
		// have 'k' columns for the k indices and distances to the k nearest
		// neighbors.
		int q = queries.length;
		int[][] indices = new int[q][k];
		double[][] distances = new double[q][k];

		// Specify metric that will be used to compute distances between
		// objects.
		Metric metric = new MetricEuclideanSquared();

		// Construct/initialize the index, and then build.
		// IndexKDTreeSingle.BuildParams buildParams = new
		// IndexKDTreeSingle.BuildParams (3, false);
		// IndexBase index = new IndexKDTreeSingle (metric, data, buildParams);
		// index.buildIndex();

		IndexKDTreeSingle.SearchParams searchParams = new IndexKDTreeSingle.SearchParams();

		// Perform kNN search.
		searchParams.eps = 0.0f;
		searchParams.maxNeighbors = k;
		// index.knnSearch (queries, indices, distances, searchParams);

		// Perform radius search.
		searchParams.radius = 0.4;
		// index.radiusSearch (queries, indices, distances, searchParams);

		// Construct and search with randomized kd-trees.
		IndexKDTree.BuildParams buildParams = new IndexKDTree.BuildParams(4);
		IndexBase index2 = new IndexKDTree(metric, data, buildParams);
		index2.buildIndex();
		IndexKDTree.SearchParams searchParams2 = new IndexKDTree.SearchParams();
		searchParams2.eps = 0.0f;
		searchParams2.maxNeighbors = k;
		searchParams2.checks = 128;
		// index2.knnSearch(queries, indices, distances, searchParams2);

		// Construct and search with IndexKMeans.
		IndexKMeans.BuildParams buildParams2 = new IndexKMeans.BuildParams();
		IndexBase index3 = new IndexKMeans(metric, data, buildParams2);
		index3.buildIndex();
		IndexKMeans.SearchParams searchParams3 = new IndexKMeans.SearchParams();
		searchParams3.maxNeighbors = k;
		searchParams3.eps = 0.0f;
		// index3.knnSearch(queries, indices, distances, searchParams3);

		// Construct and search with LSH index.
		Metric metric2 = new MetricHamming();
		IndexLSH.BuildParams buildParams3 = new IndexLSH.BuildParams();
		IndexBase index4 = new IndexLSH(metric2, dataBinary, buildParams3);
		index4.buildIndex();
		IndexLSH.SearchParams searchParams4 = new IndexLSH.SearchParams();
		searchParams4.maxNeighbors = k;
		searchParams4.eps = 0.0f;
		index4.knnSearch(queriesBinary, indices, distances, searchParams4);

		// Print the result contained in matrices 'indices' and 'distances'.
		for (int i = 0; i < q; i++) {
			for (int j = 0; j < k; j++) {
				System.out.print("(" + indices[i][j] + ", " + distances[i][j]
						+ ")  ");
			}
			System.out.println();
		}
	}
}