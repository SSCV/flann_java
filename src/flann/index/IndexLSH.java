package flann.index;

import java.util.ArrayList;

import flann.metric.Metric;
import flann.result_set.ResultSet;

public class IndexLSH extends IndexBase {
	ArrayList<LSHTable> tables;
	int tablesNumber;
	int keySize;

	// How far should we look for neighbors in multi-probe LSH.
	int multiProbeLevel;

	// The XOR masks to apply to a key to get the neighboring buckets.
	ArrayList<Integer> xorMasks;

	public static class BuildParams {
		public int tablesNumber;
		public int keySize;
		public int multiProbeLevel;

		public BuildParams() {
			this.tablesNumber = 12;
			this.keySize = 20;
			this.multiProbeLevel = 2;
		}

		public BuildParams(int tablesNumber, int keySize, int multiProbeLevel) {
			this.tablesNumber = tablesNumber;
			this.keySize = keySize;
			this.multiProbeLevel = multiProbeLevel;
		}
	}

	public static class SearchParams extends SearchParamsBase {
	}

	public IndexLSH(Metric metric, double[][] data) {
		super(metric, data);
	}

	@Override
	protected void buildIndexImpl() {
	}

	@Override
	protected void findNeighbors(ResultSet resultSet, double[] query,
			SearchParamsBase searchParams) {

	}
}
