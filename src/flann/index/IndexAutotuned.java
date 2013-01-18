package flann.index;

import java.util.ArrayList;
import java.util.Map;

import org.apache.commons.lang3.time.StopWatch;

import flann.metric.Metric;
import flann.result_set.ResultSet;

public class IndexAutotuned extends IndexBase {
	float targetPrecision;
	float buildWeight;
	float memoryWeight;
	float sampleFraction;

	float speedup;

	IndexBase bestIndex;

	double[][] sampledData;
	double[][] testData;

	Map<String, Object> bestBuildParams;
	Map<String, Object> bestSearchParams;

	int[][] gtMatches;

	public static class BuildParams {
		public float targetPrecision;
		public float buildWeight;
		public float memoryWeight;
		public float sampleFraction;

		public BuildParams() {
			this.targetPrecision = 0.8f;
			this.buildWeight = 0.01f;
			this.memoryWeight = 0.0f;
			this.sampleFraction = 0.1f;
		}

		public BuildParams(float targetPrecision, float buildWeight,
				float memoryWeight, float sampleFraction) {
			this.targetPrecision = targetPrecision;
			this.buildWeight = buildWeight;
			this.memoryWeight = memoryWeight;
			this.sampleFraction = sampleFraction;
		}
	}

	public static class SearchParams extends SearchParamsBase {
	}

	public IndexAutotuned(Metric metric, double[][] data,
			BuildParams buildParams) {
		super(metric, data);

		this.targetPrecision = buildParams.targetPrecision;
		this.buildWeight = buildParams.buildWeight;
		this.memoryWeight = buildParams.memoryWeight;
		this.sampleFraction = buildParams.sampleFraction;

		this.speedup = 0;
		this.bestIndex = null;
	}

	private class CostData {
		float searchTimeCost;
		float buildTimeCost;
		float memoryCost;
		float totalCost;
		Map<String, Object> params;
	}

	private void evaluateKMeans(CostData cost) {
		StopWatch timer = new StopWatch();
		int checks;
		int nn = 1;

		int iterations = (Integer) cost.params.get("iterations");
		int branching = (Integer) cost.params.get("branching");
		CenterChooser.Algorithm centersInit = (CenterChooser.Algorithm) cost.params
				.get("centersInit");
		System.out.printf(
				"KMeansTree using params: max_iterations=%d, branching=%d\n",
				iterations, branching);

		IndexKMeans.BuildParams buildParams = new IndexKMeans.BuildParams();
		buildParams.branching = branching;
		buildParams.iterations = iterations;
		buildParams.centersInit = centersInit;
		IndexKMeans kmeans = new IndexKMeans(metric, sampledData, buildParams);

		// Measure index build time.
		timer.start();
		kmeans.buildIndex();
		timer.stop();
		float buildTime = timer.getTime() / 1000.0f;

		// Measure search time.
		float searchTime = testIndexPrecision();

		int sampledDataRows = sampledData.length;
		int sampledDataCols = sampledData[0].length;
		float datasetMemory = sampledDataRows * sampledDataCols * 4;
		cost.memoryCost = (1 + datasetMemory) / datasetMemory;
		cost.searchTimeCost = searchTime;
		cost.buildTimeCost = buildTime;

		System.out.printf(
				"KMeansTree buildTime=%f, searchTime=%f, build_weight=%f\n",
				buildTime, searchTime, buildWeight);
	}

	private void evaluateKDTree(CostData cost) {
		StopWatch timer = new StopWatch();
		int checks;
		int nn = 1;

		int trees = (Integer) cost.params.get("trees");
		System.out.printf("KDTree using params: trees=%d\n", trees);

		IndexKDTree.BuildParams buildParams = new IndexKDTree.BuildParams();
		buildParams.trees = trees;
		IndexKDTree kdtree = new IndexKDTree(metric, sampledData, buildParams);

		// Measure index build time.
		timer.start();
		kdtree.buildIndex();
		timer.stop();
		float buildTime = timer.getTime() / 1000.0f;

		// Measure search time.
		float searchTime = testIndexPrecision();

		int sampledDataRows = sampledData.length;
		int sampledDataCols = sampledData[0].length;
		float datasetMemory = sampledDataRows * sampledDataCols * 4;
		cost.memoryCost = (1 + datasetMemory) / datasetMemory;
		cost.searchTimeCost = searchTime;
		cost.buildTimeCost = buildTime;

		System.out.printf("KDTree buildTime=%f, searchTime=%f\n", buildTime,
				searchTime);
	}

	private void optimizeKDTree(ArrayList<CostData> costs) {
		System.out.println("KD-TREE, Step 1: Exploring parameter space");

		// Explor KD-Tree parameters space using the parameters below.
		int testTrees[] = { 1, 4, 8, 16, 32 };

		// Evaluate KD-Tree for all parameter combinations.
		for (int i = 0; i < testTrees.length; i++) {
			CostData cost = new CostData();
			cost.params.put("trees", testTrees[i]);
			evaluateKDTree(cost);
			costs.add(cost);
		}
	}

	private void optimizeKMeans(ArrayList<CostData> costs) {
		System.out.println("KMEANS, Step 1: Exploring parameter space");

		// Explore K-Means parameters space using combinations of parameters.
		int[] maxIterations = { 1, 5, 10, 15 };
		int[] branchingFactors = { 16, 32, 64, 128, 256 };

		int kMeansParamSpaceSize = maxIterations.length
				* branchingFactors.length;
		// Evaluate K-Means for all parameters combinations.
		for (int i = 0; i < maxIterations.length; i++) {
			for (int j = 0; j < branchingFactors.length; j++) {
				CostData cost = new CostData();

				cost.params.put("centersInit",
						CenterChooser.Algorithm.FLANN_CENTERS_RANDOM);
				cost.params.put("iterations", maxIterations[i]);
				cost.params.put("branching", branchingFactors[i]);

				evaluateKMeans(cost);
				costs.add(cost);
			}
		}
	}

	@Override
	private void knnSearch(double[][] queries, int[][] indices,
			double[][] distances, SearchParamsBase searchParams) {
		if (searchParams.checks == 1) {
			// bestIndex.knnSearch(queries, indices, distances, searchParams);
			// bestSearchParams_);
		} else {
			// bestIndex_->knnSearch(queries, indices, dists, knn, params);
		}
	}

	@Override
	private int radiusSearch(double[][] queries, int[][] indices,
			double[][] distances, SearchParamsBase searchParams) {
		if (searchParams.checks == 1) {
			// return bestIndex_->radiusSearch(queries, indices, dists, radius,
			// bestSearchParams_);
		} else {
			// return bestIndex_->radiusSearch(queries, indices, dists, radius,
			// params);
		}
	}

	private void estimateBuildParams() {
		ArrayList<CostData> costs = new ArrayList<CostData>();

		int sampleSize = (int) (sampleFraction * data.length);
		int testSampleSize = Math.min(sampleSize / 10, 1000);

		System.out
				.printf("Entering autotuning, dataset size: %d, sampleSize: %d, testSampleSize: %d, target precision: %f\n",
						data.length, sampleSize, testSampleSize,
						targetPrecision);

		// For a very small dataset, it makes no sense to build any fancy index,
		// just use linear search.
		if (testSampleSize < 10) {
			System.out.println("Choosing linear, dataset too small");
			// return LinearIndexParams();
		}

		// We use a fraction of the original dataset to speedup the autotune
		// algorithm.
		sampledData = randomSample(data, sampleSize);
		testData = randomSample(sampledData, testSampleSize, true);

		// We compute the ground truth using linear search.
		System.out.println("Computing ground truth...");
		gtMatches = new int[testData.length][1];

		StopWatch timer = new StopWatch();
		int repeats = 0;
		timer.reset();
		while (timer.getTime() / 1000.0f < 0.2f) {
			repeats++;
			timer.start();
			computeGroundTruth(sampledData, testData, gtMatches, 0, metric);
			timer.stop();
		}

		CostData linearCost = new CostData();
		linearCost.searchTimeCost = (float) timer.getTime() / repeats;
		linearCost.buildTimeCost = 0;
		linearCost.memoryCost = 0;
		costs.add(linearCost);

		// Start parameter autotune process.
		System.out.println("Autotuning parameters...");

		optimizeKMeans(costs);
		optimizeKDTree(costs);

		float bestTimeCost = costs.get(0).buildTimeCost * buildWeight
				+ costs.get(0).searchTimeCost;
		for (int i = 0; i < costs.size(); i++) {
			float timeCost = costs.get(i).buildTimeCost * buildWeight
					+ costs.get(i).searchTimeCost;
			System.out.printf("Time cost: %f", timeCost);
			if (timeCost < bestTimeCost) {
				bestTimeCost = timeCost;
			}
		}
		System.out.printf("Best time cost: %f", bestTimeCost);

		bestBuildParams = costs.get(0).params;
		if (bestTimeCost > 0) {
			float bestCost = costs.get(0).buildTimeCost * buildWeight
					+ costs.get(0).searchTimeCost / bestTimeCost;
			for (int i = 0; i < costs.size(); i++) {
				float crtCost = costs.get(i).buildTimeCost * buildWeight
						+ costs.get(i).searchTimeCost / bestTimeCost
						+ memoryWeight * costs.get(i).memoryCost;
				System.out.printf("Cost: %f\n", crtCost);
				if (crtCost < bestCost) {
					bestCost = crtCost;
					bestBuildParams = costs.get(i).params;
				}
			}
			System.out.printf("Best cost: %f\n", bestCost);
		}
	}

	private float estimateSearchParams(SearchParamsBase searchParams) {
		int nn = 1;
		int SAMPLE_COUNT = 1000;
		
		assert(bestIndex != null);
		
		float speedup = 0;
		
		int samples = Math.min(data.length / 10, SAMPLE_COUNT);
		if(samples > 0) {
			double[][] testDataset = randomSample(data, samples);
			
			System.out.println("Computing ground truth.\n");
			int[][] gtMatches = new int[testDataset.length][1];
			StopWatch timer = new StopWatch();
			int repeats = 0;
			timer.reset();
			while(timer.getTime() / 1000.0f < 0.2) {
				repeats++;
				timer.start();
				computeGroundTruth(data, testDataset, gtMatches, 1, metric);
				timer.stop();
			}
			float linear = timer.getTime() / 1000.0f / repeats;
			int checks;
			System.out.println("Estimating number of checks.");
			float searchTime;
			float cbIndex;
			if (bestIndex_->getType() == FLANN_INDEX_KMEANS) {
				System.out.println("KMeans algorithm, estimating cluster border factor");
				IndexKMeans kmeans = (IndexKMeans)bestIndex;
				float bestSearchTime = -1;
                float best_cb_index = -1;
                int best_checks = -1;
                for (cbIndex = 0; cbIndex < 1.1f; cbIndex += 0.2f) {
                	//kmeans->set_cb_index(cb_index);
                	//                    searchTime = test_index_precision(*kmeans, dataset_, testDataset, gt_matches, target_precision_, checks, distance_, nn, 1);

                	if(searchTime < bestSearchTime || bestSearchTIme == -1) {
                		bestSearchTime = searchTime;
                		best_cb_index = cb_index;
                		best_checks = checks;
                	}
                }
  
                searchTime = bestSearchTime;
                cb_index = best_cb_index;
                checks = best_checks;

                kmeans->set_cb_index(best_cb_index);
                Logger::info("Optimum cb_index: %g\n", cb_index);
                bestParams_["cb_index"] = cb_index;
			} else {
                //searchTime = test_index_precision(*bestIndex_, dataset_, testDataset, gt_matches, target_precision_, checks, distance_, nn, 1);
			}
			
			System.out.printf("Required number of checks: %d \n", checks);
			searchParams.checks = checks;
			speedup = linear / searchTime;
		}
		
		return speedup;
	}

	@Override
	void buildIndex()
    {
        bestParams_ = estimateBuildParams();

        Logger::info("----------------------------------------------------\n");
        Logger::info("Autotuned parameters:\n");
        if (Logger::getLevel()>=FLANN_LOG_INFO)
        	print_params(bestParams_);
        Logger::info("----------------------------------------------------\n");

        flann_algorithm_t index_type = get_param<flann_algorithm_t>(bestParams_,"algorithm");

        bestIndex_ = create_index_by_type(index_type, dataset_, bestParams_, distance_);

        bestIndex_->buildIndex();

        speedup_ = estimateSearchParams(bestSearchParams_);

        Logger::info("----------------------------------------------------\n");
        Logger::info("Search parameters:\n");
        if (Logger::getLevel()>=FLANN_LOG_INFO)
        	print_params(bestSearchParams_);
        Logger::info("----------------------------------------------------\n");

        bestParams_["search_params"] = bestSearchParams_;
        bestParams_["speedup"] = speedup_;
    }

	@Override
	protected void buildIndexImpl() {
		// TODO Auto-generated method stub

	}

	@Override
	protected void findNeighbors(ResultSet resultSet, double[] query,
			SearchParamsBase searchParams) {
		// TODO Auto-generated method stub

	}

	@Override
	protected void findNeighbors(ResultSet resultSet, int[] query,
			SearchParamsBase searchParams) {
		// TODO Auto-generated method stub

	}

}