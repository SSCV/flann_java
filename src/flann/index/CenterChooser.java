package flann.index;

import java.util.ArrayList;

import flann.metric.Metric;
import flann.util.UniqueRandom;
import flann.util.Utils;

public class CenterChooser {
	public static void Random(Metric metric, double[][] data, int k,
			ArrayList<Integer> objectsIndices, int start, int count,
			ArrayList<Integer> centers) {
		UniqueRandom r = new UniqueRandom(count);

		int index;
		for (index = 0; index < k; index++) {
			boolean duplicate = true;
			int rnd;
			while (duplicate) {
				duplicate = false;
				rnd = r.next();
				if (rnd < 0)
					return;

				centers.set(index, objectsIndices.get(start + rnd));

				for (int j = 0; j < index; j++) {
					double sq = metric.distance(data[centers.get(index)],
							data[centers.get(j)]);
					if (sq < 1E-16) {
						duplicate = true;
					}
				}
			}
		}
	}

	public static void Gonzales(Metric metric, double[][] data, int k,
			ArrayList<Integer> objectsIndices, int start, int count,
			ArrayList<Integer> centers) {
		// Pick the first center randomly.
		int n = count;
		int rnd = Utils.genRandomNumberInRange(0, n - 1);
		centers.set(0, objectsIndices.get(start + rnd));

		int index;
		for (index = 1; index < k; index++) {
			int bestIndex = -1;
			double bestValue = 0;
			for (int j = 0; j < n; j++) {
				double dist = metric.distance(data[centers.get(0)],
						data[objectsIndices.get(j)]);
				for (int i = 1; i < index; i++) {
					double tmpDist = metric.distance(data[centers.get(i)],
							data[objectsIndices.get(j)]);
					if (tmpDist < dist) {
						dist = tmpDist;
					}
				}
				if (dist > bestValue) {
					bestValue = dist;
					bestIndex = j;
				}
			}
			if (bestIndex != -1) {
				centers.set(index, objectsIndices.get(bestIndex));
			} else {
				break;
			}
		}
	}

	public static void KMeansPP(Metric metric, double[][] data, int k,
			ArrayList<Integer> objectsIndices, int start, int count,
			ArrayList<Integer> centers) {
		int n = count;
		double currentPotential = 0;
		double[] closestDistSq = new double[n];

		// Chose one random center and set the closestDistSq values.
		int index = Utils.genRandomNumberInRange(0, n - 1);
		centers.set(0, objectsIndices.get(start + index));

		for (int i = 0; i < n; i++) {
			closestDistSq[i] = metric.distance(data[objectsIndices.get(i)],
					data[objectsIndices.get(index)]);
			currentPotential += closestDistSq[i];
		}

		int NUM_LOCAL_TRIES = 1;

		// Choose each center.
		int centerCount;
		for (centerCount = 1; centerCount < k; centerCount++) {
			double bestNewPotential = -1;
			int bestNewIndex = 0;

			for (int localTrial = 0; localTrial < NUM_LOCAL_TRIES; localTrial++) {
				double randVal = Utils.genRandomNumberInRange(0.0,
						currentPotential);
				for (index = 0; index < n - 1; index++) {
					if (randVal <= closestDistSq[index])
						break;
					else
						randVal -= closestDistSq[index];
				}

				// Compute the new potential.
				double newPotential = 0;
				for (int i = 0; i < n; i++) {
					double d = metric.distance(data[objectsIndices.get(i)],
							data[objectsIndices.get(index)]);
					newPotential += Math.min(d, closestDistSq[i]);
				}

				// Store the best result.
				if (bestNewPotential < 0 || newPotential < bestNewPotential) {
					bestNewPotential = newPotential;
					bestNewIndex = index;
				}
			}

			centers.set(centerCount, objectsIndices.get(bestNewIndex));
			currentPotential = bestNewPotential;
			for (int i = 0; i < n; i++) {
				double d = metric.distance(data[objectsIndices.get(i)],
						data[objectsIndices.get(bestNewIndex)]);
				closestDistSq[i] = Math.min(d, closestDistSq[i]);
			}
		}
	}
}