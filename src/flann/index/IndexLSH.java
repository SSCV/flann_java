package flann.index;

import java.util.ArrayList;
import java.util.BitSet;
import java.util.HashMap;
import java.util.Map.Entry;
import java.util.Set;

import flann.metric.Metric;
import flann.result_set.ResultSet;

public class IndexLSH extends IndexBase {
	// Hash tables, one optimized for speed, the other for space.
	BucketsSpeed bucketsSpeed;
	BucketsSpace bucketsSpace;

	// Decides what will be used to store/manage the data.
	SpeedLevel speedLevel;

	// A speedup so that we don't look in the hash table,
	// which can be much slower than checking a bitset.
	BitSet keyBitSet;

	// Size of the key in bits, used in the hash tables.
	int keySize;

	// The mask to apply to a point to get the hask key.
	ArrayList<Integer> mask;

	/**
	 * Objects of this type are stored in an LSH bucket.
	 */
	private class PointIndex {
		public int val;

		public PointIndex(int val) {
			this.val = val;
		}
	}

	/**
	 * The id from which we can get a bucket in an LSH table.
	 */
	private class BucketKey {
		public int val;

		public BucketKey(int val) {
			this.val = val;
		}
	}

	/**
	 * A bucket in an LSH table.
	 */
	private class Bucket {
		public ArrayList<PointIndex> points = new ArrayList<PointIndex>();
	}

	/**
	 * A hash table which contains buckets, and each bucket is indexed by a
	 * bucket key. Optimized for space.
	 */
	private class BucketsSpace {
		HashMap<BucketKey, Bucket> buckets = new HashMap<BucketKey, Bucket>();
	}

	/**
	 * A hash table which contains buckets, and each bucket is indexed by a
	 * bucket key. Optimized for speed.
	 */
	private class BucketsSpeed {
		ArrayList<Bucket> buckets = new ArrayList<Bucket>();
	}

	/**
	 * Defines the speed for the implementation. kArray uses an ArrayList for
	 * storing data. kBitsetHash uses a hash map but checks for the validity of
	 * a key with a bitset. kHash uses a hash map only.
	 */
	enum SpeedLevel {
		kArray, kBitsetHash, kHash
	};

	public void add(int pointIndex, int[] point) {
		BucketKey key = new BucketKey(1);
		int k = key.val;
		PointIndex p = new PointIndex(k);

		switch (speedLevel) {
		case kArray:
			bucketsSpeed.buckets.get(k).points.add(p);
			break;
		case kBitsetHash:
			keyBitSet.set(k);
			bucketsSpace.buckets.get(k).points.add(p);
			break;
		case kHash:
			bucketsSpace.buckets.get(k).points.add(p);
			break;
		}
	}

	private class PointIndexPointPair {
		int pointIndex;
		int[] point;
	}

	public void add(ArrayList<PointIndexPointPair> points) {
		for (int i = 0; i < points.size(); i++) {
			add(points.get(i).pointIndex, points.get(i).point);
		}
		// optimize();
	}

	private Bucket getBucketFromKey(BucketKey key) {
		int k = key.val;

		switch (speedLevel) {
		case kArray:
			return bucketsSpeed.buckets.get(k);
		case kBitsetHash:
			if (keyBitSet.get(k))
				return bucketsSpace.buckets.get(k);
			else
				return null;
		case kHash:
			return bucketsSpace.buckets.get(k);
		}

		return null;
	}

	public int getKey(int[] point) {
		return 1;
	}

	public IndexLSH(Metric metric, double[][] data) {
		super(metric, data);

		// Initialize bucketsSpeed.
		int bucketsSpeedSize = (1 << keySize);
		for (int i = 0; i < bucketsSpeedSize; i++) {
			bucketsSpeed.buckets.add(new Bucket());
		}
	}

	private void optimize() {
		// If we are already using fast storage, don't do anything.
		if (speedLevel == SpeedLevel.kArray)
			return;

		if (bucketsSpace.buckets.size() > (1 << keySize) / 2) {
			speedLevel = SpeedLevel.kArray;
			int size = 1 << keySize;
			bucketsSpeed.buckets.clear();
			for (int i = 0; i < size; i++) {
				bucketsSpeed.buckets.add(null);
			}
			Set<Entry<BucketKey, Bucket>> s = bucketsSpace.buckets.entrySet();
			for (Entry<BucketKey, Bucket> e : s) {
				int key = e.getKey().val;
				Bucket b = e.getValue();
				bucketsSpeed.buckets.set(key, b);
			}
			bucketsSpace.buckets.clear();
			return;
		}

		if (false) {
			// TODO: Not sure about this...
		} else {
			speedLevel = SpeedLevel.kHash;
			keyBitSet.clear();
		}
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

}
