package flann.index;

import java.util.ArrayList;

public class Bucket {
	ArrayList<Integer> points = new ArrayList<Integer>();

	public void add(int pointIndex) {
		points.add(pointIndex);
	}
}
