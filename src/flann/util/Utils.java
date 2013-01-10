package flann.util;

public class Utils {
	// Generate a random number min <= x <= max.
	public static int genRandomNumberInRange(int min, int max) {
		return min + (int) (Math.random() * ((max - min) + 1));
	}

	// Generate a random number min <= x < max.
	public static double genRandomNumberInRange(double min, double max) {
		return min + Math.random() * (max - min);
	}

	public static boolean swapArray(int[] array, int i, int j) {
		int size = array.length;
		if (size < 2 || i == j || i < 0 || i >= size || j < 0 || j >= size) {
			return false;
		}
		int temp = array[i];
		array[i] = array[j];
		array[j] = temp;
		return true;
	}
}