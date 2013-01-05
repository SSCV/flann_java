package flann.util;



public class Utils {
	public static int genRandomNumberInRange (int min, int max) {
		return min + (int)(Math.random() * ((max - min) + 1));
	}
}