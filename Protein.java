import org.apache.commons.math3.ml.clustering.Clusterable;
import java.util.Random;

/*
This class contains a protein in a cell. The protein is represented by x and y coordinates.
 */
public class Protein implements Clusterable {
    // The x and y positions of the protein (in nm)
    public double x;
    public double y;

    // Measures the distance (in nm) between this protein and a given protein
    public double distanceFrom(Protein p2)
    {
        return Math.sqrt(Math.pow(this.x-p2.x, 2)+Math.pow(this.y-p2.y, 2));
    }

    @Override
    public double[] getPoint() {
        return new double[]{x, y};
    }
}
