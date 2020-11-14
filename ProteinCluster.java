import org.apache.commons.math3.ml.clustering.Cluster;
import org.apache.commons.math3.ml.clustering.DoublePoint;
import java.util.List;

/*
This class contains a Cluster of Proteins, and its center. Cluster is a class produced by a library called "apache clustering".
More explanation about "apache clustering" is in ClusterFinder.
 */
public class ProteinCluster {
    public Cluster<Protein> cluster = new Cluster<>();
    // DoublePoint is a Class createds by "apache clusterer" containing a point whose coordinates are stores in double variables.
    private DoublePoint center;

    // returns the center of a cluster
    public DoublePoint getCenter()
    {
        // Lazy function that finds the center of the cluster by making an average of the coordinates of all proteins in the cluster
        if(center == null)
        {
            final List<Protein> points = cluster.getPoints();
            if (points.isEmpty()) {
                return null;
            }

            final int dimension = points.get(0).getPoint().length;
            final double[] centroid = new double[dimension];
            for (final Protein p : points) {
                final double[] point = p.getPoint();
                for (int i = 0; i < centroid.length; i++) {
                    centroid[i] += point[i];
                }
            }
            for (int i = 0; i < centroid.length; i++) {
                centroid[i] /= points.size();
            }
            center = new DoublePoint(centroid);
        }
        return center;
    }
}
