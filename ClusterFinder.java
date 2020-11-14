import org.apache.commons.math3.ml.clustering.Cluster;
import org.apache.commons.math3.ml.clustering.DBSCANClusterer;
import java.util.ArrayList;
import java.util.List;

/*
This class uses a library called "apache clustering", that can be found in the adress:
https://commons.apache.org/proper/commons-math/javadocs/api-3.6/org/apache/commons/math3/ml/clustering/package-summary.html#package_description
"apache clusterer" is a library containing different clustering algorithms.
This class contains a DBSCANClusterer, an object that is used for finding clusters in a list of points.
DBSCANClusterer finds clusters through usage of the DBSCAN algorithm, an algorithm that doesn't require knowing ahead of time the amount of clusters.
*/
public class ClusterFinder {
    // DBSCAN parameters
    private static int EPS = 150;
    private static int MIN_PTS = 25;

    public DBSCANClusterer<Protein> clusterer = new DBSCANClusterer<>(EPS, MIN_PTS);

    // This function gets an array list of proteins, runs a DBSCAN analysis on it, and return the results as an ArrayList of ProteinClusters (see ProteinCluster)
    public ArrayList<ProteinCluster> findClusters(ArrayList<Protein> proteins)
    {
        // Runs the DBSCAN analysis
        List<Cluster<Protein>> clusters = clusterer.cluster(proteins);
        // creates the ArrayList and transfers the information from the results given by clusterer to the ArrayList
        ArrayList<ProteinCluster> proteinClusters = new ArrayList<>();
        for(int i = 0; i < clusters.size(); i++)
        {
            ProteinCluster proteinCluster = new ProteinCluster();
            proteinCluster.cluster = clusters.get(i);
            proteinClusters.add(proteinCluster);
        }
        return proteinClusters;
    }
}
