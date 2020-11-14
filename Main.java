import com.opencsv.CSVReader;
import com.opencsv.CSVWriter;
import org.apache.commons.io.FileUtils;
import org.apache.commons.math3.ml.clustering.*;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/*
This class gets files that specify the coordinates of the proteins in the cells from the experiment.
Each cell's information is divided into 2 files: one for Grb2 (red) proteins, and one for Ras (green) proteins.
The class creates a file for each cell, containing 2 histograms for each cell - one for red proteins and one for green proteins.
The histogram specifies for each distance between - and 1500 nm (in jumps o 5 nm), how many proteins in the cluster are within that distance from the center of the cluster.
In addition, this class creates a file for each cell that contains for each cluster the density based radii of the green and red proteins in the cluster, based on the previously calculated histograms.
 */
public class Main {
    public static void main(String[] args) throws IOException {
        boolean debug = true;

        // Finds the protein files
        List<File> files = findFilesInDirectory();
        // This loop runs on every Grb2 (red) file, finds its Ras (green) counterpart, and uses both
        for(File file : files) {
            if (file.getName().contains("PAmCgrb2")) {

                if (debug)
                    System.out.println(file.getName());

                File redData = file;
                ArrayList<Protein> redProteins = createArray(redData);

                String greenFileName = file.getName().replace("PAmCgrb2", "PAGFPNRAS");
                File greenData = new File("C:\\Users\\aviga\\Documents\\\u05D0\u05DC\u05E4\u05D0\\Avigail_108T_Ras_grb2\\\u05D8\u05D1\u05DC\u05D0\u05D5\u05EA \u05E1\u05D5\u05E4\u05E8 \u05E8\u05D6\u05D5\u05DC\u05D5\u05E6\u05D9\u05D9\u05D4\\" + greenFileName);
                ArrayList<Protein> greenProteins = createArray(greenData);

                ClusterFinder clusterFinder = new ClusterFinder();

                // Uses clusterFinder to find all clusters in the arrays
                ArrayList<ProteinCluster> redClusters = clusterFinder.findClusters(redProteins);
                if (debug)
                    System.out.println(redClusters.size());

                ArrayList<ProteinCluster> greenClusters = clusterFinder.findClusters(greenProteins);
                if (debug)
                    System.out.println(greenClusters.size());

                // Creates clusters of proteins containing all proteins from all clusters
                Cluster<Protein> redAllClusters = new Cluster<>();
                for (ProteinCluster redCluster : redClusters) {
                    for (Protein redProtein : redCluster.cluster.getPoints())
                        redAllClusters.addPoint(redProtein);
                }

                Cluster<Protein> greenAllClusters = new Cluster<>();
                for (ProteinCluster greenCluster : greenClusters) {
                    for (Protein greenProtein : greenCluster.cluster.getPoints())
                        greenAllClusters.addPoint(greenProtein);
                }

                makeFile(redAllClusters, "red clusters " + redData.getName());
                makeFile(greenAllClusters, "green clusters" + greenData.getName());

                if (debug)
                {
                    System.out.println("*********");
                    System.out.println("overall num of red proteins: " + redAllClusters.getPoints().size());
                    System.out.println("overall num of green proteins: " + greenAllClusters.getPoints().size());
                }

                // Creates ArrayLists containing the centers of all clusters
                // DoublePoint is a class from apache clustering containing a point, whose coordinates are double values
                ArrayList<DoublePoint> redCenters = new ArrayList<>();
                for (ProteinCluster redCluster : redClusters)
                    redCenters.add(redCluster.getCenter());

                ArrayList<DoublePoint> greenCenters = new ArrayList<>();
                for (ProteinCluster greenCluster : greenClusters)
                    greenCenters.add(greenCluster.getCenter());

                HashMap<ProteinCluster, DoubleHistogram> clusterHistogram = new HashMap<>();
                // For every red cluster, creates a new entry with the cluster as the key and in the value histograms that contains how many proteins from each color are within a certain distance of the center of the cluster, for distances between 0 and 1500 nm, measuring every 5 nm.
                for (ProteinCluster redCluster : redClusters) {
                    //creates a new DoubleHistogram, containing 2 histograms - red and green
                    DoubleHistogram histograms = new DoubleHistogram();
                    // Runs on all green proteins in the cluster
                    for (Protein greenProtein : greenProteins) {
                        // Finds the distance between the protein and the center of the cluster
                        DoublePoint gp = new DoublePoint(greenProtein.getPoint());
                        double d = distance(gp, redCluster.getCenter());
                        // Runs on all buckets in the histogram that represent a distance bigger or equal to the distance from the protein to the center of the cluster and adds 1 to them
                        for (int i = (int) Math.ceil((d / 5)); i < 301; i++) {
                            histograms.greenHistogram[i]++;
                        }
                    }
                    // Runs on all red proteins in the cluster and deos the same thing as it did with the green proteins, in a different histogram
                    for (Protein redProtein : redProteins) {
                        DoublePoint rp = new DoublePoint(redProtein.getPoint());
                        double d = distance(rp, redCluster.getCenter());
                        for (int i = (int) Math.ceil((d / 5)); i < 301; i++) {
                            histograms.redHistogram[i]++;
                        }
                    }
                    // Inserts the values calculated into the HashMap
                    clusterHistogram.put(redCluster, histograms);
                }


                // Creates a new file
                File file1 = new File("C:\\Users\\aviga\\Documents\\Avigail_108T_Ras_grb2\\Surrouning Analysis\\" + "Histograms" + redData.getName());
                FileWriter outputfile = new FileWriter(file1);
                CSVWriter writer = new CSVWriter(outputfile);

                // For each bucket in the histogram, for all double histograms, writes in a line the bucket's value in the green histogram and in the red histogram in the double histogram
                for (int j = 0; j < 301; j++) {
                    String[] line = new String[302];
                    line[0] = String.valueOf((j * 5)); // adds a header to the line
                    int i = 1;
                    // Runs on all doubleHistograms in the HashMap and adds to the line the values in that bucket in both histograms (red and green)
                    for (DoubleHistogram doubleHistogram : clusterHistogram.values()) {
                        line[i] = String.valueOf(doubleHistogram.greenHistogram[j]);
                        line[i + 1] = String.valueOf(doubleHistogram.redHistogram[j]);
                        i = i + 2;
                    }
                    // Writes the line
                    writer.writeNext(line);
                }
                writer.close();

                File radius = new File("C:\\Users\\aviga\\Documents\\Avigail_108T_Ras_grb2\\Surrouning Analysis\\" + "Density Radiuses" + redData.getName());
                FileWriter outputRadius = new FileWriter(radius);
                CSVWriter writerRadius = new CSVWriter(outputRadius);

                // Creates a file that contains for each cluster the coordinates of its center as well as the density based radii of the green and red proteins in the cluster
                // For further explanation on density-based radius look at class "DoubleHistogram"
                String[] header = {"red density radius [nm]", "[green density radius [nm]", "x value of center [nm]", "y value of center [nm]"};
                writerRadius.writeNext(header);
                // For each cluster, writes a line containing the density based radii of the two types of proteins in the cluster and the coordinates for the cluster's center, and writes the line into the file
                for (Map.Entry<ProteinCluster, DoubleHistogram> entry : clusterHistogram.entrySet()) {
                    String[] line = new String[4];
                    line[0] = String.valueOf(entry.getValue().getRedDensityRadius());
                    line[1] = String.valueOf(entry.getValue().getGreenDensityRadius());
                    line[2] = String.valueOf(entry.getKey().getCenter().getPoint()[0]);
                    line[3] = String.valueOf(entry.getKey().getCenter().getPoint()[1]);
                    writerRadius.writeNext(line);
                }
                writerRadius.close();


            }
        }
    }

    // This function gets and reads a protein file from the experiment, and returns an ArrayList containing all proteins in the file
    public static ArrayList<Protein> createArray(File file) throws IOException {
        ArrayList<Protein> proteins = new ArrayList<>();

        FileReader filereader = new FileReader(file);
        CSVReader csvReader = new CSVReader(filereader);
        String[] nextRecord;
        csvReader.readNext(); // Read titles

        while ((nextRecord = csvReader.readNext()) != null) {
            Protein protein = new Protein();
            protein.x = Double.parseDouble(nextRecord[1]);
            protein.y = Double.parseDouble(nextRecord[2]);
            proteins.add(protein);
        }
        return proteins;
    }

    // This function gets a cluster of proteins and a name and creates a file under that name which contains a list of all proteins in the cluster and their coordinates
    public static void makeFile(Cluster<Protein> cluster, String name) throws IOException
    {
        List<Protein> proteinsInCluster = cluster.getPoints();

        // Creates a new file using the name parameter and a way to write into it
        File file = new File("C:\\Users\\aviga\\Documents\\Avigail_108T_Ras_grb2\\Surrouning Analysis\\" + name + ".csv");
        FileWriter outputfile = new FileWriter(file);
        CSVWriter writer = new CSVWriter(outputfile);

        // Writes headers
        String[] header = {"x[nm]", "y[nm]"};
        writer.writeNext(header);

        // Runs on all proteins in the list of proteins and writes them into the array
        for(Protein protein : proteinsInCluster)
        {
            String[] proteinInCluster = {Double.toString(protein.getPoint()[0]), Double.toString(protein.getPoint()[1])};
            writer.writeNext(proteinInCluster);
        }
        writer.close();
    }

    // This function gets 2 DoublePoints and returns the distance between them using the pythagorean theorem
    public static double distance(DoublePoint p1, DoublePoint p2)
    {
        return Math.sqrt(Math.pow(p1.getPoint()[0]-p2.getPoint()[0], 2)+Math.pow(p1.getPoint()[1]-p2.getPoint()[1], 2));
    }

    // This function returns all file in the directory where the experiment's result files are stored
    public static List<File> findFilesInDirectory() throws IOException {
        File dir = new File("C:\\Users\\aviga\\Documents\\Avigail_108T_Ras_grb2");
        String[] extensions = new String[] {"csv"};
        List<File> files = (List<File>) FileUtils.listFiles(dir, extensions, true);
        return files;
    }
}