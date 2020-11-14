/*
This class contains the two histograms (for green proteins and red proteins) for a center of a cluster.
the histograms are calculated in Main.
For each distance measurement, each histogram contains the amount of proteins in the cluster that are within that distance from the center of the cluster.
*/
public class DoubleHistogram {
    // the measurement for the histograms are going to be every 5 nm between 0-1500 nm (including both edges)
    private static int BUCKETS = 1500 / 5 + 1;

    // the histograms
    public int[] redHistogram = new int[BUCKETS];
    public int[] greenHistogram = new int[BUCKETS];

    // the radius of the cluster (in that type of proteins) based on the histogram
    private int redDensityRadius;
    private int greenDensityRadius;

    // lazy functions that return the histogram-based radius, or the density-based radius
    public int getRedDensityRadius()
    {
        if(redDensityRadius == 0)
            redDensityRadius = (getDensityRadius(redHistogram)*5);
        return redDensityRadius;
    }

    public int getGreenDensityRadius()
    {
        if(greenDensityRadius == 0)
            greenDensityRadius = (getDensityRadius(greenHistogram)*5);
        return greenDensityRadius;
    }

    // a function that calculated a density based radius based on a given histogram
    // the function looks for a significant drop in the density of the proteins in the cluster the further you go from the middle, checks to see that the amounts of proteins don't significantly change for the next couple of measurements, and then returns the value of the density based radius.
    private int getDensityRadius(int[] histogram)
    {
        int biggestDifference = 0;
        int taperingCounter = 0;
        for(int i = 1; i < histogram.length; i++)
        {
            if((histogram[i] - histogram[i-1]) > biggestDifference)
                biggestDifference = (histogram[i] - histogram[i-1]);

            if((histogram[i] - histogram[i-1]) < biggestDifference * (15.0/100.0))
            {
                taperingCounter++;
                if(taperingCounter >= 5)
                    return (i-5);
            }
            else
            {
                taperingCounter = 0;
            }
        }
        return -1;
    }
}
