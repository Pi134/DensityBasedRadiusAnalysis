import com.opencsv.CSVReader;
import com.opencsv.CSVWriter;
import org.apache.commons.io.FileUtils;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.List;

/*
This class takes the result files created by Main, and uses them to create one file.
For each cell, the class writes in the file in how many clusters the density based radius of the green proteins is bigger than
the density based radius of the red proteins, in how many clusters it's smaller, and in how many clusters it's equal.
 */
public class MergeResultFiles {
    public static void main(String[] args) throws IOException
    {
        mergeResultFiles();
    }

    public static void mergeResultFiles() throws IOException {
        // Creates a new file and a way to write into it
        File file = new File("C:\\Users\\aviga\\Documents\\Avigail_108T_Ras_grb2\\Surrouning Analysis\\" + "Density Radius Results Summary" + ".csv");
        FileWriter outputfile = new FileWriter(file);
        CSVWriter writer = new CSVWriter(outputfile);

        // Finds the result files created by Main and puts them in "files"
        File dir = new File("C:\\Users\\aviga\\Documents\\\\Avigail_108T_Ras_grb2\\Surrouning Analysis");
        String[] extensions = new String[] {"csv"};
        List<File> files = (List<File>) FileUtils.listFiles(dir, extensions, true);

        // Writes the titles
        writer.writeNext(new String[]{"protein", "red bigger than green clusters", "green bigger than red clusters", "equal sizes clusters"});

        for(File resultsFile : files)
        {
            //Finds only the files with the density radius results
            if (resultsFile.getName().contains("Density Radiuses"))
            {
                CSVReader reader = new CSVReader(new FileReader(resultsFile));
                // Density radius comparison buckets
                int redBigger = 0;
                int greenBigger = 0;
                int bothEqual = 0;
                reader.readNext(); //reads titles
                String[] nextLine = new String[2];
                while ((nextLine = reader.readNext()) != null) {
                    // nextLine[] is an array of values read from the line
                    // The result files contain "-5" in case a density based radius could not be found and the line should be skipped
                    if(nextLine[0] == "-5" || nextLine[1] == "-5")
                        continue;
                    // Adding 1 to the appropriate bucket
                    if(Integer.parseInt(nextLine[0]) > Integer.parseInt(nextLine[1]))
                        redBigger++;
                    else if (Integer.parseInt(nextLine[1]) > Integer.parseInt(nextLine[0]))
                        greenBigger++;
                    else
                        bothEqual++;
                }
                // Writes down the results from the buckets, in addition to the cell name found in the name of the checked file
                String header = resultsFile.getName();
                header.replace("Density Radiuses", "");
                header.replace("EGF1_PAmCgrb2", "EGF1_PAGFPNRAS_PAmCgrb2_");
                writer.writeNext(new String[]{header, String.valueOf(redBigger), String.valueOf(greenBigger), String.valueOf(bothEqual)});
            }
        }
        writer.close();
    }
}