import java.io.IOException;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Map;
import java.util.List;
import java.util.Date;
import java.text.SimpleDateFormat;

import java.util.Collection;
import java.util.stream.Collectors;
import java.util.Iterator;
import java.io.File;
import java.io.FileWriter;
import net.imagej.ImgPlus;
import net.imagej.ops.OpService;
import net.imglib2.RandomAccessibleInterval;

import net.imagej.ImageJ;
import net.imagej.DefaultDataset;

import net.imglib2.IterableInterval;
import net.imglib2.algorithm.labeling.ConnectedComponents.StructuringElement;
import net.imglib2.roi.labeling.LabelRegion;
import net.imglib2.roi.labeling.LabelRegions;
import net.imglib2.roi.labeling.ImgLabeling;
import net.imglib2.type.numeric.real.FloatType;
import net.imglib2.type.numeric.real.DoubleType;
import net.imglib2.img.Img;

public class NucFinder {

  private static final String DATA_PATH = "/home/aaron/Documents/idorsia/testData/dapi";
  private static final int SIZE_THRESHOLD = 50;
  public static void main(String[] args) {

    //Create an ImageJ instance.
    final ImageJ ij = new ImageJ();

    //Read the data files from the target directory.
    List<File> files = Arrays.asList(new File(args[0]).listFiles());

    //Stream the files through the function which performs the analysis and collect the results.
    Map<String, Integer> results2 = files.stream()
            .map(file -> count_nuclei(ij, file.getPath()))
            .map(Map::entrySet)
            .flatMap(Collection::stream)
            .collect(Collectors.toMap(Map.Entry::getKey, Map.Entry::getValue, Integer::max));

    //Write the results to a simple csv file.
    simpleMapToCSV(results2);

    //Exit the application.
    System.exit(0);
    }


  /**

   Perform some basic segmentation using and ImageJ and return a count of nuclei for an image file.

   @param ij An ImageJ context for doing analysis operations

   @param filepath A path to the image file which we will attempt to count the nuclei.

   **/
  private static HashMap<String, Integer> count_nuclei(ImageJ ij, String filepath) {

    ij.log().info("Reading file: " + filepath);
    OpService ops = ij.op();

    try {
      // load the dataset
      ImgPlus img = ((DefaultDataset) ij.io().open(filepath)).getImgPlus();

      //Convert the image to float32.  
      Img nuclei = ops.convert().float32(img);
      Img background = (Img) ij.op().filter().gauss(nuclei,100);

      float f = ops.stats().mean((Iterable) nuclei).getRealFloat();
      FloatType mean = new FloatType(f);

      Img bkgCorrection = (Img) ops.math().add(background, new FloatType(f));
      Img corrNuclei = (Img) ops.math().subtract(nuclei, (IterableInterval) bkgCorrection);
      Img smoothCorrNuclei = (Img) ops.filter().gauss(corrNuclei, 2);
      Img threshed = (Img) ops.threshold().otsu(smoothCorrNuclei);


      ImgLabeling labels = ops.labeling().cca(threshed, StructuringElement.FOUR_CONNECTED);
      LabelRegions regions = new LabelRegions(labels);

      int count = 0;
      Iterator rIterator = regions.iterator();

      //Count only if the region is large enough to realistically be a nucleus
      while (rIterator.hasNext()) {
        LabelRegion region = (LabelRegion) rIterator.next();
        DoubleType size = ops.geom().size(region);
        if (size.getRealDouble() >= SIZE_THRESHOLD) {count++;}
      }

      HashMap <String, Integer> counts = new HashMap<String, Integer>();
      String fileName = new File(filepath).getName();
      counts.put(fileName, count);
      return counts;
    } catch (IOException e) {
      e.printStackTrace();
    }
    return new HashMap<String, Integer>();
  }

  private static void simpleMapToCSV(Map<String,Integer> map){
    SimpleDateFormat dateFormat = new SimpleDateFormat("HH_mm_ss");
    Date date = new Date();
    String ts = dateFormat.format(date);
    map.entrySet().forEach(entry -> System.out.println(entry));

    String lsep = System.getProperty("line.separator");
    String psep = System.getProperty("path.separator");

    try (FileWriter writer = new FileWriter("results" + psep + ts + ".csv")) {
      for (Map.Entry<String, Integer> entry : map.entrySet()) {
        writer.append(entry.getKey())
                .append(',')
                .append(entry.getValue().toString())
                .append(lsep);
      }
    } catch (IOException e) {
      e.printStackTrace();
    }
  }
}
//EOF