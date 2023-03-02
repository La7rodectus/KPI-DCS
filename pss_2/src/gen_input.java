import java.util.Random;
import java.io.PrintWriter;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.io.File;  // Import the File class
import java.io.FileNotFoundException;
import java.io.IOException;  // Import the IOException class to handle errors

public class gen_input {

    private static int INPUT_SIZE = 1000;

    public static void main(String[] args) throws FileNotFoundException {
        genInputFolderWithFiles(new String[] { "md", "mt", "mz", "me", "mm", "b", "d" });
    }

    public static List<File> genInputFolderWithFiles(String[] fNames) throws FileNotFoundException {
        Path currentRelativePath = Paths.get("");
        String s = currentRelativePath.toAbsolutePath().toString();
        
        mkDir(s + "/input");

        File f = null;
        double[][] matrix = null;
        double[] vector = null;
        
        List<File> res = new ArrayList<>();

        for (String fn : fNames) {
            f = mkFile(s + "/input/" + fn);

            if (fn.length() == 1) {
                vector = genNewVector(INPUT_SIZE);
                saveToFile(f, Arrays.toString(vector));
            } else {
                matrix = genNewMatrix(INPUT_SIZE);
                saveToFile(f, matrixToString(matrix));
            }
            
            res.add(f);
        }

        return res;
    }

    public static void saveToFile(File f, String content) throws FileNotFoundException {
        try (PrintWriter out = new PrintWriter(f.getAbsolutePath())) {
            out.println(content);
        }
    }

    public static void mkDir(String path) {
        File theDir = new File(path);
        if (!theDir.exists()) theDir.mkdirs();
    }

    public static File mkFile(String fName) {
        try {
            File f = new File(fName);
            if (f.createNewFile()) {
              System.out.println("File created: " + f.getName());
              return f;
            } else {
              System.out.println("File already exists.");
              return f;
            }
          } catch (IOException e) {
            System.out.println("An error occurred.");
            e.printStackTrace();
          }
        return null;
    }


	public static double[][] genNewMatrix(int n) { 		
		double [][]A=new double[n][n]; 	
		Random R=new Random(); 	
		
		
		int i,j; 	
		for( i=0; i < n ; i++ ) {
			for( j=0; j < n ; j++ ) {
				A[i][j]=R.nextInt(256); 	
			}	
		}	
		return A; 	
		
	}

    public static double[] genNewVector(int len) { 		
		double []A=new double[len]; 	
		Random R=new Random(); 	
		
		int i; 	
		for( i=0; i < len ; i++ ) {
		    A[i]=R.nextInt(256); 	
		
		}	
		return A; 	
		
	}

    private static String matrixToString(double[][] matrix) {
        StringBuilder sb = new StringBuilder();
        for (int i = 0; i < matrix.length; i++) {
            sb.append("[ " + Arrays.toString(matrix[i]) + "]\n");
        }
        return sb.toString();
    }
}


