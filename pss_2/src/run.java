import java.util.Random;
import java.util.concurrent.Executors;
import java.util.concurrent.ThreadPoolExecutor;
import java.util.concurrent.locks.Lock;
import java.util.concurrent.locks.ReadWriteLock;
import java.util.concurrent.locks.ReentrantLock;
import java.util.concurrent.locks.ReentrantReadWriteLock;
import java.io.PrintWriter;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Map;
import java.io.BufferedReader;
import java.io.File; 
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException; 

public class run {
    private static double[] a;
    private static double[][] ma;
    private static int INPUT_SIZE;

    private static long tA;
    private static long tB;

    static double[][] md;
    static double[][] mt;
    static double[][] mz;
    static double[][] me;
    static double[][] mm;
    
    static double[][] mdRes;
    static double[][] meRes; 

    static double[] d;
    static double[] b;

    static Map<String, double[][]> inputMatrix = new HashMap<String, double[][]>() {{
        put("md", null);
        put("mt", null);
        put("mz", null);
        put("me", null);
        put("mm", null);
    }};

    static Map<String, double[]> inputVectors = new HashMap<String, double[]>() {{
        put("b", null);
        put("d", null);
    }};
 
    private static ThreadPoolExecutor executorA;
	private static ThreadPoolExecutor executorB;
	static {
		executorA = (ThreadPoolExecutor) Executors.newFixedThreadPool(8);
		executorB = (ThreadPoolExecutor) Executors.newFixedThreadPool(8);
	}

    public static void main(String[] args) throws IOException {
        readInput();
        meRes = new double[INPUT_SIZE][INPUT_SIZE];
        mdRes = new double[INPUT_SIZE][INPUT_SIZE];


        Thread tread = new Thread(() -> {
            calcA();
        });

		tread.start();

        calcB();

        try {
            tread.join();
        } catch (Exception e) {
            System.out.println(e);
        }

        System.out.printf("Bougth multiplications: %d ms\n", tA + tB);

        saveToFile("calc_A.txt", matrixToString(ma));
        saveToFile("calc_B.txt", Arrays.toString(a));
        
        executorA.shutdown();
        executorB.shutdown();
        
    }

    public static void readInput() throws IOException {
        Path currentRelativePath = Paths.get("");
        String s = currentRelativePath.toAbsolutePath().toString();
        
        for (String var : inputMatrix.keySet()) {
            double[][] matrix = parseMatrixFromString(readFileToString(s + "/input/" + var));   
            inputMatrix.replace(var, matrix);
        }
        
        for (String var : inputVectors.keySet()) {
            double[] vector = parseVectorFromString(readFileToString(s + "/input/" + var));   
            inputVectors.replace(var, vector);
        }

        mt = inputMatrix.get("mt");
        md = inputMatrix.get("md");
        mz = inputMatrix.get("mz");
        me = inputMatrix.get("me");
        mm = inputMatrix.get("mm");
        
        b = inputVectors.get("b");
        d = inputVectors.get("d");

        INPUT_SIZE = d.length;
    }
    
    public static String readFileToString(String fName) throws IOException {
        BufferedReader reader = new BufferedReader(new FileReader(fName));
        StringBuilder stringBuilder = new StringBuilder();
        String line = null;
        String ls = System.getProperty("line.separator");
        while ((line = reader.readLine()) != null) {
            stringBuilder.append(line);
            stringBuilder.append(ls);
        }  

        stringBuilder.deleteCharAt(stringBuilder.length() - 1);
        reader.close();

        String content = stringBuilder.toString();
        return content;
    }

    public static double[] parseVectorFromString(String str) {
        str = str.replaceAll("[\\[\\]]", " "); 
        String[] stringVector = str.split(",");
        for (int i = 0; i < stringVector.length; i++) {
            stringVector[i] = stringVector[i].trim();
        }
        double[] doubleValues = Arrays.stream(stringVector)
            .mapToDouble(Double::parseDouble)
            .toArray();
        return doubleValues;
    }

    public static double[][] parseMatrixFromString(String str) {
        String ls = System.getProperty("line.separator");
        
        str = str.replaceAll("[\\[\\]]", " ");  
        String[] rows = str.split(ls);

        double[][] res = new double[rows.length - 1][rows[0].length()];
        for (int i = 0; i < rows.length - 1; i++) {
            rows[i] = rows[i].trim();
            String[] r = rows[i].split(",");
            double[] doubleValues = Arrays.stream(r)
                        .mapToDouble(Double::parseDouble)
                        .toArray();
            res[i] = doubleValues;
        }

        return res;
    }

    public static void saveToFile(String fName, String content) throws FileNotFoundException {
        mkFile(fName);
        try (PrintWriter out = new PrintWriter(fName)) {
            out.println(content);
        }
    }

    public static void calcA() { 			
        // MА = MD*MT+MZ-ME*MM
        
        long startTime = System.currentTimeMillis();
        mulMatrixParallel(executorA, me, mm, meRes);
        mulMatrixParallel(executorA, md, mt, mdRes);
        
        w8Executopr(executorA);

        mz = subtract(mz, meRes);
        ma = sum(mdRes, mz);

        
        long endTime = System.currentTimeMillis();
        tA = endTime - startTime;
        System.out.printf("First multiplication: %d ms\n", tA);
	}

    public static void w8Executopr(ThreadPoolExecutor exe) {
        while (exe.getQueue().size() > 0);
    }

    public static void calcB() { 		
        // A=D*MT-max(D)*B
        
        long startTime = System.currentTimeMillis();
        double maxD = Arrays.stream(d).max().getAsDouble();

        double[] tmp = multiplyMatrixByVector(mt, d);
        a = subtractVectors(tmp, multiplyVectorByScalar(maxD, b));

        long endTime = System.currentTimeMillis();
        tB = endTime - startTime;
        System.out.printf("Second multiplication: %d ms\n", tB);
	}

    public static File mkFile(String fName) {
        try {
            File f = new File(fName);
            if (f.createNewFile()) {
              System.out.println("File created: " + f.getName());
              return f;
            } else {
              System.out.println("File already exists.");
            }
          } catch (IOException e) {
            System.out.println("An error occurred.");
            e.printStackTrace();
          }
        return null;
    }

    public static double[] multiplyMatrixByVector(double[][] matrix, double[] vector) {
        int rows = matrix.length;
        int columns = matrix[0].length;

        double[] result = new double[rows];

        for (int row = 0; row < rows; row++) {
            double sum = 0;
            for (int column = 0; column < columns; column++) {
                double newPart = matrix[row][column] * vector[column];
                sum = sumByKahanAlgo(sum, newPart);
            }
            result[row] = sum;
        }
        return result;
    }

    public static double[] subtractVectors(final double[] a, final double[] b) {
        double[] result = new double[a.length];
        for (int i = 0; i < a.length; i++) {
            result[i] = a[i] - b[i];
        }
        return result;
    }

    public static double[] multiplyVectorByScalar(double scalar, double[] vector) {
        double[] ret = new double[vector.length];
        for (int i = 0; i < vector.length; ++i) {
            ret[i] = scalar * vector[i];
        }
        return ret;
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

    public static double sumByKahanAlgo(double...n) {
        double sum = 0.0; 
          
            // Variable to store the error 
            double x = 0.0; 
      
            // Loop to iterate over the array 
            for (int i=0;i<n.length;i++) { 
      
                double y = n[i] - x; 
                double z = sum + y;
                x = (z - sum) - y;
                sum = z; 
            }
        return sum; 
      
    }

    private static void mulMatrixParallel(ThreadPoolExecutor executor, double[][] matA, double[][] matB, double[][] mat) {
        Lock resultMatrixLock = new ReentrantLock();
        
        for(int i = 0; i < matA.length; i++) {
			for(int j = 0; j < matB[0].length; j++) {
				executor.submit(new JobMul(matA, matB, mat, i, j, resultMatrixLock));
			}
		}
	}

    private static double[][] sum(double[][] first, double[][] second) {
        int row = first.length;
        int column = first[0].length;
        double[][] sum = new double[row][column];
    
        for (int r = 0; r < row; r++) {
            for (int c = 0; c < column; c++) {
                sum[r][c] = sumByKahanAlgo(first[r][c], second[r][c]);
            }
        }
        
        return sum;
    }

    private static double[][] subtract(double[][] first, double[][] second) {
        int row = first.length;
        int column = first[0].length;
        double[][] mat = new double[row][column];
    
        for (int r = 0; r < row; r++) {
            for (int c = 0; c < column; c++) {
                mat[r][c] = first[r][c] - second[r][c];
            }
        }

        return mat;
    }

    private static void print2dArray(double[][] matrix) {
		for (int r = 0; r < matrix.length; r++) {
			for (int c = 0; c < matrix[0].length; c++) {
				System.out.print(matrix[r][c] + "\t");
			}
		}
	}

    private static String matrixToString(double[][] matrix) {
        StringBuilder sb = new StringBuilder();
        for (int i = 0; i < matrix.length; i++) {
            sb.append("[ " + Arrays.toString(matrix[i]) + "]\n");
        }
        return sb.toString();
    }

    private static class JobMul implements Runnable {
		private double[][] matA, matB, mat;
		private int i, j;
		private Lock l;

		public JobMul(double[][] matA, double[][] matB, double[][] mat, int i, int j, Lock l) {
			this.matA = matA;
			this.matB = matB;
			this.mat = mat;
			this.i = i;
			this.j = j;
            this.l = l;
		}

		@Override
		public void run() {
			double sum = 0;
			for(int k = 0; k < matB.length; k++) {
				sum = sumByKahanAlgo(sum, matA[i][k] * matB[k][j]);
			}
            l.lock();
            try {
                mat[i][j] = sum;
            } finally {
                l.unlock();
            }
			
		}
	}

}


