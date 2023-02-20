import java.util.Random;
import java.util.concurrent.Executors;
import java.util.concurrent.ThreadPoolExecutor;
import java.io.PrintWriter;
import java.util.Arrays;
import java.io.File; 
import java.io.FileNotFoundException;
import java.io.IOException; 

public class paralel_2 {
    private static double[] a;
    private static double[][] ma;
    private static int MATRIX_SIZE = 1;

    volatile static double[][] md = genNewMatrix(MATRIX_SIZE);
    volatile static double[][] mt = genNewMatrix(MATRIX_SIZE);
    volatile static double[][] mz = genNewMatrix(MATRIX_SIZE);
    volatile static double[][] me = genNewMatrix(MATRIX_SIZE);
    volatile static double[][] mm = genNewMatrix(MATRIX_SIZE);
   
    volatile private static double[] d = genNewVector(MATRIX_SIZE);
    volatile private static double[] b = genNewVector(MATRIX_SIZE);

	private static ThreadPoolExecutor executorA;
	private static ThreadPoolExecutor executorB;
	static {
		executorA = (ThreadPoolExecutor) Executors.newFixedThreadPool(8);
		executorB = (ThreadPoolExecutor) Executors.newFixedThreadPool(8);
	}

    private static long tA;
    private static long tB;

    public static void main(String[] args) throws FileNotFoundException {
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

        mkFile("2_calcA.txt");
        mkFile("2_calcB.txt");

        try (PrintWriter out = new PrintWriter("2_calcA.txt")) {
            out.println(matrixToString(ma));
        }
        try (PrintWriter out = new PrintWriter("2_calcB.txt")) {
            out.println(Arrays.toString(a));
        }        

        
    }

    public static void calcA() { 			
        // MА = MD*MT+MZ-ME*MM
        
        long startTime = System.currentTimeMillis();
        mulMatrixParallel(executorA, me, mm, me);
        mulMatrixParallel(executorA, md, mt, md);
        
        w8Executopr(executorA);

        mz = subtract(mz, me);
        ma = sum(md, mz);

        
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
		for(int i = 0; i < matA.length; i++) {
			for(int j = 0; j < matB[0].length; j++) {
				executor.submit(new JobMul(matA, matB, mat, i, j));
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

		public JobMul(double[][] matA, double[][] matB, double[][] mat, int i, int j) {
			this.matA = matA;
			this.matB = matB;
			this.mat = mat;
			this.i = i;
			this.j = j;
		}

		@Override
		public void run() {
			double sum = 0;
			for(int k = 0; k < matB.length; k++) {
				sum = sumByKahanAlgo(sum, matA[i][k] * matB[k][j]);
			}
			mat[i][j] = sum;
		}
	}

}


