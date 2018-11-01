package units.geometria.matrici;

import java.util.Arrays;

/******************************************************************************
 *  Compilation:  javac Matrix.java
 *  Execution:    java Matrix
 *
 *  A bare-bones collection of static methods for manipulating
 *  matrices.
 *
 ******************************************************************************/

public class Matrix {

	/**
	 * Return a random m-by-n matrix with values between 0 and 1.
	 * 
	 * @param m number of rows
	 * @param n number of columns
	 * @return random generated matrix
	 */
    public static double[][] random(int m, int n) {
        double[][] a = new double[m][n];
        for (int i = 0; i < m; i++)
            for (int j = 0; j < n; j++)
                a[i][j] = Math.random();
        return a;
    }

    /**
     * Return n-by-n identity matrix I.
     * 
     * @param n number of rows and columns
     * @return identity matrix I
     */
    public static double[][] identity(int n) {
        double[][] a = new double[n][n];
        for (int i = 0; i < n; i++)
            a[i][i] = 1;
        return a;
    }
    
    /**
     * Return x^T y.
     * 
     * @param x
     * @param y
     * @return
     */
    public static double dot(double[] x, double[] y) {
        if (x.length != y.length) throw new RuntimeException("Illegal vector dimensions.");
        double sum = 0.0;
        for (int i = 0; i < x.length; i++)
            sum += x[i] * y[i];
        return sum;
    }

    /**
     * Let A matrix m-by-n, return B as the trasposed of A.
     * 
     * @param a matrix
     * @return a^T trasposed of a
     */
    public static double[][] transpose(double[][] a) {
        int m = a.length;
        int n = a[0].length;
        double[][] b = new double[n][m];
        for (int i = 0; i < m; i++)
            for (int j = 0; j < n; j++)
                b[j][i] = a[i][j];
        return b;
    }

    /**
     * Let A, B matrix m-by-n, return a matrix obtained with the sum element by element.
     * 
     * @param a first matrix
     * @param b second matrix
     * @return sum
     */
    public static double[][] add(double[][] a, double[][] b) {
        int m = a.length;
        int n = a[0].length;
        double[][] c = new double[m][n];
        for (int i = 0; i < m; i++)
            for (int j = 0; j < n; j++)
                c[i][j] = a[i][j] + b[i][j];
        return c;
    }

    /**
     * @see above method "add"
     * 
     * @param a
     * @param b
     * @return
     */
    public static double[][] subtract(double[][] a, double[][] b) {
        int m = a.length;
        int n = a[0].length;
        double[][] c = new double[m][n];
        for (int i = 0; i < m; i++)
            for (int j = 0; j < n; j++)
                c[i][j] = a[i][j] - b[i][j];
        return c;
    }

    /**
     * Let A matrix m-by-n, B matrix n-by-p, return C matrix obtained from
     * product rows-by-columns. Throws an exception if The number of rows of the second
     * matrix are not equal to the number of Columns of the first one.
     * 
     * @param a first matrix
     * @param b second matrix
     * @return
     * @throws RuntimeException
     */
    public static double[][] multiply(double[][] a, double[][] b) {
        int m = a.length;
        int n1 = a[0].length;
        int n2 = b.length;
        int p = b[0].length;
        if (n1 != n2) throw new RuntimeException("Illegal matrix dimensions.");
        double[][] c = new double[m][p];
        for (int i = 0; i < m; i++)
            for (int j = 0; j < p; j++)
                for (int k = 0; k < n1; k++)
                    c[i][j] += a[i][k] * b[k][j];
        return c;
    }

    /**
     * Matrix-vector multiplication (y = A * x)
     * 
     * @param a matrix
     * @param x vector
     * @return
     * @throws RuntimeException
     */
    public static double[] multiply(double[][] a, double[] x) {
        int m = a.length;
        int n = a[0].length;
        if (x.length != n) throw new RuntimeException("Illegal matrix dimensions.");
        double[] y = new double[m];
        for (int i = 0; i < m; i++)
            for (int j = 0; j < n; j++)
                y[i] += a[i][j] * x[j];
        return y;
    }

    /**
     * Vector-matrix multiplication (y = x^T A)
     * @see above method "multiply"
     * 
     * @param x
     * @param a
     * @return
     * @throws RuntimeException
     */
    public static double[] multiply(double[] x, double[][] a) {
        int m = a.length;
        int n = a[0].length;
        if (x.length != m) throw new RuntimeException("Illegal matrix dimensions.");
        double[] y = new double[n];
        for (int j = 0; j < n; j++)
            for (int i = 0; i < m; i++)
                y[j] += a[i][j] * x[i];
        return y;
    }
    
    /**
     * Scalar-matrix multiplication (y = a * A)
     * 
     * @param a matrix
     * @param x scalar
     * @return
     */
    public static double[][] multiply(double[][] a, double x) {
    	int m = a.length;
        int n = a[0].length;
        double[][] y = new double[m][n];
        for (int i = 0; i < m; i++)
            for (int j = 0; j < n; j++)
                y[i][j] = a[i][j] * x;
        return y;
    }
    
    /**
     * Let A matrix m-by-n, check if is symmetrical (true) of not (false).
     * 
     * @author Lorenzo Elia
     * @version 1.0
     * 
     * @param a matrix
     * @return true if the given matrix is symmetrical
     */
    public static boolean isSymmetrical(double[][] a) {
    	double[][] ta = Matrix.transpose(a);
    	if(Arrays.deepEquals(a, ta))
    		return true;
		return false;
    }
    
    /**
     * Let A matrix m-by-n, check if is diagonal (true) of not (false).
     * 
     * @author Lorenzo Elia
     * @version 1.0
     * 
     * @param a matrix
     * @return true if the given matrix is diagonal
     * @throws RuntimeException
     */
    public static boolean isDiagonal(double[][] a) {
    	int m = a.length;
    	int n = a[0].length;
    	if (m != n) throw new RuntimeException("Not square matrix");
    	for (int i = 0; i < m; i++)
    		for (int j = 0; j < n; j++)
    			if(i != j)
    				if(a[i][j] != 0)
    					return false;
		return true;
    }
    
	/**
	 * Let A a square matrix, let 1 <= i <= n and 1 <= j <= n,
	 * return the minor (i,j) of it, that is the
	 * given matrix without the i-row and the j-column.
	 * 
	 * @author Lorenzo Elia
	 * @version 1.0
	 * @since 2018-11-1
	 * 
	 * @param a square matrix
	 * @param i index on row
	 * @param j index on column
	 * @return the minor (i,j) matrix of a
	 */
	public static double[][] minor(double[][] a, int i, int j) {
		double[][] b = new double[a[0].length - 1][a.length - 1];
		for(int c = 1, e = 1; c <= a[0].length; c++, e++) {
			if(c != i) {
				for(int d = 1, f = 1; d <= a.length; d++, f++) {
					if(d != j) {
						b[e - 1][f - 1] = a[c - 1][d - 1];
					} else {
						f--;
					}
				}
			} else {
				e--;
			}
		}
		return b;
	}

	/**
	 * Let A a square matrix, return the determinant of the given matrix in a recursive manner;
	 * Throws and exception if the square matrix is 0x0.
	 * 
	 * @author Lorenzo Elia
	 * @version 1.0
	 * @since 2018-11-1
	 * 
	 * @param a square matrix
	 * @return determinant of a
	 * @throws Exception
	 */
	public static double det(double[][] a) throws Exception {
		double det = 0;
		if(a == null)
			throw new Exception("Null matrix given");
		if(a.length == 1)
			return a[0][0];
		if(a.length > 1) {
			for(int i = 1; i <= a.length; i++) {
				det += Math.pow(-1, i + 1) * a[i - 1][0] * det(minor(a, i, 1));
			}
		}
		return det;
	}
	
	/**
	 * Let A square matrix, A is reversable if and only if det(A) != 0.
	 * 
	 * @author Lorenzo Elia
	 * @version 1.0
	 * @since 2018-1-11
	 * 
	 * @param a matrix n-by-n
	 * @return true if is reversable
	 * @throws Exception 
	 */
	public static boolean isReversable(double[][] a) throws Exception {
		if(Matrix.det(a) != 0)
			return true;
		return false;
	}

}
