package com.tomfielden.util;

import java.util.*;
import java.io.*;

/*  The default interpretation of the storage is compressed diagonal
    The first index i represents the diagonal offset (0 => THE diagonal)
    The second index j represents the element along the diagonal with
    respect to the left or top edge of the matrix.

    No Zeros Allowed.
 */

public class SparseMatrix {

    static private final String SPACE         = " ";
    static private final String NEWLINE       = "\n";
    static private final String OPEN_BRACKET  = "[";
    static private final String CLOSE_BRACKET = "]";
    
    private int m_rows;
    private int m_cols;
    private HashMap<Integer, HashMap<Integer, Double>> m_matrix = new HashMap<Integer, HashMap<Integer, Double>>();

    // linear access cache. Needed to ensure synchronization and is more efficient.
    private int   [] m_row_indicies = null;
    private int   [] m_col_indicies = null;
    private double[] m_values       = null;  // No Zero's allowed

    public SparseMatrix() {
	m_rows = 0;
	m_cols = 0;
    }

    public SparseMatrix(double[] vector) {
	m_rows = vector.length;
	m_cols = vector.length;
	populateDiagonal(0, vector);
    }

    public SparseMatrix(int rows, int cols) {
	m_rows = rows;
	m_cols = cols;
    }

    public SparseMatrix(int rows, int cols, int sideband) {
	// diagonal units of width of sideband, 0 => identity
	m_rows = rows;
	m_cols = cols;
	if(sideband < 0) for(int i = sideband; i <= 0; i++) populateDiagonal(i);
	else             for(int i = 0; i <= sideband; i++) populateDiagonal(i);
    }

    public SparseMatrix(SparseMatrix A) {
	m_rows = A.m_rows;
	m_cols = A.m_cols;
	for(Integer i: A.m_matrix.keySet()) {
	    HashMap<Integer, Double> diagA = A.m_matrix.get(i);
	    HashMap<Integer, Double> diagB = getDiagonalFailsafe(i);
	    for(Integer j: diagA.keySet()) diagB.put(j, diagA.get(j));
	}
    }

    public SparseMatrix copy() { return new SparseMatrix(this); }

    public void setDiagonal(int i, Double[] V) {
	m_rows = Math.max(m_rows, V.length + Math.max(i,0));
	m_cols = Math.max(m_cols, V.length - Math.min(i,0));

	double[] v = new double[V.length];
	for(int j = 0; j < V.length; j++) v[j] = V[j];
	populateDiagonal(i, v);
    }

    public SparseMatrix transpose() {
	SparseMatrix B = new SparseMatrix(m_cols, m_rows);               // reversed
	for(Integer i: m_matrix.keySet()) {
	    HashMap<Integer, Double> diagA = m_matrix.get(i);
	    HashMap<Integer, Double> diagB = B.getDiagonalFailsafe(-i);  // key step
	    for(Integer j: diagA.keySet()) diagB.put(j, diagA.get(j));
	}
	return B;
    }

    public int size() {
	int s = 0;
	for(Integer i: m_matrix.keySet()) s += m_matrix.get(i).size();
	return s;
    }

    public int rows() { return m_rows; }
    public int cols() { return m_cols; }

    // Linear Access is auto-called and cached structures auto-reset on update.
    private void resetLinearAccess() {
	m_row_indicies = null;
	m_col_indicies = null;
	m_values       = null;
    }
    private void prepareLinearAccess() {
	int N          = size();
	m_row_indicies = new int   [N];
	m_col_indicies = new int   [N];
	m_values       = new double[N];
	int k          = 0;
	for(Integer i: m_matrix.keySet()) {
	    HashMap<Integer, Double> diag = m_matrix.get(i);
	    for(Integer j: diag.keySet()) {
		m_row_indicies[k] = getRow(i,j);
		m_col_indicies[k] = getCol(i,j);
		m_values      [k] = diag.get(j);
		k++;
	    }
	}
    }
    public int[] getRowIndicies() {
	if(m_row_indicies == null) prepareLinearAccess();
	return m_row_indicies;
    }
    public int[] getColIndicies() {
	if(m_col_indicies == null) prepareLinearAccess();
	return m_col_indicies;
    }
    public double[] getValues() {
	if(m_values == null) prepareLinearAccess();
	return m_values;
    }

    public static SparseMatrix add(SparseMatrix A, SparseMatrix B) {
	SparseMatrix C = new SparseMatrix(A.m_rows, A.m_cols);
	Set<Integer> iA  = A.m_matrix.keySet();
	Set<Integer> iB  = B.m_matrix.keySet();
	Set<Integer> iAB = new HashSet<Integer>();
	iAB.addAll(iA);
	iAB.addAll(iB);
	for(Integer i: iAB) {
	    Boolean A_hit = A.m_matrix.containsKey(i);
	    Boolean B_hit = B.m_matrix.containsKey(i);
	    if(A_hit && B_hit) {
		HashMap<Integer, Double> diagA = A.m_matrix.get(i);
		HashMap<Integer, Double> diagB = B.m_matrix.get(i);
		HashMap<Integer, Double> diagC = C.getDiagonalFailsafe(i);
		Set<Integer> jA  = diagA.keySet();
		Set<Integer> jB  = diagB.keySet();
		Set<Integer> jAB = new HashSet<Integer>();
		jAB.addAll(jA);
		jAB.addAll(jB);
		for(Integer j: jAB) {
		    Boolean Ai_hit = diagA.containsKey(j);
		    Boolean Bi_hit = diagB.containsKey(j);
		    if(Ai_hit && Bi_hit) diagC.put(j, diagA.get(j) + diagB.get(j));
		    else if(Ai_hit)      diagC.put(j, diagA.get(j));
		    else if(Bi_hit)      diagC.put(j, diagB.get(j));
		}
	    } else if(A_hit) {
		HashMap<Integer, Double> diagA = A.m_matrix.get(i);
		HashMap<Integer, Double> diagC = C.getDiagonalFailsafe(i);
		for(Integer j: diagA.keySet()) diagC.put(j, diagA.get(j));
	    } else if(B_hit) {
		HashMap<Integer, Double> diagB = B.m_matrix.get(i);
		HashMap<Integer, Double> diagC = C.getDiagonalFailsafe(i);
		for(Integer j: diagB.keySet()) diagC.put(j, diagB.get(j));
	    }
	}
	return C;
    }

    public static SparseMatrix mult(SparseMatrix A, double b) {
	SparseMatrix C = new SparseMatrix(A.m_rows, A.m_cols);
	if(b == 0.0) return C;  // all zeros special case
	for(Integer i: A.m_matrix.keySet()) {
	    HashMap<Integer, Double> diagA = A.m_matrix.get(i);
	    HashMap<Integer, Double> diagC = C.getDiagonalFailsafe(i);
	    for(Integer j: diagA.keySet()) diagC.put(j, diagA.get(j) * b);
	}
	return C;
    }

    public static Double[] mult(SparseMatrix A, Double[] V) {
	// Xr += Aij * Vc where r = getRow(i,j), c = getCol(i,j)
	if(V.length != A.m_cols) {
	    System.out.println("SparseMatrix("+A.m_rows+", "+A.m_cols+")"+
			       " and Vector("+V.length+") not conformal.");
	    return new Double[0];
	}
	Double[] vector = new Double[A.m_rows];
	for(int r = 0; r < A.m_rows; r++) vector[r] = 0.0;
	for(Integer i: A.m_matrix.keySet()) {
	    HashMap<Integer, Double> diag = A.m_matrix.get(i);
	    for(Integer j: diag.keySet()) 
		vector[getRow(i, j)] += diag.get(j) * V[getCol(i,j)];
	}
	return vector;
    }

    public static Double[] mult(Double[] V, SparseMatrix A) {
	// Xc += Vr * Aij where r = getRow(i,j), c = getCol(i,j)
	if(V.length != A.m_rows) {
	    System.out.println("SparseMatrix("+A.m_rows+", "+A.m_cols+")"+
			       " and Vector("+V.length+") not conformal.");
	    return new Double[0];
	}
	Double[] vector = new Double[A.m_cols];
	for(int r = 0; r < A.m_cols; r++) vector[r] = 0.0;
	for(Integer i: A.m_matrix.keySet()) {
	    HashMap<Integer, Double> diag = A.m_matrix.get(i);
	    for(Integer j: diag.keySet()) 
		vector[getCol(i, j)] += V[getRow(i,j)] * diag.get(j);
	}
	return vector;
    }

    public static SparseMatrix mult(SparseMatrix A, SparseMatrix B) {
	// We multiply each A.i diagonal by each B.i diagonal with the result being
	//   stored in the result A.i+B.i diagonal.
	// Within diagonal multiplies we need to line up the elements by translating
	//   A.j indicies into corresponding B.j indicies for lookup.
	SparseMatrix C = new SparseMatrix(A.m_rows, B.m_cols);
	for(Integer Ai: A.m_matrix.keySet()) {
	    HashMap<Integer, Double> Adiag = A.m_matrix.get(Ai);
	    for(Integer Bi: B.m_matrix.keySet()) {
		int Ci = Ai + Bi;
		if(!C.isValidDiagonal(Ci)) continue;
		HashMap<Integer, Double> Bdiag = B.m_matrix.get(Bi);
		for(Integer Aj: Adiag.keySet()) {
		    int Bj = Aj - Ai - Bi;
		    int Cj = Ai + Bj;
		    if(!C.isValidIdx(Ci, Cj)) continue;
		    if(Bdiag.containsKey(Bj)) {
			// System.out.println("A"+Ai+Aj+" * B"+Bi+Bj+" -> C"+Ci+Cj); // DEBUG
			C.accum(Ci, Cj, Adiag.get(Aj) * Bdiag.get(Bj));
		    }
		}
	    }
	}
	return C;
    }

    public static SparseMatrix hadamard(Double[] V, SparseMatrix A) {
	if(V.length != A.m_cols) {
	    System.out.println("SparseMatrix("+A.m_rows+", "+A.m_cols+")"+
			       " and Vector("+V.length+") not conformal.");
	    return new SparseMatrix();
	}
	SparseMatrix C = new SparseMatrix(A.m_rows, A.m_cols);
	for(Integer i: A.m_matrix.keySet()) {
	    HashMap<Integer, Double> diagA = A.m_matrix.get(i);
	    HashMap<Integer, Double> diagC = C.getDiagonalFailsafe(i);
	    for(Integer j: diagA.keySet()) diagC.put(j, diagA.get(j) * V[getCol(i,j)]);
	}
	return C;
    }

    public static SparseMatrix hadamard(SparseMatrix A, Double[] V) {
	if(V.length != A.m_rows) {
	    System.out.println("SparseMatrix("+A.m_rows+", "+A.m_cols+")"+
			       " and Vector("+V.length+") not conformal.");
	    return new SparseMatrix();
	}
	SparseMatrix C = new SparseMatrix(A.m_rows, A.m_cols);
	for(Integer i: A.m_matrix.keySet()) {
	    HashMap<Integer, Double> diagA = A.m_matrix.get(i);
	    HashMap<Integer, Double> diagC = C.getDiagonalFailsafe(i);
	    for(Integer j: diagA.keySet()) diagC.put(j, diagA.get(j) * V[getRow(i,j)]);
	}
	return C;
    }

    public static SparseMatrix pow(SparseMatrix A, int n) {
	SparseMatrix B = A.copy();
	for(int i = 1; i < n; i++) B = mult(B, A);
	return B;
    }

    // NB: A ** <blank> = 0. If you put in a zero then A**0 = 1 else blank. This is painful.
    public static SparseMatrix pow(double A, SparseMatrix B) {
	SparseMatrix C = new SparseMatrix(B.m_rows, B.m_cols);
	for(Integer i: B.m_matrix.keySet()) {
	    HashMap<Integer, Double> diagB = B.m_matrix.get(i);
	    HashMap<Integer, Double> diagC = C.getDiagonalFailsafe(i);
	    for(Integer j: diagB.keySet()) diagC.put(j, Math.pow(A, diagB.get(j)));
	}
	return C;
    }
    
    public static SparseMatrix createHorizontal(SparseMatrix[] A) {
	int col_offset = 0;
	SparseMatrix C = new SparseMatrix();
	for(int i = 0; i < A.length; i++) {
	    C.overlay(A[i], 0, col_offset);
	    col_offset += A[i].m_cols;
	}
	return C;
    }

    public static SparseMatrix createVertical(SparseMatrix[] A) {
	int row_offset = 0;
	SparseMatrix C = new SparseMatrix();
	for(int i = 0; i < A.length; i++) {
	    C.overlay(A[i], row_offset, 0);
	    row_offset += A[i].m_rows;
	}
	return C;
    }

    public static SparseMatrix createDiagonal(SparseMatrix[] A) {
	int row_offset = 0;
	int col_offset = 0;
	SparseMatrix C = new SparseMatrix();
	for(int i = 0; i < A.length; i++) {
	    C.overlay(A[i], row_offset, col_offset);
	    row_offset += A[i].m_rows;
	    col_offset += A[i].m_cols;
	}
	return C;
    }

    //////////////////////////////////////////////////////////////////////////////

    private void overlay(SparseMatrix A, int row_offset, int col_offset) {
	m_rows = Math.max(m_rows, A.m_rows + row_offset);
	m_cols = Math.max(m_cols, A.m_cols + col_offset);

	for(Integer i: A.m_matrix.keySet()) {
	    HashMap<Integer, Double> diagA = A.m_matrix.get(i);
	    for(Integer j: diagA.keySet()) {
		Double value = diagA.get(j);
		int rp = getRow(i,j)+row_offset;
		int cp = getCol(i,j)+col_offset;
		int ip = getDiag(rp, cp);
		int jp = getIdx(rp, cp);
		setDiagIdx(ip, jp, value);
	    }
	}
    }

    private boolean isValidDiagonal(int i) {
	if(i < 0 && -i >= m_cols) return false;
	if(i >=0 &&  i >= m_rows) return false;
	return true;
    }

    private boolean isValidIdx(int i, int j) {
	int row = getRow(i, j);
	if(row < 0 || row >= m_rows) return false;
	int col = getCol(i, j);
	if(col < 0 || col >= m_cols) return false;
	return true;
    }
    
    private static int getWidth(int i, int rows, int cols) {
	return i > 0 ? Math.min(rows-i,cols) : Math.min(rows, cols+i);
    }

    private HashMap<Integer, Double> getDiagonalFailsafe(int i) {
	if(!m_matrix.containsKey(i)) m_matrix.put(i, new HashMap<Integer, Double>());
	return m_matrix.get(i);
    }

    public Double getLinear(int k) {
	int r = k / m_cols;
	int c = k % m_cols;
	int i = getDiag(r, c);
	int j = getIdx (r, c);
	return get(i, j);
    }

    public Double[] get(int r) {
	Double[] row = new Double[m_cols];
	for(int c = 0; c < m_cols; c++) row[c] = get(getDiag(r,c), getIdx(r,c));
	return row;
    }

    private Double get(int i, int j) {
	if(!m_matrix.containsKey(i)) return 0.0;
	HashMap<Integer, Double> diag = m_matrix.get(i);
	if(!diag.containsKey(j)) return 0.0;
	return diag.get(j);
    }

    public void set(int k, Double value) {   // number of columns must already be known
	int r = k / m_cols;
	int c = k % m_cols;
	set(r, c, value);
    }

    public void set(int r, int c, Double value) { // (rows, cols) are updated.
	m_rows = Math.max(r+1, m_rows);
	m_cols = Math.max(c+1, m_cols);
	int i = getDiag(r,c);
	int j = getIdx (r,c);
	setDiagIdx(i, j, value);
    }

    private void setDiagIdx(int i, int j, Double value) {  // actual set() function.
	if(m_values != null) resetLinearAccess();
	if(value.equals(0.0)) return;                // no zeros allowed or needed.
	HashMap<Integer, Double> diag = getDiagonalFailsafe(i);
	diag.put(j, value);
    }

    private void populateDiagonal(int i) {  // alternative set() function
	if(m_values != null) resetLinearAccess();
	HashMap<Integer, Double> diag = getDiagonalFailsafe(i);
	int n = getWidth(i, m_rows, m_cols);
	int s = i < 0 ? -i : i;
	for(int j = 0; j < n; j++) diag.put(2*j+s, 1.0);  // Assume a 1
    }

    private void populateDiagonal(int i, double[] v) {  // alternative set() function
	if(m_values != null) resetLinearAccess();
	HashMap<Integer, Double> diag = getDiagonalFailsafe(i);
	int n = getWidth(i, m_rows, m_cols);
	int s = i < 0 ? -i : i;
	for(int j = 0; j < n; j++) diag.put(2*j+s, v[j]); // assume v large enough
    }

    private void accum(int i, int j, Double value) {
	Double old_value = get(i, j);
	setDiagIdx(i, j, value + old_value);
    }

    public boolean hasDiagonal(int i) {
	return m_matrix.containsKey(i);
    }

    // The Diag is the distance along the anti-diagonal, the Idx is distance from anti-diagonal
    private static int getRow (int i, int j) { return (j + i)/2; }
    private static int getCol (int i, int j) { return (j - i)/2; }
    private static int getDiag(int r, int c) { return r - c; }
    private static int getIdx (int r, int c) { return r + c; }  

    public Double[][] getDenseMatrix() {
	Double[][] matrix = new Double[m_rows][];
	for(int row = 0; row < m_rows; row++) {
	    Double[] Row = new Double[m_cols];
	    matrix[row] = Row;
	    for(int col = 0; col < m_cols; col++) Row[col] = 0.0;
	}
	for(Integer i: m_matrix.keySet()) {
	    HashMap<Integer, Double> diag = m_matrix.get(i);
	    for(Integer j: diag.keySet()) {
		//System.out.println("M("+i+", "+j+")["+getRow(i,j)+
		//		   ", "+getCol(i,j)+"] = "+diag.get(j));
		matrix[getRow(i,j)][getCol(i,j)] = diag.get(j);
	    }
	}
	return matrix;
    }

    public String toString() {
	String ZERO  = "   "; // TEMP
	StringBuffer line = new StringBuffer();
	Double[][] matrix = getDenseMatrix();
	for(int row = 0; row < m_rows; row++) {
	    StringBuffer item = new StringBuffer();
	    for(int col = 0; col < m_cols; col++) {
		Double value = matrix[row][col];
		if(value.equals(0.0)) item.append(ZERO);
		else                  item.append(value.toString());
		if(col < m_cols - 1) item.append(SPACE);
	    }
	    line.append(item.toString());
	    if(row < m_rows - 1) line.append(NEWLINE);
	}
	return line.toString();
    }
}
