package com.tomfielden.util;

import org.python.core.*;

import java.util.*;
import java.io.*;

import org.gnu.glpk.*;

/* Jython Gass Example
from com.tomfielden.util import *
p  = (45,80)
A  = ((5,20), (10,15))
b  = (400,450)
LP = Glpk(p, A, b).maximize()
# Objective = 2200.0: (24.0, 14.0)
 */

public class Glpk {

    static { System.loadLibrary("glpk_java"); }

    static private final String OPEN_PAREN    = "(";
    static private final String CLOSE_PAREN   = ")";
    static private final String COMMA_SPACE   = ", ";
    static private final String COLON_SPACE   = ": ";

    private Variant    m_p;
    private Variant    m_A;
    private Variant    m_b;

    private int        m_direction;
    private int        m_return;
    private double     m_objective;
    private double[]   m_x;

    public Glpk(Variant p, Variant A, Variant b) {
	m_p = p;
	m_A = A;
	m_b = b;
    }

    public Glpk(PyObject p, PyObject A, PyObject b) {
	m_p = new Variant(p);
	m_A = new Variant(A);
	m_b = new Variant(b);
    }

    public Glpk maximize() {
	m_direction = GLPKConstants.GLP_MAX;
	run(m_p.getDoubleArray(), m_A.getSparseMatrix(), m_b.getDoubleArray());
	return this;
    }

    public Glpk minimize() {
	m_direction = GLPKConstants.GLP_MIN;
	run(m_p.getDoubleArray(), m_A.getSparseMatrix(), m_b.getDoubleArray());
	return this;
    }

    private void run(double[] p, SparseMatrix A, double[] b) {
	// NB: GLPK is 1-based indexing. This is a pain and wastes the 0th element of arrays.

	int   [] irows = A.getRowIndicies();
	int   [] icols = A.getColIndicies();
	double[] ivals = A.getValues();      // NO ZEROs ALLOWED!

	int N    = irows.length;             // irows, icols, ivals all same length. Parallel.
	int rows = A.rows();
	int cols = A.cols();

	if(cols != p.length) {
	    System.out.println("Glpk::p " + cols +" != "+p.length);
	    return;
	}
	    
	if(rows != b.length) {
	    System.out.println("Glpk::b " + rows +" != "+b.length);
	    return;
	}

	glp_prob lp = GLPK.glp_create_prob();

	SWIGTYPE_p_int    ia = GLPK.new_intArray   (N+1); // rows
	SWIGTYPE_p_int    ja = GLPK.new_intArray   (N+1); // cols
	SWIGTYPE_p_double ar = GLPK.new_doubleArray(N+1); // values

	for(int i = 0; i < N; i++) {
	    // System.out.println("ia["+(i+1)+"] = "+(irows[i]+1));    // DEBUG
	    GLPK.intArray_setitem   (ia, i+1, irows[i]+1); 
	}

	for(int i = 0; i < N; i++) {
	    //System.out.println("ja["+(i+1)+"] = "+(icols[i]+1));    // DEBUG
	    GLPK.intArray_setitem   (ja, i+1, icols[i]+1); 
	}
	for(int i = 0; i < N; i++) {
	    //System.out.println("ar["+(i+1)+"] = "+ivals[i]);        // DEBUG
	    GLPK.doubleArray_setitem(ar, i+1, ivals[i]);
	}

	GLPK.glp_set_obj_dir(lp, m_direction);

	GLPK.glp_add_rows(lp, rows);
	for(int i = 0; i < rows; i++) {
	    //System.out.println("addRows["+(i+1)+"] = " + b[i]);      // DEBUG
	    GLPK.glp_set_row_bnds(lp, i+1, GLPKConstants.GLP_UP, 0.0, b[i]);
	}

	// Require 0 <= X if GLP_LO, free if GLP_FR 
	GLPK.glp_add_cols(lp, cols);
	for(int i = 0; i < cols; i++) {
	    //System.out.println("addCols["+(i+1)+"] = " + p[i]);     // DEBUG
	    //GLPK.glp_set_col_bnds(lp, i+1, GLPKConstants.GLP_LO, 0.0, 0.0);
	    GLPK.glp_set_col_bnds(lp, i+1, GLPKConstants.GLP_FR, 0.0, 0.0);
	    GLPK.glp_set_obj_coef(lp, i+1, p[i]);
	}
	    
	GLPK.glp_load_matrix(lp, N, ia, ja, ar);
	GLPK.glp_simplex(lp, null);               // <<<<< THIS IS IT! THE HEART OF THE SYSTEM
	    
	m_objective = GLPK.glp_get_obj_val(lp);
	
	m_x = new double[cols];
	for(int i = 0; i < cols; i++) m_x[i] = GLPK.glp_get_col_prim(lp, i+1);
	
	GLPK.glp_delete_prob(lp);
    }

    public double[] getX() { return m_x; }

    public String toString() {
	StringBuffer buf = new StringBuffer();
	buf.append("Objective = " + m_objective + COLON_SPACE);
	buf.append(OPEN_PAREN);
	for(int i = 0; i < m_x.length; i++) {
	    buf.append(""+m_x[i]);
	    if(i < m_x.length - 1) buf.append(COMMA_SPACE);
	}
	buf.append(CLOSE_PAREN);
	return buf.toString();
    }


    /////////////////////////////////////////////////////////////////////////////////////////

    public static void main(String[] args) {
	int    size = 1000+1;
	SWIGTYPE_p_int    ia = GLPK.new_intArray(size);
	SWIGTYPE_p_int    ja = GLPK.new_intArray(size);
	SWIGTYPE_p_double ar = GLPK.new_doubleArray(size);

	glp_prob lp   = GLPK.glp_create_prob();

	GLPK.glp_set_prob_name(lp, "sample");
	GLPK.glp_set_obj_dir(lp, GLPKConstants.GLP_MAX);
	GLPK.glp_add_rows(lp, 3);
	GLPK.glp_set_row_name(lp, 1, "p");
	GLPK.glp_set_row_bnds(lp, 1, GLPKConstants.GLP_UP, 0.0, 100.0);
	GLPK.glp_set_row_name(lp, 2, "q");
	GLPK.glp_set_row_bnds(lp, 2, GLPKConstants.GLP_UP, 0.0, 600.0);
	GLPK.glp_set_row_name(lp, 3, "r");
	GLPK.glp_set_row_bnds(lp, 3, GLPKConstants.GLP_UP, 0.0, 300.0);

	GLPK.glp_add_cols(lp, 3);
	GLPK.glp_set_col_name(lp, 1, "x1");
	GLPK.glp_set_col_bnds(lp, 1, GLPKConstants.GLP_LO, 0.0, 0.0);
	GLPK.glp_set_obj_coef(lp, 1, 10.0);
	GLPK.glp_set_col_name(lp, 2, "x2");
	GLPK.glp_set_col_bnds(lp, 2, GLPKConstants.GLP_LO, 0.0, 0.0);
	GLPK.glp_set_obj_coef(lp, 2, 6.0);
	GLPK.glp_set_col_name(lp, 3, "x3");
	GLPK.glp_set_col_bnds(lp, 3, GLPKConstants.GLP_LO, 0.0, 0.0);
	GLPK.glp_set_obj_coef(lp, 3, 4.0);

	// GLPK.intArray_setitem(SWIGTYPE_p_int ary, int index, int value) ;
	// GLPK.doubleArray_setitem(SWIGTYPE_p_double ary, int index, double value);

	GLPK.intArray_setitem   (ia, 1,  1); 
	GLPK.intArray_setitem   (ja, 1,  1); 
	GLPK.doubleArray_setitem(ar, 1, 1.0);      /* a[1,1] =  1 */
	
	GLPK.intArray_setitem   (ia, 2,  1); 
	GLPK.intArray_setitem   (ja, 2,  2); 
	GLPK.doubleArray_setitem(ar, 2, 1.0);      /* a[1,2] =  1 */

	GLPK.intArray_setitem   (ia, 3,  1); 
	GLPK.intArray_setitem   (ja, 3,  3); 
	GLPK.doubleArray_setitem(ar, 3, 1.0);      /* a[1,3] =  1 */

	GLPK.intArray_setitem   (ia, 4,  2); 
	GLPK.intArray_setitem   (ja, 4,  1); 
	GLPK.doubleArray_setitem(ar, 4, 10.0);      /* a[2,1] = 10 */
	
	GLPK.intArray_setitem   (ia, 5,  3); 
	GLPK.intArray_setitem   (ja, 5,  1); 
	GLPK.doubleArray_setitem(ar, 5, 2.0);      /* a[3,1] =  2 */

	GLPK.intArray_setitem   (ia, 6,  2); 
	GLPK.intArray_setitem   (ja, 6,  2); 
	GLPK.doubleArray_setitem(ar, 6, 4.0);      /* a[2,2] =  4 */

	GLPK.intArray_setitem   (ia, 7,  3); 
	GLPK.intArray_setitem   (ja, 7,  2); 
	GLPK.doubleArray_setitem(ar, 7, 2.0);      /* a[3,2] =  2 */
	
	GLPK.intArray_setitem   (ia, 8,  2); 
	GLPK.intArray_setitem   (ja, 8,  3); 
	GLPK.doubleArray_setitem(ar, 8, 5.0);      /* a[2,3] =  5 */

	GLPK.intArray_setitem   (ia, 9,  3); 
	GLPK.intArray_setitem   (ja, 9,  3); 
	GLPK.doubleArray_setitem(ar, 9, 6.0);      /* a[3,3] =  6 */

	GLPK.glp_load_matrix(lp, 9, ia, ja, ar);
	GLPK.glp_simplex(lp, null);
	double Z = GLPK.glp_get_obj_val(lp);
	double x1 = GLPK.glp_get_col_prim(lp, 1);
	double x2 = GLPK.glp_get_col_prim(lp, 2);
	double x3 = GLPK.glp_get_col_prim(lp, 3);
	System.out.println("\nZ = "+Z+" x1 = "+x1+" x2 = "+x2+" x3 = "+x3);
	GLPK.glp_delete_prob(lp);

	System.out.println("Finished GLPK test.");
    }

}