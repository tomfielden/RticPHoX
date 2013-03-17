package com.tomfielden.util;

import org.python.core.*;
import org.apache.commons.math.util.*;

import java.util.*;
import java.io.*;

import com.tomfielden.random.LazyRandomVariable;
import com.tomfielden.random.NumericRandomVariable;
import com.tomfielden.random.RandomSample;

@SuppressWarnings("unchecked")
public class Variant {

    static private final String OPEN_BRACE    = "{";
    static private final String CLOSE_BRACE   = "}";
    static private final String OPEN_PAREN    = "(";
    static private final String CLOSE_PAREN   = ")";
    static private final String OPEN_BRACKET  = "[";
    static private final String CLOSE_BRACKET = "]";
    static private final String COMMA_SPACE   = ", ";
    static private final String COLON_SPACE   = ": ";
    static private final String NEWLINE       = "\n";

    ////////////////////////////////////////////////////////////////////////

    public enum VTYPE {NULL, BOOL, INT, DOUBLE, STRING, MASK, ARRAY, INDEX, INT_MAP, 
	    STRING_MAP, INT_SET, STRING_SET, SLICE, VECTOR, SPMATRIX, FUNCTION, LRV, NRV, NRS};

    public VTYPE  m_type;
    public Object m_object;

    ///////////////////////////////////////////////////////////////////////////////

    public Variant() {
	m_object = null;
	m_type   = VTYPE.NULL;
    }

    // Jython doesn't see this, sets Integer even if boolean. Used internally.
    public Variant(Boolean x) {
	m_object = x;
	m_type = VTYPE.BOOL;
    }

    public Variant(Integer x) {
	m_object = x;
	m_type = VTYPE.INT;
    }

    public Variant(int x) {
	m_object = new Integer(x);
	m_type = VTYPE.INT;
    }

    public Variant(Double x) {
	m_object = x;
	m_type   = VTYPE.DOUBLE;
    }

    public Variant(double x) {
	m_object = new Double(x);
	m_type   = VTYPE.DOUBLE;
    }

    public Variant(String x) {
	m_object = x;
	m_type   = VTYPE.STRING;
    }

    // NOT a copy constructor. See Functions.java for reason. Obscures new construction and renaming.
    public Variant(Variant x) {
	setVariant(x);
    }

    private void setVariant(Variant x) {
	m_object = x.m_object;
	m_type   = x.m_type;
    }

    public Variant(Integer[] x) {
	m_object = x;
	m_type   = VTYPE.INDEX;
    }

    /* NB: Can't have this method since Python will stuff integer arrays in here.
    public Variant(Double[] x) {

	System.out.println("Variant::Double[] " + x);  // DEBUG

	Variant[] v = new Variant[x.length];
	for(int i = 0; i < x.length; i++) v[i] = new Variant(x[i]);
	m_object = v;
	m_type   = VTYPE.ARRAY;
    }
    */

    private Variant(Boolean[] x) {
	m_object = x;
	m_type   = VTYPE.MASK;
    }

    public Variant(LazyRandomVariable x) {
	m_object = x;
	m_type   = VTYPE.LRV;
    }

    public Variant(NumericRandomVariable x) {
	m_object = x;
	m_type   = VTYPE.NRV;
    }

    public Variant(RandomSample x) {
	m_object = x;
	m_type   = VTYPE.NRS;
    }

    ///////////////////////////////////////////////////////////////////////
    // Static Creators

    // NB: used to accept Double[]. Not sure the impact. Jython seems OK
    public static Variant createDoubleArray(double[] x) {
	Variant[] v = new Variant[x.length];
	for(int i = 0; i < x.length; i++) v[i] = new Variant(x[i]);
	Variant y = new Variant();
	y.m_object = v;
	y.m_type   = VTYPE.ARRAY;
	return y;
    }

    private static Variant createDoubleArray(Double[] x) {
	Variant[] v = new Variant[x.length];
	for(int i = 0; i < x.length; i++) v[i] = new Variant(x[i]);
	Variant y = new Variant();
	y.m_object = v;
	y.m_type   = VTYPE.ARRAY;
	return y;
    }

    public static Variant createDoubleArray(double[][] x) {
	Variant[] v = new Variant[x.length];
	for(int i = 0; i < x.length; i++) v[i] = Variant.createDoubleArray(x[i]);
	Variant y = new Variant();
	y.m_object = v;
	y.m_type   = VTYPE.ARRAY;
	return y;
    }

    public static Variant linspace(Variant start, Variant stop, Variant count) {
	return linspace(start.toDouble(), stop.toDouble(), count.toInteger());
    }

    public static Variant linspace(double start, double stop, int count) {
	double[] x = new double[count];
	double delta  = (stop-start)/(count-1);
	for(int i = 0; i < count; i++) x[i] = start + delta*i;
	return createDoubleArray(x);
    }

    public static Variant createIntMap(Variant A, Variant B) {
	// Create consistent map from A to B.
	HashMap<Integer,  Variant> int_map = new HashMap<Integer,  Variant>();
	for(int i = 0; i < A.size(); i++) int_map.put(A.get(i).toInteger(), B.get(i));
	return new Variant(int_map, VTYPE.INT_MAP);
    }

    public static Variant createStringMap() {
	HashMap<String,  Variant> str_map = new HashMap<String,  Variant>();
	return new Variant(str_map, VTYPE.STRING_MAP);
    }

    public static Variant createStringMap(HashMap<String, Variant> str_map) {
	return new Variant(str_map, VTYPE.STRING_MAP);
    }

    public static Variant createStringMap(Variant A, Variant B) {
	// Create consistent map from A to B.
	HashMap<String,  Variant> str_map = new HashMap<String,  Variant>();
	for(int i = 0; i < A.size(); i++) str_map.put(A.get(i).toString(), B.get(i));
	return new Variant(str_map, VTYPE.STRING_MAP);
    }

    public static Variant createIntSet() {
	HashSet<Integer> set = new HashSet<Integer>();
	return new Variant(set, VTYPE.INT_SET);
    }

    public static Variant createIntSet(Variant B) {
	Variant A = createIntSet();
	A.append(B);
	return A;
    }

    public static Variant createStringSet() {
	HashSet<String> set = new HashSet<String>();
	return createStringSet(set);
    }

    public static Variant createStringSet(HashSet<String> set) {
	return new Variant(set, VTYPE.STRING_SET);
    }

    public static Variant createStringSet(Variant B) {
	Variant A = createStringSet();
	A.append(B);
	return A;
    }

    public static Variant createSlice(Variant s) {
	return createSlice(s.toString());
    }

    public static Variant createSlice(String s) {
	Integer[] ii = new Integer[3];
	Variant   v  = new Variant();
	v.m_object   = ii;
	v.m_type     = VTYPE.SLICE;

	ii[0] = 0;   // start at beginning
	ii[1] = 0;   // means "the end". Should be 0?
	ii[2] = 1;   // increment

	String[] S = s.split(":");
	if(S.length == 1) {
	    int d = 0;
	    try {
		if(S[0].length() > 0) d = Integer.parseInt(S[0]);
	    } catch(Exception e) {
		System.out.println("Variant::createSlice failed for: " + s + "\n" + e);
		return new Variant();
	    }
	    int index = s.indexOf(':');
	    if(index == 0) ii[1] = d;
	    else           ii[0] = d;
	} else if(S.length == 2) {
	    try {
		if(S[0].length() > 0) ii[0] = Integer.parseInt(S[0]);
		if(S[1].length() > 0) ii[1] = Integer.parseInt(S[1]);
	    } catch(Exception e) {
		System.out.println("Variant::createSlice failed for: " + s + "\n" + e);
		return new Variant();
	    }
	} else if(S.length == 3) {
	    try {
		if(S[0].length() > 0) ii[0] = Integer.parseInt(S[0]);
		if(S[1].length() > 0) ii[1] = Integer.parseInt(S[1]);
		if(S[2].length() > 0) ii[2] = Integer.parseInt(S[2]);
	    } catch(Exception e) {
		System.out.println("Variant::createSlice failed for: " + s + "\n" + e);
		return new Variant();
	    }
	} else if(S.length > 3) {
	    System.out.println("Variant::createSlice unsupported slice: " + s);
	    return new Variant();
	}

	return v;
    }

    public static Variant createMask(Boolean[] x) {
	return new Variant(x);
    }

    public static Variant createMask(Object o) {
	return createMask(new Variant(o));
    }

    public static Variant createMask(int n) {
	return createMask(n, true);
    }

    public static Variant createMask(int n, boolean val) {
	Boolean[] B = new Boolean[n];
	for(int i = 0; i < n; i++) B[i] = val;
	return new Variant(B);
    }

    public static Variant createMask(Variant x) {
	switch(x.m_type) {
	case BOOL:
	    Boolean[] B3 = new Boolean[1];
	    B3[0]        = x.isTrue();
	    return new Variant(B3);
	case INT:
	case DOUBLE:
	    return createMask(x.toInteger(), true);    // default is "true"
	case MASK:
	    return x.copy();
	case INDEX:                                    // ARRAY-like
	    Integer[] X4 = (Integer[])x.m_object;
	    Boolean[] B4 = new Boolean[X4.length];
	    for(int i = 0; i < X4.length; i++) B4[i] = X4[i] > 0;
	    return new Variant(B4);
	case ARRAY:
	    Variant[] X2 = (Variant[])x.m_object;
	    Boolean[] B2 = new Boolean[X2.length];
	    for(int i = 0; i < X2.length; i++) B2[i] = X2[i].isTrue();
	    return new Variant(B2);
	}
	System.out.println("Variant::createMask not supported for type: " + x.m_type);
	return new Variant();
    }

    // Used to create mask based on a query of B into A. 
    public static Variant createMask(Variant A, Variant B) {
	return createMaskV(A, B);
    }

    private static Variant createMaskV(Variant A, Variant B) {
	// B contains generalized key-values.
	switch(A.m_type) {
	case ARRAY:
	    switch(B.m_type) {
	    case STRING_SET:
		Variant mask = new Variant();
		for(String key: B.toStringSet()) 
		    mask = add(mask,cmp(A, new Variant(key), CMP_FLAVOR.EQ));
		return mask;
	    }
	    break;
	case STRING_MAP:
	    HashMap<String, Variant> Amap = (HashMap<String, Variant>)A.m_object;
	    switch(B.m_type) {
	    case STRING_MAP:
		HashMap<String, Variant> Bmap = (HashMap<String, Variant>)B.m_object;
		Variant mask = new Variant();
		for(String Bkey: Bmap.keySet()) {  // loop over all elements in outer MAP
		    //System.out.println("Bkey: " + Bkey); // DEBUG
		    Variant C = Bmap.get(Bkey);    // the element in the MAP, peek at type.
		    Variant mask2 = new Variant();
		    switch(C.m_type) {
		    case STRING:
			mask2 = cmp(Amap.get(Bkey), C, CMP_FLAVOR.EQ);
			//System.out.println("mask2.STRING " + mask2);  // DEBUG
			break;
		    case ARRAY:
			mask2 = createMaskV(Amap.get(Bkey), C.getStringSet());
			// System.out.println("mask2.ARRAY " + mask2); // DEBUG 
			break;
		    case STRING_SET:
			mask2 = createMaskV(Amap.get(Bkey), C);
			//System.out.println("mask2.STRING_SET " + mask2); // DEBUG
			break;
		    case STRING_MAP:
			HashMap<String, Variant> Cmap = (HashMap<String, Variant>)C.m_object;
			for(String Ckey: Cmap.keySet()) {
			    Variant D = Cmap.get(Ckey);
			    //System.out.println("(inner)" + Ckey + ", " + D); // DEBUG
			    Variant mask3 = cmp(Amap.get(Bkey),new Variant(Ckey),CMP_FLAVOR.EQ);
			    if(D.m_type == VTYPE.STRING_MAP) 
				mask2 = add(mask2,mult(mask3,createMaskV(A, D)));
			    else 
				mask2 = add(mask2, mask3);
			}
			break;
		    }
		    mask = add(mask, mask2);  // outer loop OR
		}
		return mask;
	    }
	    break;
	}

	System.out.println("unhandled createMask(" + A.m_type +", "+ B.m_type + ")"); // DEBUG
	return new Variant();
    }

    public static Variant createPartition(Variant A) {
	// An array of ordered parallel non-intersecting masks.
	// uses .toInteger() to find individual sizes. Allows "examples" instead of just ints.
	switch(A.m_type) {
	case ARRAY:
	    Variant[] A2 = (Variant[])A.m_object;
	    Variant[] M  = new Variant[A2.length];   // mask for A[i]
	    int    [] n  = new int    [A2.length];   // size of A[i]
	    int    [] a  = new int    [A2.length+1]; // cumsize of A[i] shifted by 1 so a can
	    a      [0]   = 0;                        // now serve as starting index into M[i]
	    int       N = 0;
	    for(int i = 0; i < A2.length; i++) {n[i] = A2[i].toInteger(); N += n[i]; a[i+1] = N;}
	    for(int i = 0; i < A2.length; i++) {
		M[i] = createMask(N, false);
		Boolean[] b = (Boolean[])M[i].m_object;
		for(int j = a[i]; j < a[i+1]; j++) b[j] = true;
	    }
	    return new Variant(M);
	}
	System.out.println("Variant::createMaskSet["+A.m_type+"] unsupported.");
	return new Variant();
    }

    public static Variant createIndex(Variant x) {
	switch(x.m_type) { 
	case INDEX: 
	    return x.copy();
	case INT:
	    Integer[] I = new Integer[1];
	    I[0] = x.toInteger();
	    return new Variant(I);
	case STRING:
	case MASK:
	case VECTOR:
	case ARRAY:
	case INT_MAP:     
	case STRING_MAP:
	    return new Variant(x.getIndex());
	}
	System.out.println("Variant::createIndex("+x.m_type+") unsupported.");
	return new Variant();
    }

    public static Variant createVector(Object o) {
	return createVector(new Variant(o));
    }

    public static Variant createVector(Double[] x) {
	Variant v  = new Variant();
	v.m_object = x;
	v.m_type   = VTYPE.VECTOR;
	return v;
    }

    public static Variant createVector(Variant x) {
	Variant v  = new Variant();
	v.m_object = x.getVector();
	v.m_type   = VTYPE.VECTOR;
	return v;
    }

    public static Variant createSparseMatrix(int rows, int cols) {
	return createSparseMatrix(rows, cols, 0);
    }

    public static Variant createSparseMatrix(int rows, int cols, int sideband) {
	return createSparseMatrix(new SparseMatrix(rows, cols, sideband));
    }

    public static Variant createSparseMatrix(Variant rows, Variant cols, Variant sideband) {
	switch(rows.m_type) {
	case INT:
	    Integer irows = (Integer)rows.m_object;
	    switch(cols.m_type) {
	    case INT:
		Integer icols = (Integer)cols.m_object;
		switch(sideband.m_type) {
		case INT:
		case DOUBLE:
		    return createSparseMatrix(irows, icols, sideband.toInteger());
		case ARRAY:
		    Variant[] S = (Variant[])sideband.m_object;
		    Variant[] C = new Variant[S.length];
		    for(int i = 0; i < S.length; i++) C[i] = createSparseMatrix(rows,cols,S[i]);
		    return new Variant(C);
		case VECTOR:
		    return createSparseMatrix(rows, cols, sideband.toArray());
		}
		break;
	    }
	    break;
	}
	System.out.println("createSparseMatrix("+rows.m_type+", "+
			   cols.m_type+", "+sideband.m_type+") unsupported.");
	return new Variant();
    }

    public static Variant createSparseMatrix(SparseMatrix A) {
	Variant V  = new Variant();
	V.m_object = A;
	V.m_type   = VTYPE.SPMATRIX;
	return V;
    }

    public static Variant createSparseDiagonal(Variant A) {
	switch(A.m_type) {
	case INT:
	case DOUBLE:
	    int s = A.toInteger();
	    return createSparseMatrix(new SparseMatrix(s,s));
	case INDEX:
	case VECTOR:
	case ARRAY:
	    Double[] v = A.getVector();
	    double[] v2 = new double[v.length];
	    for(int i = 0; i < v.length; i++) v2[i] = v[i];
	    return createSparseMatrix(new SparseMatrix(v2));
	}
	System.out.println("Varaint::createSparseDiagonal["+A.m_type+"] unsupported.");
	return new Variant();
    }

    public static Variant createSparseDiagonal(Variant A, Variant B) {
	switch(A.m_type) {
	case INT:
	case DOUBLE:
	    switch(B.m_type) {
	    case INDEX:
	    case ARRAY:
		SparseMatrix C = new SparseMatrix();
		C.setDiagonal(A.toInteger(), B.getVector());
		return new Variant(C);
	    }
	    break;
	case INDEX:
	    switch(B.m_type) {
	    case INDEX:
		SparseMatrix C2 = new SparseMatrix();
		Integer[]    A2 = (Integer[])A.m_object;
		for(int i = 0; i < A2.length; i++) C2.setDiagonal(A2[i], B.getVector());
		return new Variant(C2);
	    case ARRAY:
		SparseMatrix C1 = new SparseMatrix();
		Integer[]    A1 = (Integer[])A.m_object;
		for(int i = 0; i < A1.length; i++) C1.setDiagonal(A1[i], B.get(i).getVector());
		return new Variant(C1);
	    }
	    break;
	}
	System.out.println("Varaint::createSparseDiag["+A.m_type+", "+B.m_type+"] unsupported.");
	return new Variant();
    }

    private enum SPARSEFLAVOR {HORIZONTAL, VERTICAL, DIAGONAL};

    public static Variant createSparseMatrixHorizontal(Variant A) {
	return createSparseMatrixByFlavor(A, SPARSEFLAVOR.HORIZONTAL);
    }

    public static Variant createSparseMatrixVertical(Variant A) {
	return createSparseMatrixByFlavor(A, SPARSEFLAVOR.VERTICAL);
    }

    public static Variant createSparseMatrixDiagonal(Variant A) {
	return createSparseMatrixByFlavor(A, SPARSEFLAVOR.DIAGONAL);
    }

    private static Variant createSparseMatrixByFlavor(Variant A, SPARSEFLAVOR F) {
	switch(A.m_type) {
	case ARRAY:
	    Variant[] A2 = (Variant[])A.m_object;
	    SparseMatrix[] C = new SparseMatrix[A2.length];
	    for(int i = 0; i < A2.length; i++) {
		switch(A2[i].m_type) {
		case SPMATRIX:
		    C[i] = (SparseMatrix)(A2[i].m_object);
		    break;
		default:
		    C[i] = new SparseMatrix();
		}
	    }
	    switch(F) {
	    case HORIZONTAL:
		return new Variant(SparseMatrix.createHorizontal(C));
	    case VERTICAL:
		return new Variant(SparseMatrix.createVertical(C));
	    case DIAGONAL:
		return new Variant(SparseMatrix.createDiagonal(C));
	    }
	}
	System.out.println("createSparseMatrixHorizontal("+A.m_type+") unsupported.)");
	return new Variant();
    }

    public Variant(String[] x) {
	Variant[] v = new Variant[x.length];
	for(int i = 0; i < x.length; i++) v[i] = new Variant(x[i]);
	m_object = v;
	m_type   = VTYPE.ARRAY;
    }

    // We can't expose this function else Jython will use it instead of Integer[]
    // This means that Boolean ARRAY's can only be generated indirectly. That's probably OK.
    /*
    public Variant(Boolean[] x) {
	System.out.println("Variant(Boolean[])");
	Variant[] v = new Variant[x.length];
	for(int i = 0; i < x.length; i++) v[i] = new Variant(x[i]);
	m_object = v;
	m_type   = VTYPE.ARRAY;
    }
    */

    public Variant(SparseMatrix x) {
	m_object = x;
	m_type   = VTYPE.SPMATRIX;
    }

    public Variant(Variant[] x) {
	m_object = x;
	m_type   = VTYPE.ARRAY; 
    }

    private Variant(Map<?, ?> x, VTYPE type) {
	m_object = x;
	m_type   = type;
    }

    private Variant(Set<?> x, VTYPE type) {
	m_object = x;
	m_type   = type;
    }

    private void setObjectArray(Object[] x) {
	Variant[] v = new Variant[x.length];
	for(int i = 0; i < x.length; i++) v[i] = new Variant(x[i]);
	m_object = v;
	m_type   = VTYPE.ARRAY; 
    }

    ///////////////////////////////////////////////////////////////////////////////

    public static Variant createObject(Object o) {
	Variant v = new Variant();
	v.setObject(o);
	return v;
    }

    private Variant(Object o) {
	//System.out.println("Variant::Object called with " + o.getClass().getName()); // DEBUG
	setObject(o);
    }

    private void setObject(Object o) {
	// System.out.println("setObject type: " + o.getClass().getName()); // DEBUG
	if(o instanceof Boolean) {
	    if((Boolean)o) m_object = 1;
	    else           m_object = 0;
	    m_type   = VTYPE.INT;
	} else if(o instanceof Integer) {
	    m_object = (Integer)o;
	    m_type   = VTYPE.INT;
	/* This is also dis-allowed
	} else if(o instanceof PyBoolean) {
	    System.out.println("setObject::PyBoolean"); // DEBUG
	    m_object = new Boolean(((PyBoolean)o).getValue()==1);
	    m_type   = VTYPE.BOOL;
	*/
	} else if(o instanceof PyInteger) {
	    m_object = ((PyInteger)o).getValue();
	    m_type   = VTYPE.INT;
	} else if(o instanceof Double) {
	    m_object = ((Double)o);
	    m_type   = VTYPE.DOUBLE;
	} else if(o instanceof String) { 
	    m_object = (String)o;
	    m_type   = VTYPE.STRING;
	} else if(o instanceof Variant) { 
	    setVariant((Variant)o);
	} else if(o instanceof SparseMatrix) {
	    m_object = (SparseMatrix)o;
	    m_type   = VTYPE.SPMATRIX;
	} else if(o instanceof Class) {
	    m_object = null;
	    m_type   = VTYPE.NULL;
	} else if(o instanceof PyFloat) {
	    m_object = ((PyFloat)o).getValue();
	    m_type   = VTYPE.DOUBLE;
	} else if(o instanceof PyLong) {
	    m_object = ((PyLong)o).doubleValue();
	    m_type   = VTYPE.DOUBLE;
	} else if(o instanceof PyComplex) {
	    PyComplex q = (PyComplex)o;
	    double    r = q.getReal().getValue();  // Assume nominal money
	    double    i = q.getImag().getValue();  // Assume financial year
	    System.out.println("Variant: Money not currently supported.");
	    m_object = null;
	    m_type   = VTYPE.NULL;
	    //setVariant(money(new Variant(r), new Variant(i)));
	} else if(o instanceof PyString) {
	    m_object = ((PyString)o).toString();
	    m_type   = VTYPE.STRING;
	} else if(o instanceof PyTuple) {
	    setObjectArray(((PyTuple)o).getArray());
	} else if(o instanceof PyList) {
	    setObjectArray(((PyList)o).toArray());
	} else if(o instanceof PyArray) {
	    setObjectArray(((PyList)((PyArray)o).tolist()).toArray());
	} else if(o instanceof PyDictionary) {
	    setDictionary((PyDictionary)o);
	} else if(o instanceof PyStringMap) {
	    setDictionary((PyStringMap)o);
	} else if(o instanceof PyReflectedFunction) {
	    // We picked up a function reference from a PyInterpreter. Just null it out quietly.
	    m_object = null;
	    m_type   = VTYPE.NULL;
	} else if(o instanceof PySet) {
	    PySet     s = (PySet)o;
	    int       i = 0;
	    Variant[] v = new Variant[s.size()];
	    for(Object x: s) v[i++] = new Variant(x);
	    Variant v2  = new Variant(v);
	    Variant v3  = v2.keys();
	    m_object    = v3.m_object;
	    m_type      = v3.m_type;
	} else if(o instanceof PySlice) {
	    // This is not very efficient. On the other hand it does exercise the code...
	    PySlice p = (PySlice)o;
	    String slice = "";
	    if(p.start.isNumberType()) slice = p.start + ":";
	    if(p.stop.isNumberType()) 
		if(slice.equals("")) slice = ":" + p.stop;
		else                 slice = slice + p.stop;
	    if(p.step.isNumberType())
		if(slice.equals("")) slice = "::" + p.step;
		else                 slice = slice + ":" + p.step;
	    Variant v = createSlice(slice);
	    m_object  = v.m_object;
	    m_type    = v.m_type;
	} else if(o instanceof PyObject) {
	    Object v = ((PyObject)o).__tojava__(Variant.class);  // magic step.
	    if(v instanceof Variant) setVariant((Variant)v);
	    else if(v instanceof PySingleton) { // This is failing. Not sure.
		PySingleton s = (PySingleton)v;
		System.out.println("Variant::setObject(PySingleton): "+(s.getType())); 
	    }
	    else {
		System.out.println("Variant::Object Weird Type: "+v.getClass().getName());
		System.out.println("  value = " + v);
	    }
	} else 
	    System.out.println("Variant::Object Unsupported Type: "+o.getClass().getName());
    }

    public Variant(PyObject p) {
	//System.out.println("Variant::PyObject called with " + p.getClass().getName()); //DEBUG
	Object o = p.__tojava__(Variant.class);  // magic step.

	if(o instanceof Variant) setVariant((Variant)o);
	else                      setObject ((Object)  p);
    }

    private void setDictionary(PyStringMap map) {
	HashMap<String, Variant> str_map = new HashMap<String,  Variant>();

	PyList items = map.items();
	for(int i = 0; i < items.size(); i++) {
	    PyTuple item = (PyTuple)items.get(i);
	    if(item.size() != 2) continue;   // should always be pairs

	    String key   = (String)(item.get(0));
	    Object value = item.get(1);

	    if(value == null)               continue;
	    if(value instanceof PyFunction) continue;
	    if(value instanceof Class)      continue;
	    if(key.length() >= 2 && key.substring(0,2).equals("__")) continue;
	    if(key.toUpperCase().equals(key)) continue;  // UPPER_CASE vars censored

	    str_map.put(key, new Variant(value));
	}

	m_type = VTYPE.STRING_MAP; 
	m_object = str_map;
    }

    private void setDictionary(PyDictionary dict) {
	// Watch for special cases
	boolean b_int_key = true;
	boolean b_num_val = true;

	HashMap<Integer, Variant> int_map = new HashMap<Integer, Variant>();
	HashMap<String,  Variant> str_map = new HashMap<String,  Variant>();

	Object[] keys = dict.keys().toArray();
	for(int i = 0; i < keys.length; i++) {
	    Object   key   = keys[i];
	    Variant value = new Variant(dict.get(key));
	    Variant vkey  = new Variant(key);
	    if(vkey.m_type != VTYPE.INT) b_int_key = false;
	    if(!value.isNumber())        b_num_val = false; 
	    if(b_int_key) int_map.put((Integer)vkey.m_object, value);
	    str_map.put(vkey.toString(), value);  // just in case.
	}

	if(b_int_key) {m_type = VTYPE.INT_MAP;    m_object = int_map;} 
	else          {m_type = VTYPE.STRING_MAP; m_object = str_map;}

	if(b_int_key && b_num_val) {
	    m_type   = VTYPE.INT_MAP;  // need to present an INT_MAP to the function
	    m_object = new LinearFunction(new Variant(this));
	    m_type   = VTYPE.FUNCTION; // and we need to be seen as a function.
	}
    }

    // Access ///////////////////////////////////////////////////////////////
    
    public VTYPE type() {return m_type;}

    public boolean isTrue() {
	switch(m_type) {
	case BOOL:   return ((Boolean)m_object).booleanValue();
	case INT:    return (Integer)m_object > 0;
	case DOUBLE: return (Double) m_object > 0.01;  // to allow for 0.0 ... 0.01 case.
	case STRING:
	    String s = (String)m_object;
	    return (s.equals("T")    || 
		    s.equals("t")    || 
		    s.equals("True") || 
		    s.equals("true") || 
		    s.equals("1"));
	case ARRAY:
	    Variant[] A = (Variant[])m_object;
	    for(int i = 0; i < A.length; i++) if(!A[i].isTrue()) return false;
	    return true;
	case MASK:
	    Boolean[] B2 = (Boolean[])m_object;
	    for(int i = 0; i < B2.length; i++) if(!B2[i]) return false;
	    return true;
	}
	System.out.println("Variant::isTrue["+m_type+"] unsupported.");
	return false;
    }

    public boolean isNaN() {
	switch(m_type) {
	case DOUBLE:
	    return ((Double)m_object).isNaN();
	case ARRAY:
	    Variant[] B = (Variant[])m_object;
	    for(int i = 0; i < B.length; i++) if(B[i].isNaN()) return true;
	    return false;
	default:
	    return false;
	}
    }

    public boolean isZero() {
	switch(m_type) {
	case BOOL:
	    return !isTrue();
	case INT:    
	    return (Integer)m_object == 0;
	case DOUBLE: 
	    Double d = (Double) m_object;
	    if(d < 0.000000000001 && d > -0.000000000001) return true;
	    return false;
	case MASK:
	    for(int i = 0; i < size(); i++) if(get(i).isTrue()) return false;
	    return true;
	}
	System.out.println("Variant::isZero unsupported type: " + m_type); 
	return false;
    }
    
    public boolean isNumber()            { return m_type == VTYPE.INT || m_type == VTYPE.DOUBLE;}
    public boolean isRV()                { return m_type == VTYPE.LRV || m_type == VTYPE.NRV;   }
    public boolean isRS()                { return m_type == VTYPE.NRS;                          }
    public boolean isGeneralizedNumber() { return isNumber()          || isRV();                }

    public boolean isNumbers() {
	switch(m_type) {
	case INT: 
	case DOUBLE:
	    return true;
	case ARRAY:
	    for(int i = 0; i < size(); i++) if(!get(i).isNumber()) return false;
	    return true;
	}
	return false;
    }

    public boolean isGeneralizedNumbers() {
	switch(m_type) {
	case INT: 
	case DOUBLE:
	    return true;
	case ARRAY:
	    for(int i = 0; i < size(); i++) if(!get(i).isGeneralizedNumber()) return false;
	    return true;
	}
	return false;
    }

    public boolean isIntegers() {
	switch(m_type) {
	case INT: return true;
	case ARRAY:
	    Variant[] vars = (Variant[])m_object;
	    for(int i = 0; i < vars.length; i++)
		if(vars[i].m_type != VTYPE.INT) return false;
	    return true;
	}
	return false;
    }

    public boolean isBooleans() {
	switch(m_type) {
	case BOOL: return true;
	case ARRAY:
	    Variant[] vars = (Variant[])m_object;
	    for(int i = 0; i < vars.length; i++)
		if(vars[i].m_type != VTYPE.BOOL) return false;
	    return true;
	}
	return false;
    }

    public boolean isArray() {
	return m_type == VTYPE.ARRAY;
    }

    private Variant itemshape() {   // Replace a single-types by '1' for later processing
	switch(m_type) {
	case NULL:
	case BOOL:
	case INT:
	case DOUBLE:
	case STRING:
	case LRV:
	case NRV:
	case NRS:
	    return new Variant(1);
	case MASK:
	case INDEX:
	case VECTOR:
	    Variant[] B1 = new Variant[size()];
	    for(int i = 0; i < B1.length; i++) B1[i] = new Variant(1);
	    return new Variant(B1);
	case ARRAY:
	    Variant[] A = (Variant[])m_object;
	    Variant[] B2 = new Variant[A.length];
	    for(int i = 0; i < A.length; i++) B2[i] = A[i].itemshape();
	    return new Variant(B2);
	}
	System.out.println("Variant::itemshape unsupported type: " + m_type);
	return new Variant();
    }

    private Variant runlength() {
	switch(m_type) {
	case INT:
	    Variant[] D = new Variant[1];
	    D[0] = this;
	    return new Variant(D);
	case ARRAY:
	    Variant[] A   = (Variant[])m_object;
	    for(int i = 0; i < A.length; i++) A[i] = A[i].runlength();
	    Variant last  = A[0];
	    int     count = 1;
	    ArrayList<Variant> tuples = new ArrayList<Variant>();
	    for(int i = 1; i < A.length; i++) {
		Variant cur = A[i];
		if(last.equals(cur)) count++;
		else {
		    tuples.add(runlengthGetExpression(count, last));
		    count = 1;
		    last  = cur;
		}
	    }
	    tuples.add(runlengthGetExpression(count, last));
	    if(tuples.size() == 1) return tuples.get(0);
	    Variant[] C = new Variant[tuples.size()];
	    for(int i = 0; i < tuples.size(); i++) C[i] = tuples.get(i);
	    return new Variant(C);
	}
	System.out.println("Variant::runlength unsupported type: " + m_type);
	return new Variant();
    }

    private static Variant runlengthGetExpression(int count, Variant A) {
	if(count == 1) return A;
	Variant[] A1 = (Variant[])A.m_object;
	Variant[] B  = new Variant[A1.length + 1];
	B[0]         = new Variant(count);
	for(int j = 0; j < A1.length; j++) B[j+1]  = A1[j];
	return new Variant(B);
    }

    public Variant shape() {
	switch(m_type) {
	case BOOL:
	case INT:
	case DOUBLE:
	case STRING:
	case LRV:
	case NRV:
	case NRS:
	    return new Variant(new Variant[0]);
	case VECTOR:
	case INDEX:
	    Variant[] A2 = new Variant[1];
	    A2[0] = new Variant(size());
	    return new Variant(A2);
	case ARRAY:
	    return itemshape().sum().runlength();
	case SPMATRIX:
	    Variant[] A1 = new Variant[2];
	    A1[0] = new Variant(((SparseMatrix)m_object).rows());
	    A1[1] = new Variant(((SparseMatrix)m_object).cols());
	    return new Variant(A1);
	}
	System.out.println("Variant::shape unsupported type: " + m_type);
	return new Variant();
    }

    public Variant reshape(Object o) {
	return reshape(new Variant(o));
    }

    public Variant reshape(Variant shape) {
	Variant[] A = getFlatArray();
	Variant   offset = new Variant(0);
	switch(shape.m_type) {
	case INT:
	    Variant[] S2 = new Variant[1];
	    S2[0] = shape;
	    return reshapeHelper(S2, A, offset);
	case ARRAY:
	    Variant[] S = (Variant[])shape.m_object;
	    return reshapeHelper(S, A, offset);
	}
	System.out.println("Variant::reshape["+m_type+"]("+shape.m_type+") unsupported.");
	return new Variant();
    }

    private static Variant reshapeHelper(Variant[] S, Variant[] A, Variant offset) {
	if(S.length == 0) return new Variant(new Variant[0]);
	boolean nonInt = false;
	for(int i = 0; i < S.length; i++) if(S[i].m_type != VTYPE.INT) nonInt = true;
	if(nonInt) {  // assume that elements all ARRAYs or single values we wrap into ARRAY.
	    Variant[] B = new Variant[S.length];
	    for(int i = 0; i < S.length; i++) {
		if(S[i].isArray())
		    B[i] = reshapeHelper((Variant[])S[i].m_object, A, offset);
		else if(S[i].m_type == VTYPE.INT) {
		    Variant[] S2 = new Variant[1];
		    S2[0] = S[i];
		    B[i] = reshapeHelper(S2, A, offset);
		} else {
		    System.out.println("Variant::reshapeHelper S["+i+"]= "+S[i].m_type+" BAD.");
		    return new Variant();
		}
	    }
	    return new Variant(B);
	} 
	// All elements are INTs.
	if(S.length > 1) {
	    Variant[] S2 = new Variant[S.length - 1];
	    for(int i = 1; i < S.length; i++) S2[i-1] = S[i];
	    int size = S[0].toInteger();
	    Variant[] B = new Variant[size];
	    for(int i = 0; i < size; i++) B[i] = reshapeHelper(S2, A, offset);
	    return new Variant(B);
	} else if(S.length == 1) {              // construct dimensional array.
	    int size = S[0].toInteger();
	    Variant[] B = new Variant[size];
	    for(int i = 0; i < size; i++) {
		int j = offset.toInteger();
		B[i] = A[j];
		offset.set(new Variant(j+1));
	    }
	    return new Variant(B);
	}
	System.out.println("reshapeHelper handed a shape of zero elements.");
	return new Variant();
    }

    public int size() {
	switch(m_type) {
	case NULL:   return 0;
	case BOOL:   return 1;
	case INT:    return 1;
	case DOUBLE: return 1;
	case STRING: return ((String)m_object).length();
	case ARRAY:  return ((Variant[])m_object).length;
	case MASK:   return ((Boolean[])m_object).length;
	case INT_MAP:     
	case STRING_MAP:
	    return ((HashMap<?, ?>)m_object).size();
	case INT_SET:
	case STRING_SET:
	    return ((HashSet<?>)m_object).size();
	case INDEX:
	    return ((Integer[])m_object).length;
	case VECTOR:
	    return ((Double[])m_object).length;
	case SPMATRIX:
	    return ((SparseMatrix)m_object).size();
	case LRV:
	case NRV:
	    return 1;
	case NRS:
	    return getRandomSample().size();
	}
	System.out.println("Variant::size unsupported type: " + m_type);
	return -1;
    }

    public int rows() {
	switch(m_type) {
	case ARRAY:
	case MASK:
	    return 1;
	case INDEX:
	    return 1;
	case VECTOR:
	    return size();
	case SPMATRIX:
	    return ((SparseMatrix)m_object).rows();
	}
	System.out.println("Variant::rows unsupported type: " + m_type);
	return -1;
    }

    public int cols() {
	switch(m_type) {
	case ARRAY:
	case MASK:
	    return size();
	case INDEX:
	    return size();
	case VECTOR:
	    return 1;
	case SPMATRIX:
	    return ((SparseMatrix)m_object).cols();
	}
	System.out.println("Variant::cols unsupported type: " + m_type);
	return -1;
    }

    public Variant toIntegers() {
	switch(m_type) {
	case ARRAY:
	    Variant[] A2 = (Variant[])m_object;
	    Variant[] C2 = new Variant[A2.length];
	    for(int i = 0; i < A2.length; i++) C2[i] = new Variant(A2[i].toInteger());
	    return new Variant(C2);
	}
	System.out.println("Variant::toIntegerArray unsupported type: " + m_type);
	return new Variant();
    }

    public Variant toArray() {
	return new Variant(getArray());
    }

    public Variant toArrayDeeply() {
	Variant[] A = getArray();
	Variant[] B = new Variant[A.length];
	for(int i = 0; i < A.length; i++) {
	    switch(A[i].m_type) {
	    case VECTOR:
	    case INDEX:
	    case MASK:
		B[i] = A[i].toArray();
		continue;
	    case ARRAY:
		B[i] = A[i].toArrayDeeply();
		continue;
	    default:
		B[i] = A[i];
		continue;
	    }
	}
	return new Variant(B);
    }

    public Variant[] getArray() {
	switch(m_type) {
	case NULL:
	case BOOL:
	case INT:
	case DOUBLE:
	case STRING:
	    Variant[] B = new Variant[1];
	    B[0] = this;
	    return B;
	case MASK:
	    Boolean[] B2 = (Boolean[])m_object;
	    Variant[] C  = new Variant[B2.length];
	    for(int i = 0; i < B2.length; i++) C[i] = new Variant(B2[i]);
	    return C;
	case ARRAY:
	    return (Variant[])m_object;
	case INT_MAP:     
	case STRING_MAP:
	    Set       K = ((HashMap<?, ?>)m_object).keySet();
	    int       i = 0;
	    int       N = K.size();
	    Variant[] A = new Variant[N];
	    if(m_type == VTYPE.INT_MAP) for(Object key: K) A[i++] = new Variant((Integer)key);
	    else                        for(Object key: K) A[i++] = new Variant((String) key);
	    return A;
	case INT_SET:
	case STRING_SET:
	    Set       K2 = (Set<?>)m_object;
	    int       i2 = 0;
	    int       N2 = K2.size();
	    Variant[] A2 = new Variant[N2];
	    if(m_type == VTYPE.INT_SET) for(Object key: K2) A2[i2++] = new Variant((Integer)key);
	    else                        for(Object key: K2) A2[i2++] = new Variant((String) key);
	    return A2;
	case INDEX:
	    Integer[] I2 = (Integer[])m_object;
	    Variant[] V3 = new Variant[I2.length];
	    for(int j = 0; j < I2.length; j++) V3[j] = new Variant(I2[j]);
	    return V3;
	case VECTOR:
	    Double [] D2 = (Double[])m_object;
	    Variant[] V2 = new Variant[D2.length];
	    for(int j = 0; j < D2.length; j++) V2[j] = new Variant(D2[j]);
	    return V2;
	}
	System.out.println("Variant::getArray unsupported type: " + m_type);
	return new Variant[0];
    }

    public Variant[] getFlatArray() {
	Variant[] A = getArray();
	int length = 0;
	for(int i = 0; i < A.length; i++) length += A[i].size();
	Variant[] B = new Variant[length];
	int k = 0;
	for(int i = 0; i < A.length; i++) {
	    if(A[i].isArray()) {
		Variant[] C = A[i].getFlatArray();
		for(int j = 0; j < C.length; j++) B[k++] = C[j];
	    } else 
		B[k++] = A[i];
	}
	return B;
    }

    public Integer[] getIndex() {
	switch(m_type) {
	case INT:
	    int N = toInteger();
	    Integer[] D2 = new Integer[N];
	    for(int i = 0; i < N; i++) D2[i] = 0;
	    return D2;
	case STRING:
	    Integer[] D = new Integer[1];
	    D[0] = toInteger();
	    return D;
	case MASK:
	    Boolean[] B2 = (Boolean[])m_object;
	    int       k  = 0;
	    Integer[] C  = new Integer[sum().toInteger()];
	    for(int i = 0; i < B2.length; i++) if(B2[i]) C[k++] = i;
	    return C;
	case INDEX:
	    return (Integer[])m_object;
	case VECTOR:
	    Double [] V2 = (Double[])m_object;
	    Integer[] I2 = new Integer[V2.length];
	    for(int i = 0; i < V2.length; i++) I2[i] = V2[i].intValue();
	    return I2;
	case ARRAY: // can handle 2D array flatten so far.
	    Variant[] V  = (Variant[])m_object;
	    int length = 0;
	    for(int i = 0; i < V.length; i++) length += V[i].size();
	    Integer[] C2 = new Integer[length];
	    int k2 = 0;
	    for(int i = 0; i < V.length; i++)
		if(V[i].size() == 1) 
		    C2[k2++] = V[i].toInteger();
		else {
		    Integer[] T = V[i].getIndex();
		    for(int j = 0; j < T.length; j++) 
			C2[k2++] = T[j];
		}
	    return C2;
	case INT_MAP:     
	case STRING_MAP:
	    Set       K  = ((HashMap<?, ?>)m_object).keySet();
	    int       k3 = 0;
	    Integer[] C3 = new Integer[K.size()];
	    if(m_type == VTYPE.INT_MAP) for(Object key: K) C3[k3++] = (Integer)key;
	    else                        for(Object key: K) C3[k3++] = 0; //TODO: (String)key
	    return C3;
	}
	System.out.println("Variant::getIndex["+m_type+"] unsupported.");
	return new Integer[0];
    }

    public Variant toIndex() {
	return new Variant(getIndex());
    }

    public Variant toVector() {
	return createVector(this);
    }

    public Double[] getVector() {
	switch(m_type) {
	case INT:
	    int N = toInteger();
	    Double[] D2 = new Double[N];
	    for(int i = 0; i < N; i++) D2[i] = 0.0;
	    return D2;
	case DOUBLE:
	case STRING:
	    Double[] D = new Double[1];
	    D[0] = toDouble();
	    return D;
	case MASK:
	    Boolean[] B2 = (Boolean[])m_object;
	    Double [] C  = new Double[B2.length];
	    for(int i = 0; i < B2.length; i++) C[i] = B2[i] ? 1.0 : 0.0;
	    return C;
	case ARRAY: // can handle 2D array flatten so far.
	    Variant[] V  = (Variant[])m_object;
	    int length = 0;
	    for(int i = 0; i < V.length; i++) length += V[i].size();
	    Double [] C2 = new Double[length];
	    int k = 0;
	    for(int i = 0; i < V.length; i++)
		if(V[i].size() == 1) 
		    C2[k++] = V[i].toDouble();
		else {
		    Double[] T = V[i].getVector();
		    for(int j = 0; j < T.length; j++) 
			C2[k++] = T[j];
		}
	    return C2;
	case INT_MAP:     
	case STRING_MAP:
	    Set       K = ((HashMap<?, ?>)m_object).keySet();
	    int       i = 0;
	    Double[] C3 = new Double[K.size()];
	    if(m_type == VTYPE.INT_MAP) for(Object key: K) C3[i++] = new Double((Integer)key);
	    else                        for(Object key: K) C3[i++] = 0.0; //TODO: (String)key
	    return C3;
	case INDEX:
	    Integer[] I2 = (Integer[])m_object;
	    Double [] V3 = new Double[I2.length]; 
	    for(int j = 0; j < I2.length; j++) V3[j] = new Double(I2[j]);
	    return V3;
	case VECTOR:
	    return (Double[])m_object;
	}
	System.out.println("Variant::getVector unsupported type: " + m_type);
	return new Double[0];
    }

    public double[] getDoubleArray() {
	Double[] D = getVector();
	double[] d = new double[D.length];
	for(int i = 0; i < D.length; i++) d[i] = D[i].doubleValue();
	return d;
    }

    public double[][] getMatrix() {
	switch(m_type) {
	case LRV:
	case NRV:
	case NRS:
	    return getNumericRandomVariable().getData();
	}
	System.out.println("Variant::getMatrix unsupported type: " + m_type);
	return new double[2][0];
    }

    public double[][] getCumulativeData() {
	switch(m_type) {
	case LRV:
	case NRV:
	case NRS:
	    return getNumericRandomVariable().getCumulativeData();
	}
	System.out.println("Variant::getCumulativeData unsupported type: " + m_type);
	return new double[2][0];
    }

    public SparseMatrix getSparseMatrix() {
	switch(m_type) {
	case ARRAY:
	    SparseMatrix S = new SparseMatrix();
	    Variant[] A = (Variant[])m_object;
	    for(int r = 0; r < A.length; r++) {
		if(A[r].m_type != VTYPE.ARRAY) {
		    System.out.println("getSparseMatrix["+m_type+"]["+r+"] = "+A[r].m_type+
				       " must be an array.");
		    return new SparseMatrix();
		}
		Variant[] B = (Variant[])A[r].m_object;
		for(int c = 0; c < B.length; c++) 
		    S.set(r,c,(Double)(B[c].getDouble().m_object));
	    }
	    return S;
	case SPMATRIX:
	    return (SparseMatrix)m_object;
	}
	System.out.println("Variant::getSparseMatrix["+m_type+"] unsupported.");
	return new SparseMatrix();
    }

    public LazyRandomVariable getLazyRandomVariable() {
	switch(m_type) {
	case LRV:
	    return (LazyRandomVariable)m_object;
	}
	System.out.println("Variant::getLazyRandomVariable["+m_type+"] unsupported.");
	return null;
    }

    public NumericRandomVariable getNumericRandomVariable() {
	switch(m_type) {
	case LRV:
	    return getLazyRandomVariable().get().get();
	case NRV:
	    return (NumericRandomVariable)m_object;
	case NRS:
	    return ((RandomSample)m_object).getContinuousRandomVariable(); 
	}
	System.out.println("Variant::getNumericRandomVariable["+m_type+"] unsupported.");
	return null;
    }

    public RandomSample getRandomSample() {
	switch(m_type) {
	case NRS:
	    return (RandomSample)m_object;
	}
	System.out.println("Variant::getRandomSample["+m_type+"] unsupported.");
	return null;
    }

    public Variant keys() {
	switch(m_type) {
	case ARRAY:
	    HashSet<String>  sset2 = new HashSet<String> ();
	    HashSet<Integer> iset2 = new HashSet<Integer>();
	    Variant v2 = new Variant();
	    v2.m_object = sset2;              // weak default.
	    v2.m_type   = VTYPE.STRING_SET;
	    Boolean isInt = null;         // Either pure INTs or pure STRINGs else fail. unset.
	    Variant[] array  = (Variant[])m_object;
	    for(Variant v: array)  {
		if(isInt == null) 
		    if(     v.m_type == VTYPE.INT)    isInt = true;
		    else if(v.m_type == VTYPE.STRING) isInt = false;
		    else {
			System.out.println("keys() not implemented for: " + v.m_type);
			return v2;
		    }
		if(isInt)
		    if(v.m_type == VTYPE.INT) iset2.add((Integer)v.m_object);
		    else {
			System.out.println("keys is Set<Integer>, not: " + v.m_type);
			return v2;
		    }
		else
		    if(v.m_type == VTYPE.STRING) sset2.add((String)v.m_object);
		    else {
			System.out.println("keys() is Set<String>, not: " + v.m_type);
			return v2;
		    }
	    }
	    if(isInt == null) return v2;
	    if(isInt) { v2.m_object = iset2; v2.m_type = VTYPE.INT_SET;    } 
	    else      { v2.m_object = sset2; v2.m_type = VTYPE.STRING_SET; }
	    return v2;
	case INT_SET:
	case STRING_SET:
	    return copy();
	case INT_MAP:
	    HashSet<Integer> iset = new HashSet<Integer>();
	    Set<Integer>  ikeyset = ((HashMap<Integer, Variant>)m_object).keySet();
	    for(Integer i: ikeyset) iset.add(i);
	    return new Variant(iset, VTYPE.INT_SET);
	case STRING_MAP:
	    HashSet<String> sset = new HashSet<String>();
	    Set<String>  skeyset = ((HashMap<String, Variant>)m_object).keySet();
	    for(String s: skeyset) sset.add(s);
	    return new Variant(sset, VTYPE.STRING_SET);
	case FUNCTION:
	    // not yet supported.
	}
	System.out.println("Variant::keys unsupported type: " + m_type);
	return new Variant();
    }

    public Variant toSet() {
	switch(m_type) {
	case INT:
	case DOUBLE:
	    Variant int_set = createIntSet();
	    int_set.append(toInteger());
	    return int_set;
	case STRING:
	    Variant str_set = createStringSet();
	    str_set.append(toString());
	    return str_set;
	case ARRAY:
	case INT_SET:
	case STRING_SET:
	case INT_MAP:
	case STRING_MAP:
	case FUNCTION:
	    return keys();
	}
	System.out.println("Variant::toSet unsupported type: " + m_type);
	return new Variant();
    }

    public static int getIndex(int i, int N) {
	// deal with modulo issues
	if(i > 0) return i % N;
	return ((1-i/N)*N + i) % N;
    }

    public Variant get(int i) {
	switch(m_type) {
	case BOOL:
	case INT:
	case DOUBLE:
	case STRING:
	    return this; // allows values to be viewed as arbitrarily long constant array
	case ARRAY:   
	    return ((Variant[])m_object)[getIndex(i, size())];               // return original
	case MASK:
	    return new Variant(((Boolean[])m_object)[getIndex(i, size())]);  // return a copy
	case INT_MAP: 
	    return ((HashMap<Integer, Variant>)m_object).get(i);
	case INDEX:
	    return new Variant(((Integer[])m_object)[getIndex(i, size())]);
	case VECTOR:
	    return new Variant(((Double[])m_object)[getIndex(i, size())]);
	case SPMATRIX:
	    return createVector(((SparseMatrix)m_object).get(i));
	case NRS:
	    return new Variant(getRandomSample().get(i));
	}
	System.out.println("Variant::get(index) unsupported type: " + m_type);
	return new Variant();
    }

    public Variant get(String key) {
	return get(new Variant(key));
    }

    public Variant get(Variant key) {
	return get(key, null);
    }

    public Variant get(Variant key, Variant dflt) {
	switch(m_type) {
	case MASK:
	    switch(key.m_type) {
	    case INT:     return get(key.toInteger());
	    case STRING:  return applySlice(createSlice(key));
	    case SLICE:   return applySlice(key);
	    case MASK:    return applyMask(key);
	    case INDEX:   return applyIndex(key);
	    case ARRAY:   return applyArray(key);
	    case INT_SET: return applyIntSet(key);
	    }
	    break;
	case INDEX:
	    switch(key.m_type) {
	    case INT:     return get(key.toInteger());
	    case STRING:  return applySlice(createSlice(key));
	    case SLICE:   return applySlice(key);
	    case MASK:    return applyMask(key);
	    case ARRAY:   return applyArray(key);
	    case INT_SET: return applyIntSet(key);
	    case INDEX:   return applyIndex(key);
	    }
	    break;
	case VECTOR:
	    switch(key.m_type) {
	    case INT:	  return get(key.toInteger());
	    case STRING:  return applySlice(createSlice(key));
	    case SLICE:   return applySlice(key);
	    case MASK:	  return applyMask(key);
	    case INDEX:   return applyIndex(key);
	    case ARRAY: //return applyArray(key);   // applyArray to depend on what inside?
		Variant[] K2 = (Variant[])key.m_object;
		Variant[] B2 = new Variant[K2.length];
		for(int i = 0; i < K2.length; i++) B2[i] = get(K2[i]);
		return new Variant(B2);
	    }
	    break;
	case ARRAY:
	    switch(key.m_type) {
	    case INT:     return get(key.toInteger());
	    case STRING:  return applySlice(createSlice(key));
	    case SLICE:   return applySlice(key);
	    case MASK:    return applyMask(key);
	    case INDEX:   return applyIndex(key);
	    case ARRAY:   return applyArray(key);
	    case INT_SET: return applyIntSet(key);
	    case STRING_SET: // push down a level (defer)
		Variant[] A = (Variant[])m_object;
		Variant[] B = new Variant[A.length];
		for(int i = 0; i < B.length; i++) B[i] = A[i].get(key, dflt);
		return new Variant(B);
	    }
	    break;
	case STRING_MAP:
	    HashMap<String, Variant> str_map = (HashMap<String, Variant>)m_object;
	    switch(key.m_type) {
	    case STRING: 
		String skey = (String)key.m_object;
		if(str_map.containsKey(skey)) return str_map.get(skey);
		return dflt;
	    case INT:
	    case MASK:   // push mask down a level (defer)
	    case SLICE:  // push slice down a level (defer)
		HashMap<String, Variant> smap2 = new HashMap<String, Variant>();
		for(String s: str_map.keySet()) smap2.put(s, str_map.get(s).get(key));
		return new Variant(smap2, VTYPE.STRING_MAP);
	    case ARRAY: 
		// Return an array of gets, one for each key in given key array.
		Variant[] V = new Variant[key.size()];
		Variant[] keys = key.getArray();
		if(dflt != null) {
		    Variant[] dflts = dflt.getArray();
		    for(int i = 0; i < key.size(); i++) V[i] = get(keys[i], dflts[i]);
		} else 
		    for(int i = 0; i < key.size(); i++) V[i] = get(keys[i], null);
		return new Variant(V);
	    case STRING_SET:
		HashMap<String, Variant> smap3 = new HashMap<String, Variant>();
		for(String k: key.toStringSet()) 
		    if(str_map.containsKey(k)) smap3.put(k, str_map.get(k));
		return new Variant(smap3, VTYPE.STRING_MAP);
	    case STRING_MAP:
		return createMaskV(this, key); 
	    }
	    break;
	case INT_MAP:
	    HashMap<Integer, Variant> int_map = (HashMap<Integer, Variant>)m_object;
	    switch(key.m_type) {
	    case INT:
		Integer ikey = (Integer)key.m_object;
		if(int_map.containsKey(ikey)) return int_map.get(ikey);
		return dflt;
	    case DOUBLE:
		double dkey = ((Double)key.m_object).doubleValue();
		if(dkey == (double)(int)dkey) {
		    Integer iikey = new Integer((int)dkey);
		    if(int_map.containsKey(iikey)) return int_map.get(iikey);
		}
		return dflt;
	    case MASK:   // push mask  down a level (defer)
	    case SLICE:  // push slice down a level (defer)
		HashMap<Integer, Variant> imap = new HashMap<Integer, Variant>();
		for(Integer i: int_map.keySet()) imap.put(i, int_map.get(i).get(key));
		return new Variant(imap, VTYPE.INT_MAP);
	    case ARRAY:
		// Return an array of gets, one for each key in given key array.
		Variant[] V    = new Variant[key.size()];
		Variant[] keys = key.getArray();
		if(dflt != null) {
		    Variant[] dflts = dflt.getArray();
		    for(int i = 0; i < key.size(); i++) V[i] = get(keys[i], dflts[i]);
		} else 
		    for(int i = 0; i < key.size(); i++) V[i] = get(keys[i], null);
		return new Variant(V);
	    case INT_SET:
		HashMap<Integer, Variant> imap3 = new HashMap<Integer, Variant>();
		for(Integer k: key.toIntSet()) 
		    if(int_map.containsKey(k)) imap3.put(k, int_map.get(k));
		return new Variant(imap3, VTYPE.INT_MAP);
	    }
	    break;
	case SPMATRIX:
	    switch(key.m_type) {
	    case INT:
		return get((int)(Integer)key.m_object);
	    }
	    break;
	case FUNCTION:
	    return ((Function)m_object).get(key);
	case NRS:
	    switch(key.m_type) {
	    case INT:     return get(key.toInteger());
	    case MASK:    return applyMask(key);
	    case SLICE:   return applySlice(key);
	    }
	    break;
	}
	System.out.println("Variant::get["+m_type+"]("+key.m_type+") unsupported.");
	return new Variant();
    }

    private Variant applyIndex(Variant key) {
	return applyIndex(key, null);
    }

    private Variant applyIndex(Variant key, Variant value) {
	switch(m_type) {
	case MASK:
	case VECTOR:
	case ARRAY:
	    switch(key.m_type) {
	    case INDEX:
		Integer[] K = (Integer[])key.m_object;
		if(value == null) {
		    Variant[] V = new Variant[K.length];
		    for(int i = 0; i < K.length; i++) V[i] = get(K[i]);
		    return new Variant(V);
		} else {
		    for(int i = 0; i < K.length; i++) set(K[i], value.get(i));
		    return this;
		}
	    }
	case INDEX:
	    switch(key.m_type) {
	    case INDEX:
		Integer[] K2 = (Integer[])key.m_object;
		if(value == null) {
		    Integer[] V2 = new Integer[K2.length];
		    for(int i = 0; i < K2.length; i++) V2[i] = (Integer)(get(K2[i]).m_object);
		    return new Variant(V2);
		} else {
		    for(int i = 0; i < K2.length; i++) set(K2[i], value.get(i));
		    return this;
		}
	    }	    
	}
	System.out.println("Variant::applyIndex["+m_type+"]("+key.m_type+") unsupported.");
	return new Variant();
    }

    private Variant applyArray(Variant key) {
	return applyArray(key, null);
    }

    private Variant applyArray(Variant key, Variant value) {
	switch(m_type) {
	case ARRAY:
	case MASK:
	case INDEX:
	case VECTOR:
	    switch(key.m_type) {
	    case ARRAY:
		if(value == null) {
		    Variant[] keys = (Variant[])key.m_object;
		    Variant[] V    = new Variant[key.size()];
		    for(int i = 0; i < key.size(); i++) 
			if(keys[i].m_type == VTYPE.INT) V[i] = get(keys[i].toInteger());
			else                            V[i] = get(i).get(keys[i]);
		    return new Variant(V);
		} else {
		    Variant[] keys = (Variant[])key.m_object;
		    for(int i = 0; i < key.size(); i++)
			if(keys[i].m_type == VTYPE.INT) set(keys[i].toInteger(), value.get(i));
			else                            get(i).set(keys[i],      value.get(i));
		    return this;
		}
	    }
	}
	System.out.println("Variant::applyArray unsupported types: "+m_type+", "+key.m_type);
	return new Variant();
    }

    private Variant applyMask(Variant key) {
	return applyMask(key, null);
    }

    private Variant applyMask(Variant key, Variant value) {
	switch(m_type) {
	case INDEX:          // <<<< TODO: Need to break this out 
	case VECTOR:         // <<<< TODO: Need to break this out
	case ARRAY:
	    switch(key.m_type) {
	    case MASK:
		if(size() == key.size()) {
		    Boolean[] keys = (Boolean[])key.m_object;
		    if(value == null) {  // get
			int N = 0;
			for(int i = 0; i < keys.length; i++) if(keys[i]) N++;
			Variant[] V = new Variant[N];
			int j = 0;
			for(int i = 0; i < size(); i++) if(keys[i]) V[j++] = get(i);
			return new Variant(V);
		    } else {             // set
			int j = 0;
			for(int i = 0; i < size(); i++) if(keys[i]) set(i, value.get(j++));
			return this;
		    }
		} else {
		    System.out.println("Variant::applyMask["+m_type+"]("+key.m_type+
				       ") size mismatch: "+size()+", "+key.size());
		    return new Variant();
		}
	    }
	    break;
	case MASK:
	    switch(key.m_type) {
	    case MASK:
		if(size() == key.size()) {
		    Boolean[] keys = (Boolean[])key.m_object;
		    if(value == null) {  // get
			int N = 0;
			for(int i = 0; i < keys.length; i++) if(keys[i]) N++;
			Boolean[] A = (Boolean[])m_object;
			Boolean[] B = new Boolean[N];
			int j = 0;
			for(int i = 0; i < size(); i++) if(keys[i]) B[j++] = A[i];
			return new Variant(B);
		    } else {             // set
			int j = 0;
			for(int i = 0; i < size(); i++) if(keys[i]) set(i, value.get(j++));
			return this;
		    }
		} else {
		    System.out.println("Variant::applyMask["+m_type+"]("+key.m_type+
				       ") size mismatch: "+size()+", "+key.size());
		    return new Variant();
		}
	    }
	    break;
	case NRS:
	    switch(key.m_type) {
	    case MASK:
		if(size() == key.size()) {
		    Boolean[] keys = (Boolean[])key.m_object;
		    if(value == null) return new Variant(getRandomSample().get(keys));
		    else {             // set
			System.out.println("applyMask["+m_type+"]("+key.m_type+
					   ") set() not supported yet.");
		    }
		} else {
		    System.out.println("Variant::applyMask["+m_type+"]("+key.m_type+
				       ") size mismatch: "+size()+", "+key.size());
		    return new Variant();
		}
	    }
	    break;
	}
	System.out.println("Variant::applyMask["+m_type+"]("+key.m_type+") unsupported.");
	return new Variant();
    }

    private Variant applySlice(Variant key) {
	return applySlice(key, null);
    }

    private Variant applySlice(Variant key, Variant value) {
	switch(m_type) {
	case ARRAY:
	case MASK:
	case INDEX:
	case VECTOR:
	case NRS:
	    switch(key.m_type) {
	    case SLICE:
		int N = size();

		int start = (int)((Integer[])key.m_object)[0];
		int stop  = (int)((Integer[])key.m_object)[1];
		int step  = (int)((Integer[])key.m_object)[2];

		while(stop <= start) stop += N;

		int begin = start;
		int dir   = 1;
		if(step < 0) {dir = -1; step = -step; begin = stop-1;}
		int M = Math.max((stop - start)/step, (stop-start+(step-1))/step);

		return applySliceHelper(key, value, M, begin, dir * step);
	    }
	}
	System.out.println("Variant::applySlice unsupported types: "+m_type+", "+key.m_type);
	return new Variant();
    }

    private Variant applySliceHelper(Variant key, Variant value, int M, int start, int step) {
	int N = size();
	if(value == null) {  // get
	    switch(m_type) {
	    case ARRAY:
		Variant[] V  = (Variant[])m_object;
		Variant[] V2 = new Variant[M];
		for(int i = 0; i < M; i++) V2[i] = V[getIndex(start + i*step, N)];
		return new Variant(V2);
	    case MASK:
		Boolean[] B  = (Boolean[])m_object;
		Boolean[] B2 = new Boolean[M];
		for(int i = 0; i < M; i++) B2[i] = B[getIndex(start + i*step, N)];
		return new Variant(B2);
	    case INDEX:
		Integer[] I  = (Integer[])m_object;
		Integer[] I2 = new Integer[M];
		for(int i = 0; i < M; i++) I2[i] = I[getIndex(start + i*step, N)];
		return new Variant(I2);
	    case VECTOR:
		Double[] D  = (Double[])m_object;
		Double[] D2 = new Double[M];
		for(int i = 0; i < M; i++) D2[i] = D[getIndex(start + i*step, N)];
		return createVector(D2);
	    case NRS:
		return new Variant(getRandomSample().get(M, start, step));
	    }
	} else {        // set
	    int j = 0;  // index into array value, if applicable.
	    switch(m_type) {
	    case ARRAY:
		Variant[] V  = (Variant[])m_object;
		for(int i = 0; i < M; i++) V[getIndex(start+i*step,N)] = value.get(i);
		return this;
	    case MASK:
		Boolean[] B = (Boolean[])m_object;
		for(int i = 0; i < M; i++) B[getIndex(start+i*step,N)] = value.get(j++).isTrue();
		return this;
	    case INDEX:
		Integer[] I  = (Integer[])m_object;
		for(int i = 0; i < M; i++) I[getIndex(start+i*step,N)]= value.get(i).toInteger();
		return this;
	    case VECTOR:
		Double[] D  = (Double[])m_object;
		for(int i = 0; i < M; i++) D[getIndex(start+i*step,N)] = value.get(i).toDouble();
		return this;
	    }
	}

	System.out.println("Variant::applySliceHelp unsupported types: "+m_type+", "+key.m_type);
	return new Variant();
    }

    private Variant applyIntSet(Variant key) {
	return applyIntSet(key, null);
    }

    private Variant applyIntSet(Variant key, Variant value) {
	switch(m_type) {
	case ARRAY:
	case MASK:
	case INDEX:
	case VECTOR:
	    switch(key.m_type) {
	    case INT_SET:
		if(value == null) {
		    int j = 0;
		    Variant[] V = new Variant[key.size()];
		    Set<Integer> keyset = key.toIntSet();
		    for(int i = 0; i < size(); i++) if(keyset.contains(i)) V[j++] = get(i);
		    return new Variant(V);
		} else {
		    for(Integer i: key.toIntSet()) set((int)i, value.get(i));  // parallel
		    return this;
		}
	    }
	}
	System.out.println("Variant::applyIntSet unsupported types: "+m_type+", "+key.m_type);
	return new Variant();
    }

    private Variant set(Variant value) {
	m_object = value.m_object;
	m_type   = value.m_type;
	return this;
    }

    private Variant set(int i, Variant value) {
	switch(m_type) {
	case ARRAY:   
	    ((Variant[])m_object)[getIndex(i, size())].set(value);
	    return this;
	case MASK:
	    ((Boolean[])m_object)[getIndex(i, size())] = value.isTrue(); 
	    return this;
	case INDEX:
	    ((Integer[])m_object)[getIndex(i, size())] = value.toInteger(); 
	    return this;
	case VECTOR:
	    ((Double[])m_object)[getIndex(i, size())] = value.toDouble(); 
	    return this;
	case INT_MAP: 
	    ((HashMap<Integer, Variant>)m_object).put(i, value);
	    return this;
	case SPMATRIX:
	    ((SparseMatrix)m_object).set(i, value.toDouble());
	    return this;
	}
	System.out.println("Variant::set["+m_type+"](INT) = "+value.m_type+"unsupported.");
	return new Variant();
    }

    public Variant set(Variant key, Variant value) {
	switch(m_type) {
	case MASK:
	case ARRAY:
	case INDEX:
	case VECTOR:
	    switch(key.m_type) {
	    case INT:     return set(key.toInteger(), value);
	    case MASK:    return applyMask(key, value);
	    case SLICE:   return applySlice(key, value);
	    case STRING:  return applySlice(createSlice(key), value);
	    case INDEX:   return applyIndex(key, value);
	    case ARRAY:   return applyArray(key, value);
	    case INT_SET: return applyIntSet(key, value);
	    }
	    break;
	case STRING_MAP:
	    HashMap<String, Variant> str_map = (HashMap<String, Variant>)m_object;
	    switch(key.m_type) {
	    case STRING: 
		str_map.put((String)key.m_object, value);
		return this;
	    case MASK:    // push down a level (defer)
	    case SLICE:   // push down a level (defer)
	    case INT_SET: // push down a level (defer)
		for(String s: str_map.keySet()) str_map.get(s).set(key, value);
		return this;
	    case ARRAY: 
		if(key.isIntegers()) { // push down a level (defer)
		    for(String s: str_map.keySet()) str_map.get(s).set(key, value);
		    return this;
		}
	    break;
	    }
	case SPMATRIX:
	    switch(key.m_type) {
	    case INT:     return set(key.toInteger(), value);
	    }
	    break;
	}
	System.out.println("Variant::set["+m_type+"]("+key.m_type+") = "+value.m_type+
			   " is unsupported.");
	return new Variant();
    }

    public void setTrue() {
	switch(m_type) {
	case BOOL:
	case INT:
	case DOUBLE:
	    m_object = true;
	    return;
	case STRING:
	    m_object = "T";
	    return;
	case MASK:
	    Boolean[] b = (Boolean[])m_object;
	    for(int i = 0; i < size(); i++) b[i] = true;
	    return;
	case ARRAY:
	    Variant[] array = (Variant[])m_object;
	    for(int i = 0; i < size(); i++) array[i].setTrue();
	    return;
	}
	System.out.println("Variant::setTrue() unsupported type: " + m_type);
    }

    public Variant setFalse() {
	switch(m_type) {
	case BOOL:
	case INT:
	case DOUBLE:
	    m_object = false;
	    return this;
	case STRING:
	    m_object = "False";
	    return this;
	case MASK:
	    Boolean[] b = (Boolean[])m_object;
	    for(int i = 0; i < size(); i++) b[i] = false;
	    return this;
	case ARRAY:
	    Variant[] array = (Variant[])m_object;
	    for(int i = 0; i < size(); i++) array[i].setFalse();
	    return this;
	case STRING_MAP: // push-down
	    HashMap<String, Variant> A2 = (HashMap<String, Variant>)m_object;
	    for(String key: A2.keySet()) A2.get(key).setFalse();
	    return this;
	}
	System.out.println("Variant::setFalse() unsupported type: " + m_type);
	return new Variant();
    }

    public void setZero() {
	switch(m_type) {
	case BOOL:
	    m_object = new Boolean(false);
	    return;
	case INT:
	    m_object = new Integer(0);
	    return;
	case DOUBLE:
	    m_object = new Double(0.0);
	    return;
	case STRING:
	    m_object = "";
	    return;
	case ARRAY:
	    Variant[] A = (Variant[])m_object;
	    for(int i = 0; i < A.length; i++) A[i].setZero();
	    return;
	case STRING_MAP:
	    HashMap<String, Variant> S = (HashMap<String, Variant>)m_object;
	    for(String s: S.keySet()) S.get(s).setZero();
	    return;
	}
	System.out.println("Variant::setZero() unsupported type: " + m_type);
    }

    public Variant setInfinite() {
	switch(m_type) {
	case BOOL:
	case INT:
	case DOUBLE:
	case STRING:
	    m_object = Double.POSITIVE_INFINITY;
	    m_type   = VTYPE.DOUBLE;
	    return this;
	case ARRAY:
	    Variant[] A = (Variant[])m_object;
	    for(int i = 0; i < A.length; i++) A[i].setInfinite();
	    return this;
	case STRING_MAP:
	    HashMap<String, Variant> S = (HashMap<String, Variant>)m_object;
	    for(String s: S.keySet()) S.get(s).setInfinite();
	    return this;
	}
	System.out.println("Variant::setInfinite() unsupported type: " + m_type);
	return new Variant();
    }

    public Variant pop() {
	switch(m_type) {
	case ARRAY:
	    Variant[] A2 = (Variant[])m_object;
	    if(A2.length == 0) return this;             // no effect. Idempotent.
	    Variant[] B2 = new Variant[A2.length - 1];
	    for(int i = 1; i < A2.length; i++) B2[i-1] = A2[i];
	    Variant B = new Variant(B2);
	    m_object = B.m_object;
	    return A2[0];
	}
	System.out.println("Variant::pop["+m_type+"] unsupported.");
	return new Variant();
    }

    public Variant lookup(Object key) {
	return lookup(new Variant(key));
    }

    public Variant lookup(Variant key) {
	// Sequential lookup versus usual [multi-]lookup. Recombining tree effect.
	switch(m_type) {
	case INT_MAP:
	case STRING_MAP:
	    switch(key.m_type) {
	    case NULL:
	    case BOOL:
	    case INT:
	    case DOUBLE:
	    case STRING:
	    case INT_MAP:
	    case STRING_MAP:
		return get(key);                   // a pass-through
	    case ARRAY:                            // sequential lookup.
		Variant key2 = key.copy();
		Variant K    = key2.pop();
		Variant L    = get(K);
		if(key2.size() == 0) return L;     // terminal
		return L.lookup(key2);             // let it ride
	    }
	    break;
	case FUNCTION:
	    switch(key.m_type) {
	    case INT:
	    case DOUBLE:
	    case STRING:
		return get(key);
	    case ARRAY:
		Variant[] K2 = (Variant[])key.m_object;  // key ARRAY's are structural here.
		if(K2.length == 1) return get(K2[0]);    // terminal
		Variant[] F2 = new Variant[K2.length];
		for(int i = 0; i < K2.length; i++) F2[i] = get(K2[i]);
		return new Variant(F2);
	    }
	    break;
	case ARRAY:
	    Variant[] A2 = (Variant[])m_object;                 // ARRAY's are structural.
	    Variant[] C2 = new Variant[A2.length];
	    switch(key.m_type) {
	    case INT:
	    case DOUBLE:
	    case STRING:
		for(int i = 0; i < A2.length; i++) C2[i] = A2[i].get(key);
		return new Variant(C2);
	    case ARRAY:
		Variant key2 = key.copy();
		Variant K    = key2.pop();
		switch(K.m_type) {
		case INT:
		case DOUBLE:
		case STRING:
		    Variant L = lookup(K);
		    if(key2.size() == 0) return L;     // terminal
		    return L.lookup(key2);          // let it ride
		case ARRAY:
		    for(int i = 0; i < A2.length; i++) C2[i] = A2[i].lookup(K.get(i));
		    Variant R = new Variant(C2);
		    if(key2.size() == 0) return R;     // terminal
		    return R.lookup(key2);          // let it ride
		}
	    }
	    break;
	}
	System.out.println("Variant::lookup["+m_type+"]("+key.m_type+") unsupported.");
	return new Variant();
    }

    public Variant copy() {
	switch(m_type) {
	case NULL:       return new Variant();
	case BOOL:       return new Variant((Boolean)m_object);
	case INT:        return new Variant((Integer)m_object);
	case DOUBLE:     return new Variant((Double) m_object);
	case STRING:     return new Variant((String) m_object);
	case MASK:
	    Boolean[] b  = (Boolean[])m_object;
	    Boolean[] b2 = new Boolean[b.length];;
	    for(int i = 0; i < b.length; i++) b2[i] = b[i];
	    return new Variant(b2);
	case INDEX:
	    Integer[] A2 = (Integer[])m_object;
	    Integer[] C2 = new Integer[A2.length];
	    for(int i = 0; i < A2.length; i++) C2[i] = A2[i];
	    return new Variant(C2);
	case VECTOR:
	    Double[] A3 = (Double[])m_object;
	    Double[] C3 = new Double[A3.length];
	    for(int i = 0; i < A3.length; i++) C3[i] = A3[i];
	    return createVector(C3);
	case ARRAY:
	    Variant[] array = new Variant[size()];
	    for(int i = 0; i < size(); i++) array[i] = get(i).copy();
	    return new Variant(array);
	case INT_SET:
	    HashSet<Integer> int_set = new HashSet<Integer>();
	    HashSet<?>          iset = (HashSet<?>)m_object;
	    for(Object key: iset) int_set.add((Integer)key);
	    return new Variant(int_set, VTYPE.INT_SET);
	case STRING_SET:
	    HashSet<String> string_set = new HashSet<String>();
	    HashSet<?>            sset = (HashSet<?>)m_object;
	    for(Object key: sset) string_set.add((String)key);
	    return new Variant(string_set, VTYPE.STRING_SET);
	case INT_MAP:
	    HashMap<Integer,  Variant> int_map = new HashMap<Integer,  Variant>();
	    HashMap<?,?>                  imap = (HashMap<?, ?>)m_object;
	    for(Object key: imap.keySet()) 
		int_map.put((Integer)key, ((Variant)imap.get(key)).copy());
	    return new Variant(int_map, VTYPE.INT_MAP);
	case STRING_MAP:
	    HashMap<String, Variant> str_map = new HashMap<String,  Variant>();
	    HashMap<?,?>                smap = (HashMap<?, ?>)m_object;
	    for(Object key: smap.keySet()) 
		str_map.put((String)key, ((Variant)smap.get(key)).copy());
	    return new Variant(str_map, VTYPE.STRING_MAP);
	case SLICE:
	    Integer[] ii  = (Integer[])m_object;
	    Integer[] ii2 = new Integer[ii.length];
	    for(int i = 0; i < ii.length; i++) ii2[i] = ii[i];
	    Variant v = new Variant();
	    v.m_object = ii2;
	    v.m_type   = VTYPE.SLICE;
	    return v;
	case SPMATRIX:
	    return new Variant(((SparseMatrix)m_object).copy());
	case LRV:
	    return new Variant(((LazyRandomVariable)m_object).copy());
	case FUNCTION:
	    return ((Function)m_object).copy();
	}
	System.out.println("Variant::copy unsupported for type: " + m_type);
	return new Variant();
    }

    public Variant transpose() {
	switch(m_type) {
	case INDEX:
	case VECTOR:
	case MASK:
	    return this;
	case ARRAY:
	    int Ncols = size();                   // formerly 'rows'
	    if(Ncols == 0) return new Variant();  // bail out
	    int Nrows = get(0).size();            // formerly 'cols'
	    if(Nrows == 0) return new Variant();  // bail out
	    Variant[] C = new Variant[Nrows];     // top level = rows
	    for(int i = 0; i < Nrows; i++) {
		Variant[] row = new Variant[Ncols];
		for(int j = 0; j < Ncols; j++) row[j] = get(j).get(i);
		C[i] = new Variant(row);
	    }
	    return new Variant(C);
	case SPMATRIX:
	    return new Variant(((SparseMatrix)m_object).transpose());
	}
	System.out.println("Variant::transpose["+m_type+"] not supported.");
	return new Variant();
    }

    public Variant T() {return transpose();}  // an alias

    public Variant append(Object o) {
	return append(new Variant(o));
    }

    public Variant append(Variant B) {
	switch(m_type) {
	case BOOL:
	    Variant[] B3 = B.getArray();
	    Boolean[] C3 = new Boolean[1+B3.length];
	    C3[0] = (Boolean)m_object;
	    for(int i = 0; i < B3.length; i++) C3[i+1] = B3[i].isTrue();
	    m_object = C3;
	    m_type   = VTYPE.MASK;
	    return this;
	case INT:
	    Variant[] B4 = B.getArray(); 
	    Variant[] C4 = new Variant[1+B4.length];
	    C4[0] = new Variant(this.toInteger());
	    for(int i = 0; i < B4.length; i++) C4[i+1] = new Variant(B4[i].toInteger());
	    m_object = C4;
	    m_type   = VTYPE.ARRAY;
	    return this;
	case DOUBLE:
	    Variant[] B5 = B.getArray();
	    Variant[] C5 = new Variant[1+B5.length];
	    C5[0] = new Variant(this.toDouble());
	    for(int i = 0; i < B5.length; i++) C5[i+1] = new Variant(B5[i].toDouble());
	    m_object = C5;
	    m_type   = VTYPE.ARRAY;
	    return this;
	case STRING:
	    Variant[] B6 = B.getArray();
	    Variant[] C6 = new Variant[1+B6.length];
	    C6[0] = new Variant(this.toString());
	    for(int i = 0; i < B6.length; i++) C6[i+1] = new Variant(B6[i].toString());
	    m_object = C6;
	    m_type   = VTYPE.ARRAY;
	    return this;
	case MASK:
	    switch(B.m_type) {
	    case NULL:
	    case BOOL:
	    case INT:
	    case DOUBLE:
	    case STRING:
	    case ARRAY:
		return append(createMask(B));
	    case MASK:
		Boolean[] A2 = (Boolean[])m_object;
		Boolean[] B2 = (Boolean[])B.m_object;
		Boolean[] C  = new Boolean[A2.length + B2.length];
		for(int i = 0; i < A2.length; i++) C[i]           = A2[i];
		for(int i = 0; i < B2.length; i++) C[i+A2.length] = B2[i];
		m_object = C;
		return this;
	    }
	case ARRAY:
	    switch(B.m_type) {
	    case NULL:
	    case BOOL:
	    case INT:
	    case DOUBLE:
	    case STRING:
	    case MASK:
	    case ARRAY:
		Variant[] A2 = (Variant[])m_object;
		Variant[] B2 = B.getArray(); 
		Variant[] C = new Variant[A2.length + B2.length];
		for(int i = 0; i < A2.length; i++) C[i]           = A2[i];
		for(int i = 0; i < B2.length; i++) C[i+A2.length] = B2[i];
		m_object = C;
		return this;
	    }
	    break;
	case INT_MAP:
	case STRING_MAP:
	    HashMap<String,Variant> a_map = (HashMap<String, Variant>)m_object;
	    switch(B.m_type) {
	    case BOOL:
	    case INT:
	    case DOUBLE:
	    case STRING:
	    case MASK:
	    case ARRAY:
		for(Object key: a_map.keySet()) ((Variant)a_map.get(key)).append(B);
		return this;
	    case INT_MAP:
	    case STRING_MAP:
		HashMap<?,?> b_map = (HashMap<?, ?>)B.m_object;
		for(Object key: a_map.keySet()) 
		    if(b_map.containsKey(key)) 
			((Variant)a_map.get(key)).append((Variant)b_map.get(key));
		for(Object key: b_map.keySet())
		    if(!a_map.containsKey(key)) 
			a_map.put((String)key, (Variant)b_map.get(key));
		return this;
	    }
	    break;
	case INT_SET:
	    HashSet<Integer> int_set = (HashSet<Integer>)m_object;
	    switch(B.m_type) {
	    case BOOL:
	    case INT:
	    case DOUBLE:
	    case STRING:
		int_set.add(B.toInteger());
		return this;
	    case INT_SET:
		for(Integer i: B.toIntSet()) int_set.add(i);
		return this;
	    case ARRAY:
		append(new Variant(B.toIntSet(), VTYPE.INT_SET));
		return this;
	    }
	    break;
	case STRING_SET:
	    HashSet<String> str_set = (HashSet<String>)m_object;
	    switch(B.m_type) {
	    case BOOL:
	    case INT:
	    case DOUBLE:
	    case STRING:
		str_set.add(B.toString());
		return this;
	    case STRING_SET:
		for(String s: B.toStringSet()) str_set.add(s);
		return this;
	    case ARRAY:
		append(new Variant(B.toStringSet(), VTYPE.STRING_SET));
		return this;
	    }
	    break;
	}
	System.out.println("Variant::append unsupported for types: " + m_type + ", " + B.m_type);
	return new Variant();
    }

    public static Variant cat(Variant A, Variant B) {
	Variant[] AB = new Variant[2];
	AB[0] = A;
	AB[1] = B;
	return (new Variant(AB)).cat();
    }

    public Variant cat() {
	switch(m_type) {
	case ARRAY:
	    int     N     = size();
	    Variant dims  = shapes().sizes();
	    int     max   = dims.max().toInteger();
	    Variant sizes = sizes();

	    if(max == 1) {
		int       M = sizes().sum().toInteger();
		Variant[] A = new Variant[M];
		int       j = 0;
		for(int i = 0; i < N; i++) {
		    Variant a = get(i);
		    switch(a.m_type) {
		    case BOOL:
		    case INT:
		    case DOUBLE:
		    case STRING:
		    case LRV:
		    case NRV:
			A[j++] = a;
			break;
		    case ARRAY:
		    case INDEX:
		    case VECTOR:
			for(int k = 0; k < sizes.get(i).toInteger(); k++) A[j++] = a.get(k);
			break;
		    default:
			System.out.println("Variant::cat["+a.m_type+"] unsupported.");
			return new Variant();
		    }
		}
		return new Variant(A);
	    } else {
		int       M = sizes.max().toInteger();
		Variant[] A = new Variant[M];
		for(int j = 0; j < M; j++) {
		    Variant[] B = new Variant[N];   // prepare for recursive cat()
		    for(int i = 0; i < N; i++) {
			Variant b = get(i);
			switch(b.m_type) {
			case BOOL:
			case INT:
			case DOUBLE:
			case STRING:
			case LRV:
			case NRV:
			    B[i] = b;
			    break;
			case ARRAY:
			case INDEX:
			case VECTOR:
			    int size = sizes.get(i).toInteger();
			    if(size == 1) B[i] = b.get(0);
			    else          B[i] = b.get(j);
			    break;
			default:
			    System.out.println("Variant::cat["+b.m_type+"] unsupported.");
			    return new Variant();
			}
		    }
		    A[j] = (new Variant(B)).cat();
		}
		return new Variant(A);
	    }
	}
	System.out.println("Variant::cat["+m_type+"] unsupported");
	return new Variant();
    }

    public Variant shapes() {
	switch(m_type) {
	case ARRAY:
	    int N = size();
	    Variant[] s = new Variant[N];
	    for(int i = 0; i < N; i++) s[i] = get(i).shape();
	    return new Variant(s);
	}
	System.out.println("Variant::shapes["+m_type+"] unsupported");
	return new Variant();
    }

    public Variant sizes() {
	switch(m_type) {
	case INDEX:
	case ARRAY:
	    int N = size();
	    Integer[] s = new Integer[N];
	    for(int i = 0; i < N; i++) s[i] = get(i).size();
	    return new Variant(s);
	}
	System.out.println("Variant::sizes["+m_type+"] unsupported");
	return new Variant();
    }

    public Variant sort() { // arrays not sorted in-place
	switch(m_type) {
	case INDEX:
	    Integer[] I = (Integer[])(copy().m_object);
	    Arrays.sort(I);
	    return new Variant(I);
	case ARRAY:
	    Variant[] A = (Variant[])(copy().m_object);
	    Arrays.sort(A, new VariantComparator());
	    return new Variant(A);
	case NRS:
	    return new Variant(getRandomSample().sort()); 
	}
	System.out.println("Variant::sort unsupported for type: " + m_type);
	return new Variant();
    }

    public Variant isort() { // arrays not sorted, just return sort INDEX.
	switch(m_type) {
	case INDEX:
	    Integer[] B = (Integer[])(copy().m_object);
	    Variant[] W = new Variant[B.length];
	    for(int i = 0; i < B.length; i++) W[i] = new Variant(B[i]);
	    return (new Variant(W)).isort();
	case ARRAY:
	    Variant[] A = (Variant[])(copy().m_object);
	    HashMap<Variant, Integer> M = new HashMap<Variant, Integer>();
	    for(int i = 0; i < A.length; i++) M.put(A[i],i);
	    Arrays.sort(A, new VariantComparator());
	    Integer[] I = new Integer[A.length];
	    for(int i = 0; i < A.length; i++) I[i] = M.get(A[i]);
	    return new Variant(I);
	}
	System.out.println("Variant::isort unsupported for type: " + m_type);
	return new Variant();
    }

    static class VariantComparator implements java.util.Comparator {
	public int compare(Object o1, Object o2) {return cmp((Variant)o1, (Variant)o2);}
    }

    public Variant diff()  {return diff (this);}
    public Variant diff0() {return diff0(this);}
    public Variant diff1() {return diff1(this);}

    public Variant sum() {
	switch(m_type) {
	case BOOL:
	    if(isTrue()) return new Variant(1);
	    else         return new Variant(0);
	case INT:
	case DOUBLE:
	    return this;
	case MASK:
	    int M = 0;
	    Boolean[] B2 = (Boolean[])m_object;
	    for(int i = 0; i < B2.length; i++) if(B2[i]) M++;
	    return new Variant(M);
	case ARRAY: 
	    if(isGeneralizedNumbers()) {
		Variant C = new Variant(0);
		for(int i = 0; i < size(); i++) C = add(C, get(i));
		return C;
	    } else {
		int N = size();
		Variant[] C = new Variant[N];
		for(int i = 0; i < size(); i++) C[i] = get(i).sum();
		return new Variant(C);
	    }
	case INDEX:
	    Integer[] V3 = (Integer[])m_object;
	    int total = 0;
	    for(int i = 0; i < V3.length; i++) total += V3[i];
	    return new Variant(total);
	case VECTOR:
	    Double[] V2 = (Double[])m_object;
	    Double total2 = 0.0;
	    for(int i = 0; i < V2.length; i++) total2 += V2[i];
	    return new Variant(total2);
	case INT_SET:
	    Integer N = 0;
	    for(Integer i: toIntSet()) N += i;
	    return new Variant(N);
	case LRV:
	    return new Variant(getLazyRandomVariable().P());
	case NRV:
	    return new Variant(getNumericRandomVariable().P());
	case NRS:
	    return new Variant(getRandomSample().sum());
	}
	System.out.println("Variant::sum unsupported type: " + m_type);
	return new Variant();
    }

    public Variant ave() {return mean();}

    public Variant mean() {
	switch(m_type) {
	case INT:
	    return new Variant(toInteger());
	case DOUBLE:
	    return new Variant(toDouble());
	case ARRAY: 
	    if(isGeneralizedNumbers()) {
		Variant S = new Variant(0.0);
		for(int i = 0; i < size(); i++) S = add(S,get(i).mean());
		return div(S, new Variant(size()));
	    } else {
		int N = size();
		Variant[] C = new Variant[N];
		for(int i = 0; i < size(); i++) C[i] = get(i).mean();
		return new Variant(C);
	    }
	case INDEX:
	case VECTOR:
	    return div(sum(), new Variant(size()));
	case LRV:
	    return new Variant(getLazyRandomVariable().E());
	case NRV:
	    return new Variant(getNumericRandomVariable().E());
	}
	System.out.println("Variant::mean unsupported type: " + m_type);
	return new Variant();
    }

    public Variant median() {
	switch(m_type) {
	case INT:
	    return new Variant(toInteger());
	case DOUBLE:
	    return new Variant(toDouble());
	case LRV:
	    return new Variant(getLazyRandomVariable().median());
	case NRV:
	    return new Variant(getNumericRandomVariable().median());
	}
	System.out.println("Variant::median unsupported type: " + m_type);
	return new Variant();
    }

    public Variant variance() {
	switch(m_type) {
	case ARRAY: 
	    if(isGeneralizedNumbers()) {
		Variant center = sub(this, mean());
		return div(mult(center,center).sum(), new Variant(size()-1));
	    } else {
		int N = size();
		Variant[] C = new Variant[N];
		for(int i = 0; i < size(); i++) C[i] = get(i).variance();
		return new Variant(C);
	    }
	case INDEX:
	case VECTOR:
	    Variant center = sub(this, mean());
	    return div(mult(center,center).sum(), new Variant(size()-1));
	case NRV:
	    return new Variant(getNumericRandomVariable().variance());
	}
	System.out.println("Variant::variance not supported for " + m_type);
	return new Variant();
    }
    
    public Variant stdev() {
	switch(m_type) {
	case INT:
	    return new Variant(0);
	case DOUBLE:
	    return new Variant(0.0);
	case LRV:
	    return new Variant(getLazyRandomVariable().stdev());
	case NRV:
	    return new Variant(getNumericRandomVariable().stdev());
	}
	return variance().sqrt();
    }

    public Variant correlation() {
	// I've only defined correlation for a array of parallel numerical arrays
	Variant shape = shape();
	if(shape.size() == 2) {
	    int K = shape.get(0).toInteger();
	    int N = shape.get(1).toInteger();
	    Variant standard = div(sub(this,mean()),stdev()); // should have K vectors.
	    double[][] corr = new double[K][K];
	    for(int i = 0; i < K; i++) 
		for(int j = i + 1; j < K; j++) 
		    corr[i][j] = mult(standard.get(i),standard.get(j)).sum().toDouble()/(N-1);
	    for(int i = 0; i < K; i++)
		for(int j = 0; j <= i; j++) {
		    if(i == j) corr[i][j] = 1;
		    else       corr[i][j] = corr[j][i];
		}
	    return Variant.createDoubleArray(corr);
	}
	System.out.println("Variant::correlation not supported for " + m_type);
	return new Variant();
    }

    public Variant cumsum() {
	// longitudinal (accumulates along arrays)
	switch(m_type) {
	case INT:
	case DOUBLE:
	    return copy();
	case INDEX:
	    Integer[] A2 = (Integer[])m_object;
	    Integer[] B2 = new Integer[A2.length];
	    if(A2.length == 0) return new Variant(B2);
	    B2[0] = A2[0];
	    for(int i = 1; i < A2.length; i++) B2[i] = A2[i] + B2[i-1];
	    return new Variant(B2);
	case VECTOR:
	    Double[] A3 = (Double[])m_object;
	    Double[] B3 = new Double[A3.length];
	    if(A3.length == 0) return createVector(B3);
	    B3[0] = A3[0];
	    for(int i = 1; i < A3.length; i++) B3[i] = A3[i] + B3[i-1];
	    return createVector(B3);
	case ARRAY:
	    int       t = 0;     // dual channel accum
	    double    s = 0.0;   // dual channel accum
	    Variant[] A = (Variant[])m_object;
	    Variant[] B = new Variant[A.length];
	    for(int i = 0; i < A.length; i++) {
		switch(A[i].m_type) {
		case INT:
		    t += A[i].toInteger();
		    B[i] = new Variant(t + (int)s);
		    break;
		case DOUBLE:
		    s += A[i].toDouble();
		    B[i] = new Variant(s + (double)t);
		    break;
		case ARRAY:
		    B[i] = A[i].cumsum();
		    break;
		default:
		    B[i] = A[i];
		}
	    }
	    return new Variant(B);
	case INT_MAP:
	    HashMap<Integer, Variant> int_map = new HashMap<Integer, Variant>();
	    HashMap<Integer, Variant>    imap = (HashMap<Integer, Variant>)m_object;
	    for(Integer key: imap.keySet()) 
		int_map.put(key, imap.get(key).cumsum());
	    return new Variant(int_map, VTYPE.INT_MAP);
	case STRING_MAP:
	    HashMap<String, Variant> str_map = new HashMap<String, Variant>();
	    HashMap<String, Variant>    smap = (HashMap<String, Variant>)m_object;
	    for(String key: smap.keySet()) 
		str_map.put(key, smap.get(key).cumsum());
	    return new Variant(str_map, VTYPE.STRING_MAP);
	}
	System.out.println("Variant::cumsum unsupported type: " + m_type);
	return new Variant();
    }

    public Variant xcumsum() {
	// Accumulates across arrays. Assumes rectangular 2D.
	switch(m_type) {
	case ARRAY:
	    Variant[] A = (Variant[])m_object;
	    Variant[] B = new Variant[A.length];
	    B[0] = A[0].copy();
	    for(int i = 1; i < A.length; i++) B[i] = add(A[i], B[i-1]);
	    return new Variant(B);
	}
	System.out.println("Variant::xcumsum unsupported type: " + m_type);
	return new Variant();
    }

    public Variant exp() {
	return Variant.exp(this);
    }

    public Variant log() {
	return Variant.log(this);
    }

    public Variant add() {
	switch(m_type) {
	case STRING_MAP:
	    HashMap<String, Variant> A2 = (HashMap<String, Variant>)m_object;
	    Variant C2 = new Variant();
	    for(String key: A2.keySet()) C2 = add(C2, A2.get(key)); // NULL is invisible.
	    return C2;
	}
	System.out.println("Variant::add["+m_type+"] unsupported.");
	return new Variant();
    }

    public Variant min() {
	switch(m_type) {
	case BOOL:
	    if(isTrue()) return new Variant(1);
	    else         return new Variant(0);
	case INT:
	case DOUBLE:
	    return this;
	case MASK:
	    int M = 0;
	    Boolean[] B2 = (Boolean[])m_object;
	    for(int i = 0; i < B2.length; i++) if(!B2[i]) return new Variant(0);
	    return new Variant(1);
	case INDEX:
	    Integer[] I = (Integer[])m_object;
	    if(I.length == 0) return (new Variant(0)).setInfinite();
	    int m = I[0];
	    for(int i = 1; i < I.length; i++) if(I[i] < m) m = I[i];
	    return new Variant(m);
	case ARRAY: 
	    if(isNumbers()) {
		Variant C = new Variant(0.0);
		if(size() == 0) return C.setInfinite();
		C = get(0);
		for(int i = 1; i < size(); i++) if(cmp(C, get(i)) > 0) C = get(i);
		return C;
	    } else {
		int N = size();
		Variant[] C = new Variant[N];
		for(int i = 0; i < size(); i++) C[i] = get(i).min();
		return new Variant(C);
	    }
	case INT_SET:
	    Integer N = null;
	    if(size() == 0) {Variant C = new Variant(0.0); return C.setInfinite();}
	    for(Integer i: toIntSet()) if(N == null || N < i) N = i;
	    return new Variant(N);
	}
	System.out.println("Variant::min() unsupported for type: " + m_type);
	return new Variant();
    }

    public Variant max() {
	switch(m_type) {
	case BOOL:
	    if(isTrue()) return new Variant(1);
	    else         return new Variant(0);
	case INT:
	case DOUBLE:
	    return this;
	case MASK:
	    int M = 0;
	    Boolean[] B2 = (Boolean[])m_object;
	    for(int i = 0; i < B2.length; i++) if(B2[i]) return new Variant(1);
	    return new Variant(0);
	case INDEX:
	    Integer[] I = (Integer[])m_object;
	    if(size() == 0) return neg(((new Variant(0)).setInfinite()));
	    int m = I[0];
	    for(int i = 1; i < I.length; i++) if(I[i] > m) m = I[i];
	    return new Variant(m);
	case ARRAY: 
	    if(isNumbers()) {
		Variant C = new Variant(0.0);
		if(size() == 0) return neg(C.setInfinite());
		C = get(0);
		for(int i = 1; i < size(); i++) if(cmp(C, get(i)) < 0) C = get(i);
		return C;
	    } else {
		int N = size();
		Variant[] C = new Variant[N];
		for(int i = 0; i < size(); i++) C[i] = get(i).max();
		return new Variant(C);
	    }
	case INT_SET:
	    Integer N = null;
	    if(size() == 0) {Variant C = new Variant(0.0); return neg(C.setInfinite());}
	    for(Integer i: toIntSet()) if(N == null || i < N) N = i;
	    return new Variant(N);
	}
	System.out.println("Variant::max() unsupported for type: " + m_type);
	return new Variant();
    }

    public Variant sqrt() {
	switch(m_type) {
	case INT:
	case DOUBLE:
	    return new Variant(java.lang.Math.sqrt(toDouble()));
	case ARRAY:
	case INDEX:
	case VECTOR:
	    int N = size();
	    Variant[] C = new Variant[N];
	    for(int i = 0; i < N; i++) C[i] = get(i).sqrt();
	    return new Variant(C);
	case LRV:
	    return new Variant(getLazyRandomVariable().sqrt());
	case NRV:
	    return new Variant(getNumericRandomVariable().sqrt());
	}
	System.out.println("Variant::sqrt not supported for " + m_type);
	return new Variant();
    }

    public Variant validityMask() {
	switch(m_type) {
	case NULL: 
	    return createMask(false);
	case INT:  
	case STRING:
	    return createMask(true);
	case DOUBLE:
	    return createMask(!isNaN());
	case ARRAY:
	    Variant[] B = (Variant[])m_object;
	    Boolean[] C = new Boolean[B.length];
	    for(int i = 0; i < B.length; i++) C[i] = !B[i].isNaN();
	    return new Variant(C);
	case STRING_MAP: // assume named parallel arrays. Will AND cumulatively
	    HashMap<String, Variant> A2 = (HashMap<String, Variant>)m_object;
	    Variant D = new Variant();
	    for(String s: A2.keySet()) {
		Variant V = A2.get(s).validityMask();
		if(D.m_type == VTYPE.NULL) D = V; else D = mult(D, V);
	    }
	    return D;
	}
	System.out.println("Variant::validityMask("+m_type+") unsupported.");
	return new Variant();
    }

    public Variant toMask() {
	switch(m_type) {
	case MASK: return this;
	default:   return createMask(this);
	}
    }

    public int toInteger() {
	switch(m_type) {
	case BOOL:    if(isTrue()) return 1; else return 0;
	case INT:     return (int)((Integer)m_object);
	case DOUBLE:  return ((Double) m_object).intValue();
	case STRING:  
	    try {
		return (int)Integer.parseInt(toString());
	    } catch(Exception e) {
		System.out.println("Variant::toInteger failed: " + e);
		return -1;
	    }
	case INDEX:
	case VECTOR:
	case MASK:
	case ARRAY:
	    return size();

	}
	System.out.println("Variant::toInteger unsupported type: " + m_type);
	return -1;
    }

    public Variant getDouble() {
	double d = toDouble();
	return new Variant(new Double(d));
    }

    public double toDouble() {
	switch(m_type) {
	case BOOL:
	    if(isTrue()) return 1.0;
	    else         return 0.0;
	case INT:     return (double)((Integer)m_object);
	case DOUBLE:  return (double)((Double) m_object);
	case STRING:  
	    try {
		return (double)Double.valueOf(toString());
	    } catch(Exception e) {
		System.out.println("Variant::toDouble failed: " + e);
		return -1.0;
	    }
	}
	System.out.println("Variant::toDouble unsupported type: " + m_type);
	return -1.0;
    }

    public Set<Integer> toIntSet() {
	switch(m_type) {
	case ARRAY:
	    Variant k = keys();
	    if(k.m_type == VTYPE.INT_SET) return (HashSet<Integer>)k.m_object;
	    else System.out.println("Variant::toIntSet unsupported type: " + k.m_type);
	    break;
	case INT_MAP:
	    return (HashSet<Integer>)(keys().m_object);
	case INT_SET:
	    return (HashSet<Integer>)m_object;
	}
	System.out.println("Variant::toIntSet unsupported type: " + m_type);
	return new HashSet<Integer>();
    }

    // NOTE: There is still an inconsistency about toF() and getF()...
    public Variant getStringSet() {
	switch(m_type) {
	case STRING_SET:
	    return this;
	}
	return new Variant(toStringSet(), VTYPE.STRING_SET);
    }

    public Set<String> toStringSet() {
	switch(m_type) {
	case ARRAY:
	    Variant k = keys();
	    if(k.m_type == VTYPE.STRING_SET) return (HashSet<String>)k.m_object;
	    else System.out.println("Variant::toStringSet unsupported type: " + k.m_type);
	    break;
	case STRING_MAP:
	    return (HashSet<String>)(keys().m_object);
	case STRING_SET:
	    return (HashSet<String>)m_object;
	}
	System.out.println("Variant::toStringSet unsupported type: " + m_type);
	return new HashSet<String>();
    }

    public Map<String, Variant> toStringMap() {
	switch(m_type) {
	case STRING_MAP:
	    return (HashMap<String, Variant>)m_object;
	}
	System.out.println("Variant::toStringMap unsupported type: " + m_type);
	return new HashMap<String, Variant>();
    }

    public String toString() {
	switch(m_type) {
	case NULL:    return "NULL";
	case BOOL:
	    if(isTrue()) return "T";
	    else         return "F";
	case INT:     return ((Integer)m_object).toString();
	case DOUBLE:  return ((Double) m_object).toString();
	case STRING:  return (String)m_object;
	case ARRAY: 
	    StringBuffer buf6 = new StringBuffer();
	    Variant[] variants  = (Variant[])m_object;
	    buf6.append(OPEN_PAREN);
	    for(int i = 0; i < variants.length; i++) {
		buf6.append(variants[i]);
		if(i < variants.length - 1) buf6.append(COMMA_SPACE);
	    }
	    buf6.append(CLOSE_PAREN);
	    return buf6.toString();
	case INT_MAP:
	case STRING_MAP:
	    StringBuffer buf2 = new StringBuffer();
	    HashMap<?,?> map = (HashMap<?, ?>)m_object;
	    Set keyset = map.keySet();
	    int count  = 1;
	    buf2.append(OPEN_BRACE);
	    for(Object key: keyset) {
		buf2.append(key);
		buf2.append(COLON_SPACE);
		buf2.append(map.get(key));
		if(count++ < keyset.size()) buf2.append(COMMA_SPACE);
	    }
	    buf2.append(CLOSE_BRACE);
	    return buf2.toString();
	case INT_SET:
	case STRING_SET:
	    StringBuffer buf3 = new StringBuffer();
	    HashSet<?> set = (HashSet<?>)m_object;
	    int count3  = 1;
	    buf3.append(OPEN_BRACE);
	    for(Object key: set) {
		buf3.append(key);
		if(count3++ < set.size()) buf3.append(COMMA_SPACE);
	    }
	    buf3.append(CLOSE_BRACE);
	    return buf3.toString();
	case SLICE:
	    Integer[]    ii     = (Integer[])m_object;
	    StringBuffer buf1   = new StringBuffer();
	    int          count1 = 1;
	    buf1.append(OPEN_PAREN);
	    for(int i = 0; i < ii.length; i++) {
		buf1.append(ii[i]);
		if(count1++ < ii.length) buf1.append(COMMA_SPACE);
	    }
	    buf1.append(CLOSE_PAREN);
	    return buf1.toString();
	case MASK:
	    Boolean[]    b      = (Boolean[])m_object;
	    StringBuffer buf4   = new StringBuffer();
	    int          count4 = 1;
	    buf4.append(OPEN_PAREN);
	    for(int i = 0; i < b.length; i++) {
		if(b[i]) buf4.append('T');
		else     buf4.append('F');
		if(count4++ < b.length) buf4.append(COMMA_SPACE);
	    }
	    buf4.append(CLOSE_PAREN);
	    return buf4.toString();
	case INDEX:
	    StringBuffer buf8  = new StringBuffer();
	    Integer[] values2  = (Integer[])m_object;
	    buf8.append(OPEN_BRACKET);
	    for(int i = 0; i < values2.length; i++) {
		buf8.append(values2[i]);
		if(i < values2.length - 1) buf8.append(COMMA_SPACE);
	    }
	    buf8.append(CLOSE_BRACKET);
	    return buf8.toString();
	case VECTOR:
	    StringBuffer buf7 = new StringBuffer();
	    Double[] values  = (Double[])m_object;
	    buf7.append(OPEN_BRACKET);
	    for(int i = 0; i < values.length; i++) {
		buf7.append(values[i]);
		if(i < values.length - 1) buf7.append(COMMA_SPACE);
	    }
	    buf7.append(CLOSE_BRACKET);
	    return buf7.toString();
	case SPMATRIX:
	    return ((SparseMatrix)m_object).toString();
	case FUNCTION:
	    return ((Function)m_object).toString();
	case LRV:
	    return getLazyRandomVariable().toString();
	case NRV:
	    return getNumericRandomVariable().toString();
	case NRS:
	    return getRandomSample().toString();
	}
	return null; // should never occur.
    }

    public void print() {
	if(m_type == VTYPE.LRV) {getLazyRandomVariable().print(); return;}
	System.out.println(toString());
    }

    public boolean equals(Object o) {
	return equals(new Variant(o));
    }

    public boolean equals(Variant B) {
	if(m_type != B.m_type) return false;
	switch(m_type) {
	case NULL:
	    return true;
	case BOOL:
	    return (Boolean)m_object.equals((Boolean)B.m_object);
	case INT:
	    return ((Integer)m_object).equals((Integer)B.m_object);
	case DOUBLE:
	    return ((Double)m_object).equals((Double)B.m_object);  // tolerance?
	case STRING:
	    return ((String)m_object).equals((String)B.m_object);
	case ARRAY:
	    Variant[] A2 = (Variant[])m_object;
	    Variant[] B2 = (Variant[])B.m_object;
	    if(A2.length != B2.length) return false;
	    for(int i = 0; i < A2.length; i++) if(! A2[i].equals(B2[i])) return false;
	    return true;
	case MASK:
	    Boolean[] A3 = (Boolean[])m_object;
	    Boolean[] B3 = (Boolean[])B.m_object;
	    if(A3.length != B3.length) return false;
	    for(int i = 0; i < A3.length; i++) if(!A3[i].equals(B3[i])) return false;
	    return true;
	case SLICE:
	    Integer[] A7 = (Integer[])m_object;
	    Integer[] B7 = (Integer[])B.m_object;
	    if(A7.length != B7.length) return false;
	    for(int i = 0; i < A7.length; i++) if(A7[i] != B7[i]) return false;
	    return true;
	case INT_SET:
	case STRING_SET:
	    Set<?> A4 = (Set<?>)m_object;
	    Set<?> B4 = (Set<?>)B.m_object;
	    if(A4.size() != B4.size()) return false;
	    for(Object key: A4) if(!B4.contains(key)) return false;
	    return true;
	case INT_MAP:
	    HashMap<Integer, Variant> A5 = (HashMap<Integer, Variant>)m_object;
	    HashMap<Integer, Variant> B5 = (HashMap<Integer, Variant>)B.m_object;
	    if(A5.size() != B5.size()) return false;
	    for(Integer key: A5.keySet()) {
		if(!B5.containsKey(key)) return false;
		if(!A5.get(key).equals(B5.get(key))) return false;
	    }
	    return true;
	case STRING_MAP:
	    HashMap<String, Variant> A6 = (HashMap<String, Variant>)m_object;
	    HashMap<String, Variant> B6 = (HashMap<String, Variant>)B.m_object;
	    if(A6.size() != B6.size()) return false;
	    for(String key: A6.keySet()) {
		if(!B6.containsKey(key)) return false;
		if(!A6.get(key).equals(B6.get(key))) return false;
	    }
	    return true;
	case FUNCTION:
	    // Not yet implemented.
	}
	System.out.println("Variant::equals not supported to type: " + m_type);
	return false;
    }

    public Variant blunt(double x) {
	// make sure that not values less than x.XXX
	switch(m_type) {
	case NRV:
	    return new Variant(getNumericRandomVariable().blunt(x));
	}
	System.out.println("Variant::blunt not supported to type: " + m_type);
	return null;
    }

    ///////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////// Jython Support ///////////////////////////////////

    public Variant  __getitem__(Object key) {return get(new Variant(key));}
    public void     __setitem__(Object key, Object value) 
                               {set(new Variant(key), new Variant(value));} 

    public Variant  __neg__()           {return neg(this);}
    public Variant  __abs__()           {return abs(this);}

    public Variant  __add__(Variant  b) {return add(this,b);}
    public Variant  __add__(PyObject b) {return add(this,new Variant(b));}
    public Variant __radd__(PyObject b) {return add(new Variant(b),this);}

    public Variant  __sub__(Variant  b) {return sub(this,b);}
    public Variant  __sub__(PyObject b) {return sub(this,new Variant(b));}
    public Variant __rsub__(PyObject b) {return sub(new Variant(b),this);}
    
    public Variant  __mul__(Variant  b) {return mult(this,b);}
    public Variant  __mul__(PyObject b) {return mult(this,new Variant(b));}
    public Variant __rmul__(PyObject b) {return mult(new Variant(b),this);}

    public Variant  __div__(Variant  b) {return div(this,b);}
    public Variant  __div__(PyObject b) {return div(this,new Variant(b));}
    public Variant __rdiv__(PyObject b) {return div(new Variant(b),this);}

    public Variant  __mod__(Variant  b) {return mod(this,b);}
    public Variant  __mod__(PyObject b) {return mod(this,new Variant(b));}
    public Variant __rmod__(PyObject b) {return mod(new Variant(b),this);}
    
    public Variant  __pow__(Variant  b) {return pow(this,b);}
    public Variant  __pow__(PyObject b) {return pow(this,new Variant(b));}
    public Variant __rpow__(PyObject b) {return pow(new Variant(b),this);}

    public Variant  __lt__ (Variant  b) {return cmp(this, b,              CMP_FLAVOR.LT);}
    public Variant  __lt__ (PyObject b) {return cmp(this, new Variant(b), CMP_FLAVOR.LT);}

    public Variant  __le__ (Variant  b) {return cmp(this, b,              CMP_FLAVOR.LE);}
    public Variant  __le__ (PyObject b) {return cmp(this, new Variant(b), CMP_FLAVOR.LE);}

    public Variant  __eq__ (Variant  b) {return cmp(this, b,              CMP_FLAVOR.EQ);}
    public Variant  __eq__ (PyObject b) {return cmp(this,new Variant(b),  CMP_FLAVOR.EQ);}

    public Variant  __gt__ (Variant  b) {return cmp(this, b,              CMP_FLAVOR.GT);}
    public Variant  __gt__ (PyObject b) {return cmp(this, new Variant(b), CMP_FLAVOR.GT);}

    public Variant  __ge__ (Variant  b) {return cmp(this, b,              CMP_FLAVOR.GE);}
    public Variant  __ge__ (PyObject b) {return cmp(this, new Variant(b), CMP_FLAVOR.GE);}

    ///////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////// Statics /////////////////////////////////////////
    ////////////// Unary //////////////////////////////////////////////////////////////

    public static Variant neg(Variant A) {
	switch(A.m_type) {
	case BOOL:
	    return new Variant(!(Boolean)A.m_object);
	case INT:
	    return new Variant(-(Integer)A.m_object);
	case DOUBLE:
	    return new Variant(-(Double)A.m_object);
	case ARRAY:
	    Variant[] C = new Variant[A.size()];
	    for(int i = 0; i < A.size(); i++) C[i] = neg(A.get(i));
	    return new Variant(C);
	case MASK:
	    Boolean[] A2 = (Boolean[])A.m_object;
	    Boolean[] D  = new Boolean[A.size()];
	    for(int i = 0; i < A.size(); i++) D[i] = !A2[i];
	    return new Variant(D);
	case INDEX:
	    Integer[] A4 = (Integer[])A.m_object;
	    Integer[] B4 = new Integer[A4.length];
	    for(int i = 0; i < A4.length; i++) B4[i] = -A4[i];
	    return new Variant(B4);
	case VECTOR:
	    Double[] A3 = (Double[])A.m_object;
	    Double[] B3 = new Double[A3.length];
	    for(int i = 0; i < A3.length; i++) B3[i] = -A3[i];
	    return createVector(B3);
	case SPMATRIX:
	    return new Variant(SparseMatrix.mult((SparseMatrix)A.m_object, -1.0));
	case LRV:
	    return new Variant(A.getLazyRandomVariable().neg());
	case NRV:
	    return new Variant(A.getNumericRandomVariable().neg());
	}
	System.out.println("Variant::neg not supported for " + A.m_type);
	return new Variant();
    }

    public static Variant abs(Variant A) {
	switch(A.m_type) {
	case INT:
	    return new Variant(Math.abs((Integer)A.m_object));
	case DOUBLE:
	    return new Variant(Math.abs((Double)A.m_object));
	case ARRAY:
	    int        N = A.size();
	    Variant[] C = new Variant[N];
	    for(int i = 0; i < N; i++) C[i] = abs(A.get(i));
	    return new Variant(C);
	}
	System.out.println("Variant::abs not supported for " + A.m_type);
	return new Variant();
    }

    public static Variant inv(Variant A) {
	switch(A.m_type) {
	case INT:
	    return new Variant(1.0/(Integer)A.m_object);
	case DOUBLE:
	    return new Variant(1.0/(Double)A.m_object);
	case INDEX:
	    Integer[] A4 = (Integer[])A.m_object;
	    Double [] C4 = new Double[A4.length];
	    for(int i = 0; i < A4.length; i++) C4[i] = 1.0/A4[i];
	    return createVector(C4);
	case VECTOR:
	    Double[] A3 = (Double[])A.m_object;
	    Double[] C3 = new Double[A3.length];
	    for(int i = 0; i < A3.length; i++) C3[i] = 1/A3[i];
	    return createVector(C3);
	case ARRAY:
	    Variant[] A2 = (Variant[])A.m_object;
	    Variant[] C2 = new Variant[A2.length];
	    for(int i = 0; i < A2.length; i++) C2[i] = inv(A2[i]);
	    return new Variant(C2);
	case LRV:
	    return new Variant(A.getLazyRandomVariable().reciprocal());
	case NRV:
	    return new Variant(A.getNumericRandomVariable().reciprocal());
	case NRS:
	    return new Variant(A.getRandomSample().reciprocal());
	}
	System.out.println("Variant::inv["+A.m_type+"] not supported.");
	return new Variant();
    }

    public static Variant exp(Variant A) {
	switch(A.m_type) {
	case INT:
	    return new Variant(Math.exp((Integer)A.m_object));
	case DOUBLE:
	    return new Variant(Math.exp((Double)A.m_object));
	case ARRAY:
	    int        N = A.size();
	    Variant[] C = new Variant[N];
	    for(int i = 0; i < N; i++) C[i] = exp(A.get(i));
	    return new Variant(C);
	case LRV:
	    return new Variant(A.getLazyRandomVariable().exp());
	case NRV:
	    return new Variant(A.getNumericRandomVariable().exp());
	}
	System.out.println("Variant::exp type not supported: " + A.m_type);
	return new Variant();
    }

    public static Variant log(Variant A) {
	switch(A.m_type) {
	case INT:
	    return new Variant(Math.log((Integer)A.m_object));
	case DOUBLE:
	    return new Variant(Math.log((Double)A.m_object));
	case ARRAY:
	    int        N = A.size();
	    Variant[] C = new Variant[N];
	    for(int i = 0; i < N; i++) C[i] = log(A.get(i));
	    return new Variant(C);
	case LRV:
	    return new Variant(A.getLazyRandomVariable().log());
	case NRV:
	    return new Variant(A.getNumericRandomVariable().log());
	}
	System.out.println("Variant::log type not supported: " + A.m_type);
	return new Variant();
    }

    public static Variant sqrt(Variant A) {
	switch(A.m_type) {
	case INT:
	    return new Variant(Math.sqrt((Integer)A.m_object));
	case DOUBLE:
	    return new Variant(Math.sqrt((Double)A.m_object));
	case ARRAY:
	    int        N = A.size();
	    Variant[] C = new Variant[N];
	    for(int i = 0; i < N; i++) C[i] = sqrt(A.get(i));
	    return new Variant(C);
	case LRV:
	    return new Variant(A.getLazyRandomVariable().sqrt());
	case NRV:
	    return new Variant(A.getNumericRandomVariable().sqrt());
	}
	System.out.println("Variant::sqrt type not supported: " + A.m_type);
	return new Variant();
    }

    // Returns the N-1 array of differences. Uses diff1 then chops off first element.
    static public Variant diff(Variant A) {
	Variant d = diff1(A);
	return d.get(createSlice("1:"));
    }

    // Forces first element to be zero. Maybe this should be the base case?
    static public Variant diff0(Variant A) {
	Variant d = diff1(A);
	if(d.size() == 0) return d;
	d.set(0, new Variant(0));   // this is herioc. What does "0" mean for all cases? 
	return d;
    }

    // Prepend the diff (N-1 elements) with first element of array.
    static public Variant diff1(Variant A) {
	switch(A.m_type) {
	case VECTOR:
	    Double[] A2 = (Double[])A.m_object;
	    if(A2.length == 0) return createVector(0);
	    Double[] C2 = new Double[A2.length];
	    C2[0] = A2[0];
	    for(int i = 1; i < A2.length; i++) C2[i] = A2[i] - A2[i-1];
	    return createVector(C2);
	case INDEX:
	    Integer[] A3 = (Integer[])A.m_object;
	    if(A3.length == 0) return new Variant(new Integer[0]);
	    Integer[] C3 = new Integer[A3.length];
	    C3[0] = A3[0];
	    for(int i = 1; i < A3.length; i++) C3[i] = A3[i] - A3[i-1];
	    return new Variant(C3);
	case ARRAY:
	    int N = A.size();
	    Variant[] C = new Variant[N];
	    if(N == 0) return new Variant(C);
	    if(A.isNumbers()) {
		Variant X = A.get(0).copy();
		C[0] = X;
		for(int i = 1; i < N; i++) {Variant Y = A.get(i); C[i] = sub(Y, X); X = Y;}
		return new Variant(C);
	    } else {
		for(int i = 0; i < N; i++) C[i] = diff(A.get(i));
		return new Variant(C);
	    }
	}
	System.out.println("Variant::diff type not supported: "+A.m_type);
	return new Variant();
    }

    static public Variant accum(Variant A) {
	switch(A.m_type) {
	case ARRAY:
	    int N = A.size();
	    Variant[] C = new Variant[N];
	    if(N == 0) return new Variant(C);
	    if(A.isGeneralizedNumbers()) {
		Variant X = A.get(0);
		C[0] = X;
		for(int i = 1; i < N; i++) C[i] = add(A.get(i), C[i-1]);
		return new Variant(C);
	    } else {
		for(int i = 0; i < N; i++) C[i] = accum(A.get(i));
		return new Variant(C);
	    }
	}
	System.out.println("Variant::accum type not supported: "+A.m_type);
	return new Variant();
    }

    static public Variant rshift(Variant A) {
	switch(A.m_type) {
	case INT:
	    return new Variant(0);
	case DOUBLE:
	    return new Variant(0.0);
	case MASK:
	    int N3 = A.size();
	    Boolean[] A3 = (Boolean[])A.m_object;
	    Boolean[] C3 = new Boolean[N3];
	    C3[0] = false;
	    for(int i = 1; i < N3; i++) C3[i] = A3[i-1];
	    return new Variant(C3);
	case ARRAY:
	    int N = A.size();
	    Variant[] C = new Variant[N];
	    if(N == 0) return new Variant(C);
	    if(A.isNumbers()) {
		C[0] = new Variant(0);
		for(int i = 1; i < N; i++) C[i] = A.get(i-1);
		return new Variant(C);
	    } else {
		for(int i = 0; i < N; i++) C[i] = rshift(A.get(i));
		return new Variant(C);
	    }
	case INDEX:
	    int N2 = A.size();
	    Integer[] A2 = (Integer[])A.m_object;
	    Integer[] C2 = new Integer[N2];
	    C2[0] = 0;
	    for(int i = 1; i < N2; i++) C2[i] = A2[i-1];
	    return new Variant(C2);
	}
	System.out.println("Variant::rshift type not supported: "+A.m_type);
	return new Variant();
    }

    //// Binary ////////////////////////////////////////////////////////////////////////

    public static Variant add(Variant A, Variant B) {
	switch(A.m_type) {
	case NULL:
	    return B;
	case INT:
	    Integer A2 = (Integer)A.m_object;
	    switch(B.m_type) {
	    case INT:    return new Variant(A2 + (Integer)B.m_object);
	    case DOUBLE: return new Variant(A2 + (Double) B.m_object);
	    case STRING: return new Variant(A2 + B.toInteger());
	    case ARRAY:
		int       N = B.size();
		Variant[] C = new Variant[N];
		for(int i = 0; i < N; i++) C[i] = add(A, B.get(i));
		return new Variant(C);
	    case MASK:
		Boolean[] D  = new Boolean[B.size()];
		Boolean[] B2 = (Boolean[])B.m_object;
		for(int i = 0; i < B.size(); i++) D[i] = A.isTrue() || B2[i];
		return new Variant(D);
	    case INDEX:
		Integer[] B3 = (Integer[])B.m_object;
		Integer[] C3 = new Integer[B3.length];
		for(int i = 0; i < B3.length; i++) C3[i] = A2 + B3[i];
		return new Variant(C3);
	    case VECTOR:
		Double[] B4 = (Double[])B.m_object;
		Double[] C4 = new Double[B4.length];
		for(int i = 0; i < B4.length; i++) C4[i] = A2 + B4[i];
		return createVector(C4);
	    case INT_SET:
		HashSet<Integer> S = new HashSet<Integer>();
		for(Integer i: B.toIntSet()) S.add(i + A.toInteger());
		return new Variant(S, VTYPE.INT_SET);
	    case STRING_SET: 
	    case LRV:
	    case NRV:
	    case NRS:
		return add(B, A);
	    }
	    break;
	case DOUBLE:
	    Double Ad = (Double)A.m_object;
	    switch(B.m_type) {
	    case INT:    return add(B,A);
	    case DOUBLE: return new Variant(Ad + (Double)B.m_object);
	    case STRING: return new Variant(Ad + B.toDouble());
	    case MASK: 
		Boolean[] D  = new Boolean[B.size()];
		Boolean[] B2 = (Boolean[])B.m_object;
		for(int i = 0; i < B.size(); i++) D[i] = A.isTrue() || B2[i];
		return new Variant(D);
	    case INDEX:
		Integer[] B3 = (Integer[])B.m_object;
		Double [] C3 = new Double [B3.length];
		for(int i = 0; i < B3.length; i++) C3[i] = Ad + B3[i];
		return createVector(C3);
	    case VECTOR:
		Double [] B4 = (Double[])B.m_object;
		Double [] C4 = new Double[B4.length];
		for(int i = 0; i < B4.length; i++) C4[i] = Ad + B4[i];
		return createVector(C4);
	    case ARRAY:
		int       N = B.size();
		Variant[] C = new Variant[N];
		for(int i = 0; i < N; i++) C[i] = add(A, B.get(i));
		return new Variant(C);
	    case INT_SET:
		HashSet<Integer> S = new HashSet<Integer>();
		for(Integer i: B.toIntSet()) S.add(i + A.toInteger());
		return new Variant(S, VTYPE.INT_SET);
	    case STRING_SET: 
	    case LRV:
	    case NRV:
	    case NRS:
		return add(B, A);
	    }
	    break;
	case STRING:
	    switch(B.m_type) {
	    case BOOL:
	    case INT:
	    case DOUBLE:
	    case STRING:
	    case ARRAY:
		return new Variant(A.toString() + B.toString());
	    case STRING_SET:
		return add(B,A);
	    }
	    break;
	case MASK:
	    Boolean[] A3 = (Boolean[])A.m_object;
	    switch(B.m_type) {
	    case NULL:   return A;
	    case INT:    return add(B, A);
	    case DOUBLE: return add(B, A);
	    case ARRAY:
		Boolean[] C  = new Boolean[A.size()];
		for(int i = 0; i < A.size(); i++) C[i] = A3[i] || B.get(i).isTrue();
		return new Variant(C);
	    case MASK:
		Boolean[] D  = new Boolean[A.size()];
		Boolean[] B2 = (Boolean[])B.m_object;
		for(int i = 0; i < A.size(); i++) D[i] = A3[i] || B2[i];
		return new Variant(D);
	    }
	    break;
	case INDEX:
	    Integer[] A4 = (Integer[])A.m_object;
	    switch(B.m_type) {
	    case INT:    return add(B, A);
	    case DOUBLE: return add(B, A);
	    case MASK:   return add(B, A);
	    case INDEX:
		Integer[] B2 = (Integer[])B.m_object;
		Integer[] C2 = new Integer[A4.length];
		for(int i = 0; i < A4.length; i++) C2[i] = A4[i] + B2[i];
		return new Variant(C2);
	    case ARRAY:  return add(B, A);
	    case VECTOR:
		Double[] B3 = (Double[])B.m_object;
		Double[] C3 = new Double[A4.length];
		for(int i = 0; i < A4.length; i++) C3[i] = A4[i] + B3[i];
		return createVector(C3);
	    }
	    break;
	case VECTOR:
	    Double[] A5 = (Double[])A.m_object;
	    switch(B.m_type) {
	    case INT:    return add(B, A);
	    case DOUBLE: return add(B, A);
	    case MASK:   return add(B, A);
	    case INDEX:  return add(B, A);
	    case VECTOR:
		Double[] B4 = (Double[])B.m_object;
		Double[] C4 = new Double[A5.length];
		for(int i = 0; i < A5.length; i++) C4[i] = A5[i] + B4[i];
		return createVector(C4);
	    case ARRAY:
		Double [] C5 = new Double[A5.length];
		for(int i = 0; i < A5.length; i++) C5[i] = A5[i] + B.get(i).toDouble();
		return createVector(C5);
	    }
	    break;
	case ARRAY:
	    switch(B.m_type) {
	    case INT:    return add(B, A);
	    case DOUBLE: return add(B, A);
	    case MASK:   return add(B, A);
	    case VECTOR: return add(B, A);
	    case ARRAY:
	    case INDEX:
		Variant[] C = new Variant[B.size()];
		for(int i = 0; i < C.length; i++) C[i] = add(A.get(i), B.get(i));
		return new Variant(C);
	    }
	    break;
	case INT_SET:
	    switch(B.m_type) {
	    case INT:     return add(B,A);
	    case DOUBLE:  return add(B,A);
	    case INT_SET: return A.copy().append(B); // set union
	    }
	    break;
	case STRING_SET:
	    switch(B.m_type) {
	    case BOOL:
	    case INT:
	    case DOUBLE:
	    case STRING:
	    case ARRAY:
	    case STRING_SET: 
		return A.copy().append(B); // set union
	    }
	    break;
	case STRING_MAP:
	    switch(B.m_type) {
	    case STRING_MAP: 
		return A.copy().append(B); // set union
	    }
	    break;
	case SPMATRIX:
	    SparseMatrix Asp = (SparseMatrix)A.m_object;
	    switch(B.m_type) {
	    case SPMATRIX:
		return createSparseMatrix(SparseMatrix.add(Asp, (SparseMatrix)B.m_object));
	    }
	    break;
	case LRV:
	    LazyRandomVariable Arv = A.getLazyRandomVariable();
	    switch(B.m_type) {
	    case INT:
	    case DOUBLE:
		return new Variant(Arv.add(B.toDouble()));
	    case LRV:
		return new Variant(Arv.add(B.getLazyRandomVariable()));
	    }
	    break;
	case NRV:
	    NumericRandomVariable Anrv = A.getNumericRandomVariable();
	    switch(B.m_type) {
	    case INT:
	    case DOUBLE:
		return new Variant(Anrv.add(B.toDouble()));
	    case NRV:
		return new Variant(Anrv.add(B.getNumericRandomVariable()));
	    }
	    break;
	case NRS:
	    RandomSample Ars = A.getRandomSample();
	    switch(B.m_type) {
	    case INT:
	    case DOUBLE:
		return new Variant(Ars.add(B.toDouble()));
	    case NRS:
		return new Variant(Ars.add(B.getRandomSample()));
	    }
	    break;
	}
	System.out.println("Variant::add not supported for " + A.m_type + ", " + B.m_type);
	return new Variant();
    }

    public static Variant sub(Variant A, Variant B) {
	switch(A.m_type) {
	case NULL:
	    return neg(B);
	case MASK:
	    switch(B.m_type) {
	    case MASK:
		Variant C = A.copy();
		Boolean[] C2 = (Boolean[])C.m_object;
		for(int i = 0; i < A.size(); i++) if(B.get(i).isTrue()) C2[i] = false;
		return C;
	    }
	    break;
	case INT_SET:
	    switch(B.m_type) {
	    case INT_SET:
		Variant C = A.copy();
		for(Integer i: B.toIntSet()) C.toIntSet().remove(i);
		return C;
	    }
	    break;
	case STRING_SET:
	    switch(B.m_type) {
	    case STRING_SET:
		Variant C = A.copy();
		for(String s: B.toStringSet()) C.toStringSet().remove(s);
		return C;
	    }
	    break;
	case STRING_MAP:
	    switch(B.m_type) {
	    case STRING_MAP: 
		Variant C = createStringMap();
		HashMap<String, Variant> a = (HashMap<String, Variant>)A.m_object;
		HashMap<String, Variant> b = (HashMap<String, Variant>)B.m_object;
		HashMap<String, Variant> c = (HashMap<String, Variant>)C.m_object;
		for(String k: a.keySet()) {
		    if(b.containsKey(k)) continue;
		    c.put(k, a.get(k));
		}
		return C;
	    }
	    break;
	}
	return add(A, neg(B));  // default
    }

    public static Variant mult(Variant A, Variant B) {
	switch(A.m_type) {
	case NULL:
	    return B;
	case INT:
	    switch(B.m_type) {
	    case INT:    return new Variant((Integer)A.m_object * (Integer)B.m_object);
	    case DOUBLE: return new Variant((Integer)A.m_object * (Double) B.m_object);
	    case STRING: return mult(B,A);
	    case ARRAY:
		int       N = B.size();
		Variant[] C = new Variant[N];
		for(int i = 0; i < N; i++) C[i] = mult(A, B.get(i));
		return new Variant(C);
	    case MASK:
	    case INDEX:
	    case VECTOR:
	    case SPMATRIX:
	    case LRV:
	    case NRV:
	    case NRS:
		return mult(B, A);
	    }
	    break;
	case DOUBLE:
	    switch(B.m_type) {
	    case INT:    return new Variant((Double)A.m_object * (Integer)B.m_object);
	    case DOUBLE: return new Variant((Double)A.m_object * (Double) B.m_object);
	    case ARRAY:
		int       N = B.size();
		Variant[] C = new Variant[N];
		for(int i = 0; i < N; i++) C[i] = mult(A, B.get(i));
		return new Variant(C);
	    case MASK:
	    case INDEX:
	    case VECTOR:
	    case SPMATRIX:
	    case LRV:
	    case NRV:
	    case NRS:
		return mult(B, A);
	    }
	    break;
	case STRING:
	    switch(B.m_type) {
	    case INT:
		int       N1 = B.toInteger();
		Variant[] S1 = new Variant[N1];
		for(int i = 0; i < N1; i++) S1[i] = A;
		return new Variant(S1);
	    }
	    break;
	case MASK:
	    Boolean[] A2 = (Boolean[])A.m_object;
	    switch(B.m_type) {
	    case NULL:   return A;
	    case INT:    
		Variant[] E = new Variant[A2.length];
		for(int i = 0; i < A2.length; i++) E[i] = A2[i] ? B.copy() : new Variant(0);
		return new Variant(E);
	    case DOUBLE:
		Variant[] F = new Variant[A2.length];
		for(int i = 0; i < A2.length; i++) F[i] = A2[i] ? B.copy() : new Variant(0.0);
		return new Variant(F);
	    case SPMATRIX:
		Variant[] G = new Variant[A2.length];
		SparseMatrix Bsp   = (SparseMatrix)B.m_object;
		Variant      Ctemp = new Variant(new SparseMatrix(Bsp.rows(), Bsp.cols()));
		for(int i = 0; i < A2.length; i++) G[i] = A2[i] ? B.copy() : Ctemp.copy();
		return new Variant(G);
	    case ARRAY:
		Boolean[] C  = new Boolean[A.size()];
		for(int i = 0; i < A.size(); i++) C[i] = A2[i] && B.get(i).isTrue();
		return new Variant(C);
	    case MASK:
		Boolean[] D  = new Boolean[A.size()];
		Boolean[] B2 = (Boolean[])B.m_object;
		for(int i = 0; i < A.size(); i++) D[i] = A2[i] && B2[i];
		return new Variant(D);
	    }
	    break;
	case ARRAY:
	    switch(B.m_type) {
	    case INT:      return mult(B, A);
	    case DOUBLE:   return mult(B, A);
	    case MASK:     return mult(B, A);
	    case SPMATRIX: return mult(B, A);
	    case VECTOR:   return mult(B, A);
	    case ARRAY:
		int       N = B.size();
		Variant[] C = new Variant[N];
		for(int i = 0; i < N; i++) C[i] = mult(A.get(i), B.get(i));
		return new Variant(C);
	    case INDEX:    return mult(B,A);
	    case LRV:
		Variant[] LRVC = new Variant[A.size()];
		for(int i = 0; i < A.size(); i++) LRVC[i] = mult(A.get(i), B);
		return new Variant(LRVC);
	    }
	    break;
	case INT_SET:
	    switch(B.m_type) {
	    case INT_SET:               // set intersection
		if(A.size() < B.size()) return mult(B, A);
		HashSet<Integer> C = new HashSet<Integer>();
		for(Integer i: A.toIntSet()) if(B.toIntSet().contains(i)) C.add(i);
		return new Variant(C, VTYPE.INT_SET);
	    }
	    break;
	case STRING_SET:
	    switch(B.m_type) {
	    case STRING_SET:
		if(A.size() < B.size()) return mult(B, A);
		HashSet<String> C = new HashSet<String>();
		for(String s: A.toStringSet()) if(B.toStringSet().contains(s)) C.add(s);
		return new Variant(C, VTYPE.STRING_SET);
	    }
	    break;
	case STRING_MAP:
	    HashMap<String, Variant> A3 = (HashMap<String, Variant>)A.m_object;
	    switch(B.m_type) {
	    case BOOL:
	    case INT:
	    case DOUBLE:
	    case STRING:
	    case ARRAY:
		HashMap<String, Variant> C3 = new HashMap<String, Variant>();
		for(String key: A3.keySet()) C3.put(key, mult(A3.get(key), B));
		return new Variant(C3, VTYPE.STRING_MAP);
	    }
	    break;
	case INDEX:
	    Integer[] Ai = (Integer[])A.m_object;
	    switch(B.m_type) {
	    case INT:
		Integer[] Bi = new Integer[Ai.length];
		for(int i = 0; i < Bi.length; i++) Bi[i] = Ai[i] * B.toInteger();
		return new Variant(Bi);
	    case DOUBLE:
		Double[] Bd = new Double[Ai.length];
		for(int i = 0; i < Bd.length; i++) Bd[i] = Ai[i] * B.toDouble();
		return createVector(Bd);
	    case ARRAY:
		Variant[] C = new Variant[B.size()];
		for(int i = 0; i < C.length; i++) C[i] = mult(A.get(i), B.get(i));
		return new Variant(C);
	    case INDEX: 
		Integer[] Bi2 = (Integer[])B.m_object;
		Integer[] Ci2 = new Integer[Ai.length];
		for(int i = 0; i < Ai.length; i++) Ci2[i] = Ai[i] * Bi2[i];
		return new Variant(Ci2);
	    case VECTOR: // scale
		Double[] Bv3 = (Double[])B.m_object;
		Double[] Cv3 = new Double[Ai.length];
		for(int i = 0; i < Ai.length; i++) Cv3[i] = Ai[i] * Bv3[i];
		return createVector(Cv3);
	    }
	    break;
	case VECTOR:
	    Double[] Av = (Double[])A.m_object;
	    switch(B.m_type) {
	    case INT:
	    case DOUBLE:
		Double[] Bv = new Double[Av.length];
		for(int i = 0; i < Bv.length; i++) Bv[i] = Av[i] * B.toDouble();
		return createVector(Bv);
	    case SPMATRIX:
		return createVector(SparseMatrix.mult(Av, (SparseMatrix)B.m_object));
	    case INDEX:
		return mult(B, A);
	    case VECTOR: // dot product
		Double[] Bv2 = (Double[])B.m_object;
		Double total = 0.0;
		for(int i = 0; i < Av.length; i++) total += Av[i] * Bv2[i];
		return new Variant(total);
	    case ARRAY:
		Double [] C5 = new Double[Av.length];
		for(int i = 0; i < Av.length; i++) C5[i] = Av[i] * B.get(i).toDouble();
		return createVector(C5);
	    }
	    break;
	case SPMATRIX:
	    SparseMatrix Asp = (SparseMatrix)A.m_object;
	    switch(B.m_type) {
	    case INT:
	    case DOUBLE:
		return new Variant(SparseMatrix.mult(Asp, B.toDouble())); 
	    case INDEX:
	    case ARRAY:                                    
		Variant[] B2 = B.getArray();
		Variant[] C2 = new Variant[B2.length];
		for(int i = 0; i < B2.length; i++) C2[i] = mult(A, B2[i]);
		return new Variant(C2);
	    case VECTOR:
		return createVector(SparseMatrix.mult(Asp, B.getVector()));
	    case MASK:
		return mult(B, A);
	    case SPMATRIX:
		return createSparseMatrix(SparseMatrix.mult(Asp, (SparseMatrix)B.m_object));
	    }
	    break;
	case LRV:
	    LazyRandomVariable Arv = A.getLazyRandomVariable();
	    switch(B.m_type) {
	    case INT:
	    case DOUBLE:
		return new Variant(Arv.multiply(B.toDouble()));
	    case ARRAY:
		Variant[] LRVC = new Variant[B.size()];
		for(int i = 0; i < B.size(); i++) LRVC[i] = mult(A,B.get(i));
		return new Variant(LRVC);
	    case LRV:
		return new Variant(Arv.multiply(B.getLazyRandomVariable()));
	    }
	    break;
	case NRV:
	    NumericRandomVariable Anrv = A.getNumericRandomVariable();
	    switch(B.m_type) {
	    case INT:
	    case DOUBLE:
		return new Variant(Anrv.multiply(B.toDouble()));
	    case NRV:
		return new Variant(Anrv.multiply(B.getNumericRandomVariable()));
	    }
	    break;
	case NRS:
	    RandomSample Ars = A.getRandomSample();
	    switch(B.m_type) {
	    case INT:
	    case DOUBLE:
		return new Variant(Ars.multiply(B.toDouble()));
	    case MASK:
		return new Variant(Ars.multiply(new RandomSample(B.getVector())));
	    case NRS:
		return new Variant(Ars.multiply(B.getRandomSample()));
	    }
	    break;
	}
	System.out.println("Variant::mult not supported for " + A.m_type + ", " + B.m_type);
	return new Variant();
    }

    public static Variant hadamard(Variant A, Variant B) {
	switch(A.m_type) {
	case MASK:
	case ARRAY:
	case INDEX:
	    switch(B.m_type) {
	    case MASK:
	    case ARRAY:
	    case INDEX:
	    case VECTOR:
		return mult(A,B);
	    case SPMATRIX:
		Double[] AV = A.getVector();
		return new Variant(SparseMatrix.hadamard(AV, (SparseMatrix)B.m_object));
	    }
	    break;
	case VECTOR:
	    switch(B.m_type) {
	    case NULL:
	    case INT:
	    case DOUBLE:
	    case MASK:
	    case ARRAY:
		return mult(A,B);
	    case VECTOR:
		Double[] AV = (Double[])A.m_object;
		Double[] BV = (Double[])B.m_object;
		Double[] CV = new Double[AV.length];
		for(int i = 0; i < AV.length; i++) CV[i] = AV[i] * BV[i];
		return createVector(CV);
            case SPMATRIX:
		Double[] AV2 = A.getVector();
		return new Variant(SparseMatrix.hadamard(AV2, (SparseMatrix)B.m_object));
	    }
	    break;
	case SPMATRIX:
	    SparseMatrix AS = (SparseMatrix)A.m_object;
	    switch(B.m_type) {
	    case MASK:
	    case INDEX:
	    case ARRAY:
	    case VECTOR:
		return new Variant(SparseMatrix.hadamard(AS, B.getVector()));
	    }
	    break;
	}
	System.out.println("Variant::hadamard not supported for " + A.m_type + ", " + B.m_type);
	return new Variant();
    }

    public static Variant outer(Variant A, Variant B) {
	switch(A.m_type) {
	case INT:
	case DOUBLE:
	    return mult(A,B);
	case VECTOR:
	    switch(B.m_type) {
	    case INT:
	    case DOUBLE:
		return mult(A,B);
	    case VECTOR:
		// This is inefficient, but technically correct.
		Double[] A1 = A.getVector();
		Double[] B1 = B.getVector();
		SparseMatrix S = new SparseMatrix(A1.length, B1.length);
		for(int r = 0; r < A1.length; r++)
		    for(int c = 0; c < B1.length; c++)
			S.set(r, c, A1[r] * B1[c]);
		return new Variant(S);
	    }
	    break;
	case INDEX:
	case ARRAY:
	    Variant[] C = new Variant[A.size()];
	    for(int i = 0; i < A.size(); i++) C[i] = outer(A.get(i), B);
	    return new Variant(C);
	}
	System.out.println("Variant::outer not supported for " + A.m_type + ", " + B.m_type);
	return new Variant();
    }

    public static Variant div0(Variant A, Variant B) {
	Variant z      = new Variant(0);
	Variant zero   = cmp(B, z, CMP_FLAVOR.EQ); 
	Variant result = div(A,B);
	if(zero.isZero()) return result;
	result.set(zero, new Variant(0.0));
	return result;
    }

    public static Variant div(Variant A, Variant B) {
	// TODO: Support or at least catch set operations.
	switch(A.m_type) {
	case VECTOR:
	    Double[] A2 = (Double[])A.m_object;
	    switch(B.m_type) {
	    case INDEX:
		Integer[] B2 = (Integer[])B.m_object;
		Double [] C2 = new Double[A2.length];
		for(int i = 0; i < A2.length; i++) C2[i] = A2[i] / B2[i];
		return createVector(C2);
	    case VECTOR:
		Double[] B3 = (Double[])B.m_object;
		Double[] C3 = new Double[A2.length];
		for(int i = 0; i < A2.length; i++) C3[i] = A2[i] / B3[i];
		return createVector(C3);
	    }
	}
	return mult(A, inv(B));
    }

    public static Variant mod(Variant A, Variant B) {
	switch(A.m_type) {
	case INT:
	    switch(B.m_type) {
	    case INT:
		return new Variant((Integer)A.m_object % (Integer)B.m_object);
	    }
	    break;
	case ARRAY:
	    Variant[] A2 = (Variant[])A.m_object;
	    switch(B.m_type) {
	    case INT:
		Variant[] C2 = new Variant[A2.length];
		for(int i = 0; i < A2.length; i++) C2[i] = mod(A2[i], B);
		return new Variant(C2);
	    }
	    break;
	case INDEX:
	    Integer[] A3 = (Integer[])A.m_object;
	    switch(B.m_type) {
	    case INT:
		int       B3 = B.toInteger();
		Integer[] C3 = new Integer[A3.length];
		for(int i = 0; i < A3.length; i++) C3[i] = A3[i] % B3;
		return new Variant(C3);
	    }
	    break;
	}
	System.out.println("Variant::mod not supported for " + A.m_type + ", " + B.m_type);
	return new Variant();
    }

    // NB: (3,4)**(5,6) = ((3**5, 3**6),(4**5, 4**6)) ~ 2x2 array
    public static Variant pow(Variant A, Variant B) {
	switch(A.m_type) {
	case INT:
	    Integer I2 = (Integer)A.m_object;
	    switch(B.m_type) {
	    case INT:    return new Variant(Math.pow(I2, (Integer)B.m_object));
	    case DOUBLE: return new Variant(Math.pow(I2, (Double) B.m_object));
	    case ARRAY:
		int       N = B.size();
		Variant[] C = new Variant[N];
		for(int i = 0; i < N; i++) C[i] = pow(A, B.get(i));
		return new Variant(C);
	    case INDEX:
		Integer[] B2 = (Integer[])B.m_object;
		Integer[] C2 = new Integer[B2.length];
		for(int i = 0; i < B2.length; i++) C2[i] = (int)Math.pow(I2,B2[i]);
		return new Variant(C2);
	    case VECTOR:
		Double [] B3 = (Double[])B.m_object;
		Double [] C3 = new Double[B3.length];
		for(int i = 0; i < B3.length; i++) C3[i] = Math.pow(I2,B3[i]);
		return createVector(C3);
	    case SPMATRIX:
		return createSparseMatrix(SparseMatrix.pow(A.toDouble(), (SparseMatrix)B.m_object));
	    case LRV:
		return new Variant(B.getLazyRandomVariable().rpow(A.toDouble()));
	    case NRV:
		return new Variant(B.getNumericRandomVariable().rpow(A.toDouble()));
	    }
	    break;
	case DOUBLE:
	    Double Ad = (Double)A.m_object;
	    switch(B.m_type) {
	    case INT:    return new Variant(Math.pow(Ad, (Integer)B.m_object));
	    case DOUBLE: return new Variant(Math.pow(Ad, (Double) B.m_object));
	    case ARRAY:
		int       N = B.size();
		Variant[] C = new Variant[N];
		for(int i = 0; i < N; i++) C[i] = pow(A, B.get(i));
		return new Variant(C);
	    case INDEX:
		Integer[] I3 = (Integer[])B.m_object;
		Double [] C3 = new Double[I3.length];
		for(int i = 0; i < I3.length; i++) C3[i] = Math.pow(Ad,I3[i]);
		return createVector(C3);
	    case VECTOR:
		Double [] B4 = (Double[])B.m_object;
		Double [] C4 = new Double[B4.length];
		for(int i = 0; i < B4.length; i++) C4[i] = Math.pow(Ad,B4[i]);
		return createVector(C4);
	    case SPMATRIX:
		return createSparseMatrix(SparseMatrix.pow(A.toDouble(), (SparseMatrix)B.m_object));
	    case LRV:
		return new Variant(B.getLazyRandomVariable().rpow(A.toDouble()));
	    case NRV:
		return new Variant(B.getNumericRandomVariable().rpow(A.toDouble()));
	    }
	    break;
	case ARRAY:
	    int       N = A.size();
	    Variant[] C = new Variant[N];
	    switch(B.m_type) {
	    case INT:    
	    case DOUBLE:
	    case ARRAY:
	    case INDEX:
		for(int i = 0; i < N; i++) C[i] = pow(A.get(i), B);
		return new Variant(C);
	    }
	    break;
	case INDEX:
	    Integer[] A2 = (Integer[])A.m_object;
	    switch(B.m_type) {
	    case INT:    
		Integer   B2 = (Integer)B.m_object;
		Integer[] C2 = new Integer[A2.length];
		for(int i = 0; i < A2.length; i++) C2[i] = (int)Math.pow(A2[i], B2);
		return new Variant(C2);
	    case DOUBLE:
		Double   B3 = (Double)B.m_object;
		Double[] C3 = new Double[A2.length];
		for(int i = 0; i < A2.length; i++) C3[i] = Math.pow(A2[i], B3);
		return createVector(C3);
	    case ARRAY:
	    case INDEX:
	    case VECTOR:
		Variant[] C4 = new Variant[A2.length];
		for(int i = 0; i < A2.length; i++) C4[i] = pow(A.get(i), B);
		return new Variant(C4);
	    }
	    break;
	case VECTOR:
	    Double[] A3 = (Double[])A.m_object;
	    switch(B.m_type) {
	    case INT:    
		Integer   B2 = (Integer)B.m_object;
		Double [] C2 = new Double[A3.length];
		for(int i = 0; i < A3.length; i++) C2[i] = Math.pow(A3[i], B2);
		return createVector(C2);
	    case DOUBLE:
		Double   B3 = (Double)B.m_object;
		Double[] C3 = new Double[A3.length];
		for(int i = 0; i < A3.length; i++) C3[i] = Math.pow(A3[i], B3);
		return createVector(C3);
	    case ARRAY:
	    case INDEX:
	    case VECTOR:
		Variant[] C4 = new Variant[A3.length];
		for(int i = 0; i < A3.length; i++) C4[i] = pow(A.get(i), B);
		return new Variant(C4);
	    }
	    break;
	case SPMATRIX:
	    SparseMatrix Asp = (SparseMatrix)A.m_object;
	    switch(B.m_type) {
	    case INT:
		return new Variant(SparseMatrix.pow(Asp, (Integer)B.m_object));
	    }
	    break;
	case LRV:
	    LazyRandomVariable Arv = A.getLazyRandomVariable();
	    switch(B.m_type) {
	    case INT:
		int b = B.toInteger();
		if(b == 0) return new Variant(1);
		LazyRandomVariable Crv = Arv;
		for(int i = 1; i < Math.abs(b); i++) Crv = Crv.multiply(new LazyRandomVariable(Arv));
		if(b < 0) Crv = Crv.reciprocal();
		return new Variant(Crv);
	    case DOUBLE:
		return new Variant(Arv.pow(B.toDouble()));
	    case ARRAY:
		if(B.isIntegers()) return pow(A, B.toIndex());
		break;
	    case INDEX:
		// NB: Assumes nested multiplications
		// Will compute all the intermediate cases and pick the requested ones from this list.
		Integer[] B6 = B.getIndex();
		int N7 = 0;
		for(int i = 0; i < B6.length; i++) if(Math.abs(B6[i]) > N7) N7 = Math.abs(B6[i]);
		LazyRandomVariable[] S = new LazyRandomVariable[N7]; // will only do positive cases, in order.
		S[0] = Arv;
		for(int i = 1; i < N7; i++) S[i] = S[i-1].multiply(new LazyRandomVariable(Arv));
		int N6 = B6.length;
		Variant[] C6 = new Variant[N6];
		for(int i = 0; i < N6; i++) {
		    if(B6[i] < 0)       C6[i] = new Variant(S[-B6[i]-1].reciprocal());
		    else if(B6[i] == 0) C6[i] = new Variant(1.0);
		    else                C6[i] = new Variant(S[ B6[i]-1]);
		}
		return new Variant(C6);
	    case VECTOR:
		int       N5 = B.size();
		Variant[] C5 = new Variant[N5];
		for(int i = 0; i < N5; i++) C5[i] = pow(A, B.get(i));
		return new Variant(C5);
	    case LRV:
		return new Variant(Arv.pow(B.getLazyRandomVariable()));
	    }
	    break;
	case NRV:
	    NumericRandomVariable Anrv = A.getNumericRandomVariable();
	    switch(B.m_type) {
	    case INT:
	    case DOUBLE:
		return new Variant(Anrv.pow(B.toDouble()));
	    case INDEX:
		// NB: Assumes nested multiplications (copy of LRV case above)
		// Will compute all the intermediate cases and pick the requested ones from this list.
		Integer[] B11 = B.getIndex();
		int N11 = 0;
		for(int i = 0; i < B11.length; i++) if(Math.abs(B11[i]) > N11) N11 = Math.abs(B11[i]);
		NumericRandomVariable[] S = new NumericRandomVariable[N11]; // will only do positive cases, in order.
		S[0] = Anrv;
		for(int i = 1; i < N11; i++) S[i] = S[i-1].multiply(new NumericRandomVariable(Anrv));
		Variant[] C11 = new Variant[B11.length];
		for(int i = 0; i < B11.length; i++) {
		    if(B11[i] < 0)       C11[i] = new Variant(S[-B11[i]-1].reciprocal());
		    else if(B11[i] == 0) C11[i] = new Variant(1.0);
		    else                 C11[i] = new Variant(S[ B11[i]-1]);
		}
		return new Variant(C11);
	    case ARRAY:
		if(B.isIntegers()) return pow(A, B.toIndex());
		break;
	    case VECTOR:
		Variant[] C9 = new Variant[B.size()];
		for(int i = 0; i < B.size(); i++) C9[i] = pow(A, B.get(i));
		return new Variant(C9);
	    case NRV:
		return new Variant(Anrv.pow(B.getNumericRandomVariable()));
	    }
	    break;
	}
	System.out.println("Variant::pow not supported for " + A.m_type + ", " + B.m_type);
	return new Variant();
    }

    public static int cmp(Variant A, Variant B) {
	// -1 if A <  B, 1 if A > B, 0 else (does not mean equals, just not comparable)
	// -1 if A <~ B in the case of array compares. i.e. 0 <~ (0,1,0)
	// +1 if A >~ B in the case of array compares. i.e. 1 >~ (0,1,0)
	switch(A.m_type) {
	case BOOL:
	    Boolean ab = (Boolean)A.m_object;
	    switch(B.m_type) {
	    case BOOL:
		Boolean bb = (Boolean)B.m_object;
		if(!ab &&  bb) return -1;
		if( ab && !bb) return  1;
		return 0;
	    case INT:    
		Integer bi = (Integer)B.m_object;
		if((!ab && bi > 0) || ( ab && bi > 1)) return -1;
		if(( ab && bi < 1) || (!ab && bi < 0)) return  1;
		return 0;
	    case DOUBLE:
		Double bd = (Double)B.m_object;
		if((!ab && bd > 0) || ( ab && bd > 1)) return -1;
		if(( ab && bd < 1) || (!ab && bd < 0)) return  1;
		return 0;
	    case STRING: 
		Boolean bs = B.isTrue();
		if(!ab &&  bs) return -1;
		if( ab && !bs) return  1;
		return 0;
	    case ARRAY:
		return cmpArray(A, B);
	    case MASK:  // unhandled, yet.
	    }
	    break;
	case INT:
	    Integer ai = (Integer)A.m_object;
	    switch(B.m_type) {
	    case BOOL: return -cmp(B, A);
	    case INT:
		Integer bi = (Integer)B.m_object;
		if(ai < bi) return -1; if(ai > bi) return +1; return 0;
	    case DOUBLE:
		Double bd = (Double)B.m_object;
		if(ai < bd) return -1; if(ai > bd) return +1; return 0;
	    case STRING:
		Integer bs = B.toInteger();
		if(ai < bs) return -1; if(ai > bs) return +1; return 0;
	    case ARRAY:
		return cmpArray(A, B);
	    case MASK:  // unimplemented, yet.
	    }
	    break;
	case DOUBLE:
	    Double ad = (Double)A.m_object;
	    switch(B.m_type) {
	    case BOOL:  return -cmp(B, A);
	    case INT:   return -cmp(B, A);
	    case DOUBLE:
		Double bd = (Double)B.m_object;
		if(ad < bd) return -1; if(ad > bd) return +1; return 0;
	    case STRING:
		Double bs = B.toDouble();
		if(ad < bs) return -1; if(ad > bs) return +1; return 0;
	    case ARRAY:
		return cmpArray(A, B);
	    case MASK:  // unimplemented, yet
	    }
	    break;
	case STRING:
	    String As = (String)A.m_object;
	    switch(B.m_type) {
	    case BOOL:   return -cmp(B, A);
	    case INT:    return -cmp(B, A);
	    case DOUBLE: return -cmp(B, A);
	    case STRING: 
		return (As).compareTo((String)B.m_object);
	    case ARRAY:
		return cmpArray(A, B);
	    case MASK:
	    }
	    break;
	case ARRAY:
	    switch(B.m_type) {
	    case BOOL:   return -cmpArray(B, A);
	    case INT:    return -cmpArray(B, A);
	    case DOUBLE: return -cmpArray(B, A);
	    case STRING: return -cmpArray(B, A);
	    case ARRAY:
		Variant[] av = (Variant[])A.m_object;
		Variant[] bv = (Variant[])B.m_object;
		int N = Math.min(av.length, bv.length);
		boolean flag = false;
		for(int i = 0; i < N; i++) if(cmp(av[i], bv[i]) > 0) {flag = true; break;}
		if(!flag) return -1;
		flag = false;
		for(int i = 0; i < N; i++) if(cmp(av[i], bv[i]) < 0) {flag = true; break;}
		if(!flag) return +1;
		return 0;
	    case MASK:
	    }
	    break;
	case VECTOR:
	    switch(B.m_type) {
	    case VECTOR:  
	    }
	    break;
	case MASK:
	    switch(B.m_type) {
	    case BOOL:   
	    case INT:    
	    case DOUBLE: 
	    case ARRAY:
	    case MASK:
	    }
	    break;
	}
	System.out.println("Variant::cmp not supported for " + A.m_type + ", " + B.m_type);
	return 0;
    }

    private static int cmpArray(Variant A, Variant B) {
	// Assume B is an array, but A is not.
	switch(B.m_type) {
	case ARRAY:
	    Variant[] bv = (Variant[])B.m_object;
	    boolean flag = false;
	    for(int i = 0; i < bv.length; i++) if(cmp(A, bv[i]) > 0) {flag = true; break;}
	    if(!flag) return -1;
	    flag = false;
	    for(int i = 0; i < bv.length; i++) if(cmp(A, bv[i]) < 0) {flag = true; break;}
	    if(!flag) return +1;
	    return 0;
	}
	System.out.println("Variant::cmpArray not supported for " + A.m_type + ", " + B.m_type);
	return 0;
    }

    public enum CMP_FLAVOR {LT, LE, EQ, GT, GE};

    private static Variant cmp(Variant A, Variant B, CMP_FLAVOR F) {
	// NB: A second switch for [L|N]RV appears after this one,
	switch(A.m_type) {
	case BOOL:
	case INT:
	case DOUBLE:
	case STRING:
	    switch(B.m_type) {
	    case BOOL:
	    case INT:
	    case DOUBLE:
	    case STRING:
		switch(F) {
		case LT: return new Variant(cmp(A, B) < 0);
		case LE: return new Variant(cmp(A, B) <= 0);
		case EQ: return new Variant(cmp(A, B) == 0);
		case GT: return new Variant(cmp(A, B) > 0);
		case GE: return new Variant(cmp(A, B) >= 0);
		}
		break;
	    case ARRAY:
		Variant[] B2 = (Variant[])B.m_object;
		Variant[] C2 = new Variant[B2.length];
		for(int i = 0; i < B2.length; i++) C2[i] = cmp(A, B2[i], F);
		return cmpUpgradeToMask(new Variant(C2));
	    }
	    break;
	case ARRAY:	
		Variant[] A2 = (Variant[])A.m_object;
		Variant[] C2 = new Variant[A2.length];
	    switch(B.m_type) {
	    case BOOL:
	    case INT:
	    case DOUBLE:
	    case STRING:
		for(int i = 0; i < A2.length; i++) C2[i] = cmp(A2[i], B, F);
		return cmpUpgradeToMask(new Variant(C2));
	    case ARRAY:
		for(int i = 0; i < A2.length; i++) C2[i] = cmp(A2[i], B.get(i), F);
		return cmpUpgradeToMask(new Variant(C2));
	    }
	    break;
	case INDEX:
	    Integer[] A4 = (Integer[])A.m_object;
	    Boolean[] C4 = new Boolean[A4.length];
	    switch(B.m_type) {
	    case INT:
		int b = B.toInteger();
		switch(F) { 
		case LT: for(int i = 0; i < A4.length; i++) C4[i] = A4[i] <  b;      break;
		case LE: for(int i = 0; i < A4.length; i++) C4[i] = A4[i] <= b;      break;
		case EQ: for(int i = 0; i < A4.length; i++) C4[i] = A4[i].equals(b); break;
		case GT: for(int i = 0; i < A4.length; i++) C4[i] = A4[i] >  b;      break;
		case GE: for(int i = 0; i < A4.length; i++) C4[i] = A4[i] >= b;      break;
		}
		return new Variant(C4);
	    case DOUBLE:
		Double d = B.toDouble();
		switch(F) { 
		case LT: for(int i = 0; i < A4.length; i++) C4[i] = A4[i] <  d;      break;
		case LE: for(int i = 0; i < A4.length; i++) C4[i] = A4[i] <= d;      break;
		case EQ: for(int i = 0; i < A4.length; i++) C4[i] = A4[i].equals(d); break; 
		case GT: for(int i = 0; i < A4.length; i++) C4[i] = A4[i] >  d;      break;
		case GE: for(int i = 0; i < A4.length; i++) C4[i] = A4[i] >= d;      break;
		}
		return new Variant(C4);
	    case INDEX:
		Integer[] I4 = (Integer[])B.m_object;
		switch(F) { 
		case LT: for(int i = 0; i < A4.length; i++) C4[i] = A4[i] <  I4[i];      break;
		case LE: for(int i = 0; i < A4.length; i++) C4[i] = A4[i] <= I4[i];      break;
		case EQ: for(int i = 0; i < A4.length; i++) C4[i] = A4[i].equals(I4[i]); break;
		case GT: for(int i = 0; i < A4.length; i++) C4[i] = A4[i] >  I4[i];      break;
		case GE: for(int i = 0; i < A4.length; i++) C4[i] = A4[i] >= I4[i];      break;
		}
		return new Variant(C4);
	    case VECTOR:
		Double[] B4 = (Double[])B.m_object;
		switch(F) { 
		case LT: for(int i = 0; i < A4.length; i++) C4[i] = A4[i] <  B4[i];      break;
		case LE: for(int i = 0; i < A4.length; i++) C4[i] = A4[i] <= B4[i];      break;
		case EQ: for(int i = 0; i < A4.length; i++) C4[i] = A4[i].equals(B4[i]); break;
		case GT: for(int i = 0; i < A4.length; i++) C4[i] = A4[i] >  B4[i];      break;
		case GE: for(int i = 0; i < A4.length; i++) C4[i] = A4[i] >= B4[i];      break;
		}
		return new Variant(C4);
	    }
	    break;
	case VECTOR:
	    Double [] A3 = (Double[])A.m_object;
	    Boolean[] C3 = new Boolean[A3.length];
	    switch(B.m_type) {
	    case INT:
	    case DOUBLE:
		Double d = B.toDouble();
		switch(F) { 
		case LT: for(int i = 0; i < A3.length; i++) C3[i] = A3[i] <  d;      break;
		case LE: for(int i = 0; i < A3.length; i++) C3[i] = A3[i] <= d;      break;
		case EQ: for(int i = 0; i < A3.length; i++) C3[i] = A3[i].equals(d); break;
		case GT: for(int i = 0; i < A3.length; i++) C3[i] = A3[i] >  d;      break;
		case GE: for(int i = 0; i < A3.length; i++) C3[i] = A3[i] >= d;      break;
		}
		return new Variant(C3);
	    case INDEX:
		Integer[] I3 = (Integer[])B.m_object;
		switch(F) { 
		case LT: for(int i = 0; i < A3.length; i++) C3[i] = A3[i] <  I3[i];      break;
		case LE: for(int i = 0; i < A3.length; i++) C3[i] = A3[i] <= I3[i];      break;
		case EQ: for(int i = 0; i < A3.length; i++) C3[i] = A3[i].equals(I3[i]); break;
		case GT: for(int i = 0; i < A3.length; i++) C3[i] = A3[i] >  I3[i];      break;
		case GE: for(int i = 0; i < A3.length; i++) C3[i] = A3[i] >= I3[i];      break;
		}
		return new Variant(C3);
	    case VECTOR:
		Double[] B3 = (Double[])B.m_object;
		switch(F) { 
		case LT: for(int i = 0; i < A3.length; i++) C3[i] = A3[i] <  B3[i];      break;
		case LE: for(int i = 0; i < A3.length; i++) C3[i] = A3[i] <= B3[i];      break;
		case EQ: for(int i = 0; i < A3.length; i++) C3[i] = A3[i].equals(B3[i]); break;
		case GT: for(int i = 0; i < A3.length; i++) C3[i] = A3[i] >  B3[i];      break;
		case GE: for(int i = 0; i < A3.length; i++) C3[i] = A3[i] >= B3[i];      break;
		}
		return new Variant(C3);
	    }
	    break;
	}
	// Special handling for RV's
	switch(A.m_type) {
	case INT:
	case DOUBLE:
	    switch(B.m_type) {
	    case LRV:
		switch(F) { 
		case LT: return new Variant(B.getLazyRandomVariable().GT (A.toDouble()));
		case LE: return new Variant(B.getLazyRandomVariable().GTE(A.toDouble())); 
		case EQ: return new Variant(B.getLazyRandomVariable().EQ (A.toDouble()));
		case GT: return new Variant(B.getLazyRandomVariable().LT (A.toDouble()));
		case GE: return new Variant(B.getLazyRandomVariable().LTE(A.toDouble())); 
		}
	    case NRV:
		switch(F) { 
		case LT: return new Variant(B.getNumericRandomVariable().GT (A.toDouble()));
		case LE: return new Variant(B.getNumericRandomVariable().GTE(A.toDouble())); 
		case EQ: return new Variant(B.getNumericRandomVariable().EQ (A.toDouble()));
		case GT: return new Variant(B.getNumericRandomVariable().LT (A.toDouble()));
		case GE: return new Variant(B.getNumericRandomVariable().LTE(A.toDouble())); 
		}
	    }
	    break;
	case LRV:
	    switch(B.m_type) {
	    case INT:
	    case DOUBLE:
		switch(F) { 
		case LT: return new Variant(A.getLazyRandomVariable().LT (B.toDouble()));
		case LE: return new Variant(A.getLazyRandomVariable().LTE(B.toDouble())); 
		case EQ: return new Variant(A.getLazyRandomVariable().EQ (B.toDouble()));
		case GT: return new Variant(A.getLazyRandomVariable().GT (B.toDouble()));
		case GE: return new Variant(A.getLazyRandomVariable().GTE(B.toDouble())); 
		}
	    case LRV:
		switch(F) { 
		case LT: return new Variant(A.getLazyRandomVariable().LT (B.getLazyRandomVariable()));
		case LE: return new Variant(A.getLazyRandomVariable().LTE(B.getLazyRandomVariable()));
		case EQ: return new Variant(A.getLazyRandomVariable().EQ (B.getLazyRandomVariable()));
		case GT: return new Variant(A.getLazyRandomVariable().GT (B.getLazyRandomVariable()));
		case GE: return new Variant(A.getLazyRandomVariable().GTE(B.getLazyRandomVariable()));
		}
	    }
	    break;
	case NRV:
	    switch(B.m_type) {
	    case INT:
	    case DOUBLE:
		switch(F) { 
		case LT: return new Variant(A.getNumericRandomVariable().LT (B.toDouble()));
		case LE: return new Variant(A.getNumericRandomVariable().LTE(B.toDouble())); 
		case EQ: return new Variant(A.getNumericRandomVariable().EQ (B.toDouble()));
		case GT: return new Variant(A.getNumericRandomVariable().GT (B.toDouble()));
		case GE: return new Variant(A.getNumericRandomVariable().GTE(B.toDouble())); 
		}
	    case NRV:
		switch(F) { 
		case LT: return new Variant(A.getNumericRandomVariable().LT (B.getNumericRandomVariable()));
		case LE: return new Variant(A.getNumericRandomVariable().LTE(B.getNumericRandomVariable()));
		case EQ: return new Variant(A.getNumericRandomVariable().EQ (B.getNumericRandomVariable()));
		case GT: return new Variant(A.getNumericRandomVariable().GT (B.getNumericRandomVariable()));
		case GE: return new Variant(A.getNumericRandomVariable().GTE(B.getNumericRandomVariable()));
		}
	    }
	    break;
	case NRS:
	    switch(B.m_type) {
	    case INT:
	    case DOUBLE:
		switch(F) { 
		case LT: return createMask(A.getRandomSample().LT (B.toDouble()));
		case LE: return createMask(A.getRandomSample().LTE(B.toDouble())); 
		case EQ: return createMask(A.getRandomSample().EQ (B.toDouble()));
		case GT: return createMask(A.getRandomSample().GT (B.toDouble()));
		case GE: return createMask(A.getRandomSample().GTE(B.toDouble())); 
		}
	    case NRS:
		switch(F) {  
		case LT: return createMask(A.getRandomSample().LT (B.getRandomSample()));
		case LE: return createMask(A.getRandomSample().LTE(B.getRandomSample()));
		case EQ: return createMask(A.getRandomSample().EQ (B.getRandomSample()));
		case GT: return createMask(A.getRandomSample().GT (B.getRandomSample()));
		case GE: return createMask(A.getRandomSample().GTE(B.getRandomSample()));
		}
	    }
	    break;
	}
	System.out.println("Variant::cmp.FLAVOR not supported for "+A.m_type+", "+B.m_type);
	return new Variant();
    }

    private static Variant cmpUpgradeToMask(Variant A) {
	if(A.isArray() && A.isBooleans()) {
	    Variant[] A2 = (Variant[])A.m_object;
	    Boolean[] B  = new Boolean[A2.length];
	    for(int i = 0; i < B.length; i++) B[i] = A2[i].isTrue();
	    return new Variant(B);
	} else
	    return A;
    }

    public static Variant ones(Object o) {
	return ones(new Variant(o));
    }
    public static Variant ones(Variant A) {
	switch(A.m_type) {
	case INT:
	case DOUBLE:
	    int N = A.toInteger();
	    Variant[] B = new Variant[N];
	    for(int i = 0; i < N; i++) B[i] = new Variant(1);
	    return new Variant(B);
	case ARRAY:
	    Variant[] C = new Variant[A.size()];
	    for(int i = 0; i < A.size(); i++) C[i] = ones(A.get(i));
	    return new Variant(C);
	}
	System.out.println("Variant::ones not supported for type: "+A.m_type);
	return new Variant();
    }
    public static Variant zeros(Variant A) {
	switch(A.m_type) {
	case INT:
	case DOUBLE:
	    return zeros(A.toInteger());
	case ARRAY:
	    Variant[] C = new Variant[A.size()];
	    for(int i = 0; i < A.size(); i++) C[i] = zeros(A.get(i));
	    return new Variant(C);
	}
	System.out.println("Variant::zeros not supported for type: "+A.m_type);
	return new Variant();
    }

    private static Variant zeros(int N) {
	Variant[] B = new Variant[N];
	for(int i = 0; i < N; i++) B[i] = new Variant(0);
	return new Variant(B);
    }

    public static Variant arange(Variant S, Variant D) {
	switch(S.m_type) {
	case INT:
	case DOUBLE:
	    int start = S.toInteger();
	    switch(D.m_type) {
	    case INT:
	    case DOUBLE:
		int duration = D.toInteger();
		Variant[] array = new Variant[duration];
		for(int i = 0; i < duration; i++)
		    array[i] = new Variant(i+start);
		return new Variant(array);
	    case ARRAY:
		int        N = D.size();
		Variant[] C = new Variant[N];
		for(int i = 0; i < N; i++) C[i] = arange(S, D.get(i));
		return new Variant(C);
	    }
	    break;
	case ARRAY:
	    int N = S.size();
	    Variant[] C = new Variant[N];
	    switch(D.m_type) {
	    case INT:
	    case DOUBLE:
		for(int i = 0; i < N; i++) C[i] = arange(S.get(i), D);
		return new Variant(C);
	    case ARRAY:
		for(int i = 0; i < N; i++) C[i] = arange(S.get(i), D.get(i));
		return new Variant(C);
	    }
	}
	System.out.println("Variant::arange not supported for types: "+S.m_type+", "+D.m_type);
	return new Variant();
    }

    public static Variant arange(Variant D) {
	return arange(new Variant(0), D);
    }

    public static Variant min(Variant A, Variant B) {
	switch(A.m_type) {
	case INT:
	    switch(B.m_type) {
	    case INT:    return new Variant(Math.min((Integer)A.m_object, (Integer)B.m_object));
	    case DOUBLE: return new Variant(Math.min((Integer)A.m_object, (Double) B.m_object));
	    case ARRAY:
		int        N = B.size();
		Variant[] C = new Variant[N];
		for(int i = 0; i < N; i++) C[i] = min(A, B.get(i));
		return new Variant(C);
	    }
	    break;
	case DOUBLE:
	    switch(B.m_type) {
	    case INT:    return new Variant(Math.min((Double)A.m_object, (Integer)B.m_object));
	    case DOUBLE: return new Variant(Math.min((Double)A.m_object, (Double) B.m_object));
	    case ARRAY:
		int        N = B.size();
		Variant[] C = new Variant[N];
		for(int i = 0; i < N; i++) C[i] = min(A, B.get(i));
		return new Variant(C);
	    }
	    break;
	case ARRAY:
	    int       N = A.size();
	    Variant[] C = new Variant[N];
	    switch(B.m_type) {
	    case INT:    
	    case DOUBLE:
		for(int i = 0; i < N; i++) C[i] = min(A.get(i), B);
		return new Variant(C);
	    case ARRAY:
		for(int i = 0; i < N; i++) C[i] = min(A.get(i), B.get(i));
		return new Variant(C);
	    }
	}
	System.out.println("Variant::min not supported for " + A.m_type + ", " + B.m_type);
	return new Variant();
    }

    public static Variant max(Variant A, Variant B) {
	switch(A.m_type) {
	case INT:
	    switch(B.m_type) {
	    case INT:    return new Variant(Math.max((Integer)A.m_object, (Integer)B.m_object));
	    case DOUBLE: return new Variant(Math.max((Integer)A.m_object, (Double) B.m_object));
	    case ARRAY:
		int        N = B.size();
		Variant[] C = new Variant[N];
		for(int i = 0; i < N; i++) C[i] = max(A, B.get(i));
		return new Variant(C);
	    }
	    break;
	case DOUBLE:
	    switch(B.m_type) {
	    case INT:    return new Variant(Math.max((Double)A.m_object, (Integer)B.m_object));
	    case DOUBLE: return new Variant(Math.max((Double)A.m_object, (Double) B.m_object));
	    case ARRAY:
		int        N = B.size();
		Variant[] C = new Variant[N];
		for(int i = 0; i < N; i++) C[i] = max(A, B.get(i));
		return new Variant(C);
	    }
	    break;
	case ARRAY:
	    int        N = A.size();
	    Variant[] C = new Variant[N];
	    switch(B.m_type) {
	    case INT:    
	    case DOUBLE:
		for(int i = 0; i < N; i++) C[i] = max(A.get(i), B);
		return new Variant(C);
	    case ARRAY:
		for(int i = 0; i < N; i++) C[i] = max(A.get(i), B.get(i));
		return new Variant(C);
	    }
	}
	System.out.println("Variant::max not supported for " + A.m_type + ", " + B.m_type);
	return new Variant();
    }

    public static Variant where(Variant mask, Variant A, Variant B) {
	switch(mask.m_type) {
	case BOOL:
	case INT:
	case DOUBLE:
	    if(mask.isTrue()) return A; return B;
	case ARRAY:
	    return where(mask.toMask(), A, B);
	case MASK:
	    Boolean[] M = (Boolean[])mask.m_object;
	    int       N = M.length;
	    Variant[] C = new Variant[N];
	    switch(A.m_type) {
	    case INT:
	    case DOUBLE:
		switch(B.m_type) {
		case INT:
		case DOUBLE:
		    for(int i = 0; i < N; i++) C[i] = M[i] ? A : B;
		    return new Variant(C);
		case ARRAY:
		    for(int i = 0; i < N; i++) C[i] = M[i] ? A : B.get(i);
		    return new Variant(C);
		}
	    case ARRAY:
		switch(B.m_type) {
		case INT:
		case DOUBLE:
		    for(int i = 0; i < N; i++) C[i] = M[i] ? A.get(i) : B;
		    return new Variant(C);
		case ARRAY:
		    for(int i = 0; i < N; i++) C[i] = M[i] ? A.get(i) : B.get(i);
		    return new Variant(C);
		}
	    }
	}
	System.out.println("Variant::where not supported for " + 
			   mask.m_type + ", " +A.m_type + ", " + B.m_type);
	return new Variant();
    }

    public static Variant cases(Variant key, Variant map, Variant dflt) {
	switch(map.m_type) {
	case INT_MAP:
	    HashMap<Integer, Variant> int_map = (HashMap<Integer, Variant>)map.m_object;
	    switch(key.m_type) {
	    case INT:
		Integer int_key = (Integer)key.m_object;
		if(int_map.containsKey(int_key)) return int_map.get(int_key);
		else                             return dflt;
	    case ARRAY:
		int N = key.size();
		Variant[] C = new Variant[N];
		for(int i = 0; i < N; i++) C[i] = cases(key.get(i), map, dflt.get(i)).get(i);
		return new Variant(C);
	    }
	    break;
	case STRING_MAP:
	    HashMap<String, Variant> str_map = (HashMap<String, Variant>)map.m_object;
	    switch(key.m_type) {
	    case STRING:
		String str_key = (String)key.m_object;
		if(str_map.containsKey(str_key)) return str_map.get(str_key);
		else                             return dflt;
	    case ARRAY:
		int N = key.size();
		Variant[] C = new Variant[N];
		for(int i = 0; i < N; i++) C[i] = cases(key.get(i), map, dflt.get(i)).get(i);
		return new Variant(C);
	    }
	}
	System.out.println("Variant::cases not supported for " + 
			   key.m_type + ", " +map.m_type + ", " + dflt.m_type);
	return new Variant();
    }

    public static Variant cases(Variant key, Variant map) {
	return cases(key, map, new Variant());
    }

    ///////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////// Services //////////////////////////////////////////

    public static Variant OLS(Variant vX, Variant vY) {
        // Y = (1 X)*b, solve for b .. a simple 2D system. Return B = (b0, b1).

	// System.out.println("OLS X: " + X + ", Y: " + Y); // DEBUG

	Variant    X   = vX.toVector();               // Must be DOUBLE's else possible overflow!
	Variant    Y   = vY.toVector();               //
	double     a   = Y.size();                    // <1,1>   (a b) = (1 X)'(1 X) 
	double     b   = X.sum().toDouble();          // <1,X>   (c d)
	double     c   = b;                           // <X,1>   (a b)^-1 = 1/  ( d -b)
	double     d   = mult(X,X).toDouble();        // <X,X>   (c d)      det (-c  a)
	double     det = a*d - b*c;                  
	double     e   = Y.sum().toDouble();          // <1,Y>   (e) = (1 X)'Y
	double     f   = mult(X,Y).toDouble();        // <X,Y>   (f)
	Variant[]  B   = new Variant[2];
	B[0] = new Variant((d*e - b*f)/det);          // ((1 X)'(1 X))^-1 (1 X)'Y
	B[1] = new Variant((a*f - c*e)/det);

	// System.out.println("OLS B: " + new Variant(B));  // DEBUG

	return new Variant(B);
    }

    public static Variant DCF(Variant r, Variant A) {
	// = A0 + A1/(1+r) + A2/(1+r)^2 + ...
	switch(r.m_type) {
	case INT:
	case DOUBLE:
	    Variant d = Variant.inv(Variant.add(new Variant(1),r));  // 1/(1+r)
	    switch(A.m_type) {
	    case ARRAY:
		if(A.isGeneralizedNumbers()) {
		    Variant D = arange(new Variant(A.size()));
		    return (Variant.mult(A, Variant.pow(d, D))).sum(); // (A * d**D).sum()
		} else {
		    int N = A.size();
		    Variant[] C = new Variant[N];
		    for(int i = 0; i < N; i++) C[i] = DCF(r, A.get(i));
		    return new Variant(C);
		}
	    }
	    break;
	case ARRAY:
	    switch(A.m_type) {
	    case ARRAY:
		int N = A.size();
		Variant[] C = new Variant[N];
		for(int i = 0; i < N; i++) C[i] = DCF(r.get(i), A.get(i));
		return new Variant(C);
	    }
	}
	System.out.println("Variant::DCF type not supported: "+r.m_type+", "+A.m_type);
	return new Variant();
    }
    
    public static Variant PMT(Variant r, Variant n) {
	// pmt = r/((1+r)^{-n}-1) ;  negative by convention for positive loan (assumed $1)
	switch(r.m_type) {
	case INT:
	case DOUBLE:
	    Variant one = new Variant(1);
	    switch(n.m_type) {
	    case INT:
	    case DOUBLE:
		if(r.isZero()) return new Variant(0);
		return Variant.div(r, Variant.sub(Variant.pow(Variant.add(one,r),Variant.neg(n)),one));
	    }
	    break;
	}
	System.out.println("Variant::PMT type not supported: "+r.m_type+", "+n.m_type);
	return new Variant();
    }

    public static Variant NPV(Variant R, Variant CF) {
	// CF = value: n := CF. npv = (1-d^n)/(1-d) == 1+d+d^2+...+d^(n-1);  d = 1/(1+r)
	// CF = array: x := CF. npv = x1*d+x2*d^2+...x{N-1}*d^{N-1}  (if 1-based arrays)
	switch(R.m_type) {
	case INT:
	case DOUBLE:
	    double r = R.toDouble();
	    double d = 1/(1+r);
	    switch(CF.m_type) {
	    case INT:
	    case DOUBLE:
		if(R.isZero()) return CF;
		double n = CF.toDouble();
		return new Variant( (1-Math.pow(d,n))/(1-d) );
	    case ARRAY:
		int N = CF.size();
		if(N == 0) return new Variant(0.0);
		if(CF.isNumbers()) {
		    Double[] X = CF.getVector();
		    double C = 0;
		    double di = 1;
		    for(int i = 0; i < N; i++) {di *= d; C += X[i]*di;}  // discount first value
		    return new Variant(C);
		} else {
		    Variant C = mult(CF.get(0), new Variant(d));
		    double di = d;
		    for(int i = 1; i < N; i++) {di *= d; C = add(C, mult(CF.get(i), new Variant(di)));}
		    return C;
		}
	    }
	    break;
	case ARRAY:
	    int N = CF.size();
	    Variant[] C = new Variant[N];
	    switch(CF.m_type) {
	    case INT:
	    case DOUBLE:
		for(int i = 0; i < N; i++) C[i] = NPV(R.get(i), CF);
		return new Variant(C);
	    case ARRAY:
		for(int i = 0; i < N; i++) C[i] = NPV(R.get(i), CF.get(i));
		return new Variant(C);
	    }
	}
	System.out.println("Variant::NPV type not supported: "+R.m_type+", "+CF.m_type);
	return new Variant();
    }

    public static Variant LDCF(Variant r, Variant b, Variant m, Variant n) {
	// Assume discount: r, cash flows: LDCF = DCF(r,[b+0m, b+1m, b+2m, ..., b+(n-1)m])
	// ldcf = b*[(1-d^n)/(1-d)] + m*[d*(1-d^(n-1))/(1-d)-(n-1)*d^n]/(1-d); d:=1/(1+r)
	if(r.isNumber() && b.isNumber() && m.isNumber() && n.isNumber()) {
	    Variant d = Variant.inv(Variant.add(new Variant(1),r));  // 1/(1+r)
	    Variant one  = new Variant (1);
	    Variant n1   = Variant.sub (n,    one);                  // n-1
	    Variant d1   = Variant.sub (one,  d  );                  // 1-d
	    Variant dn   = Variant.pow (d,    n  );                  // d^n
	    Variant dn1  = Variant.pow (d,    n1 );                  // d^(n-1)
	    Variant d1n  = Variant.sub (one,  dn );                  // 1-d^n
	    Variant d1n1 = Variant.sub (one,  dn1);                  // 1-d^(n-1)
	    Variant n1dn = Variant.mult(n1,   dn );                  // (n-1)*d^n
	    Variant c    = Variant.div (d1n,  d1 );                  // (1-d^n)/(1-d)
	    Variant x1   = Variant.mult(d, Variant.div(d1n1, d1) ); // d*(1-d^(n-1))/(1-d)
	    Variant x    = Variant.div (Variant.sub(x1, n1dn), d1);
	    return Variant.add(Variant.mult(b, c), Variant.mult(m, x));
	} 
	int N = 0;
	if     (r.isArray()) N = r.size();
	else if(b.isArray()) N = b.size();
	else if(m.isArray()) N = m.size();
	else if(n.isArray()) N = n.size();
	else {
	    System.out.println("Variant::LDCF type not supported: "+
			       r.m_type+", "+b.m_type+", "+m.m_type+", "+n.m_type);
	    return new Variant();
	}
	Variant[] C = new Variant[N];
	for(int i = 0; i < N; i++)
	    C[i] = LDCF(r.get(i), b.get(i), m.get(i), n.get(i));
	return new Variant(C);
    }

    public static Variant PLDCF(Variant r,Variant lifetime,Variant periods,Variant accums) {
	Variant x_un = periods.T();
	Variant y    = accums.T();
	Variant x    = min(lifetime, x_un);       // lifetime limited version
	Variant dx   = diff(x);
	Variant dy   = diff(y);
	Variant b    = rshift(y);                 // initial value of segment
	Variant m    = div(dy,dx);
	Variant ldcf = LDCF(r,b,m,dx);            // compute each segment separately
	Variant one  = new Variant(1);
        Variant d    = div(one,add(one, r));      // 1/(1+r)
	Variant e    = rshift(x);
	Variant dcf  = mult(pow(d,e),ldcf).sum(); // d^e . ldcf
	return dcf;
    }

    public static Variant CCF(Variant R, Variant A) {
	// CCF(DOUBLE, ARRAY[DOUBLES]) = (A0, A0*r + A1, ..., A0*r^{n-1}+A1*r^{n-2}+...+A{n-1})
	switch(R.m_type) {
	case INT:
	case DOUBLE:
	    double r = R.toDouble();
	    switch(A.m_type) {
	    case ARRAY:
		if(A.isNumbers()) {
		    Double[] x = A.getVector();
		    double[] y = new double[x.length];
		    y[0] = x[0].doubleValue();
		    for(int i = 1; i < x.length; i++) y[i] = y[i-1]*r + x[i].doubleValue();  // recursion
		    return createDoubleArray(y);
		} else {
		    int N = A.size();
		    Variant[] C = new Variant[N];
		    for(int i = 0; i < N; i++) C[i] = CCF(R, A.get(i));
		    return new Variant(C);
		}
	    }
	    break;
	}
	System.out.println("Variant::CCF type not supported: "+R.m_type+", "+A.m_type);
	return new Variant();
    }

    public static Variant convolve(Variant A, Variant B) {
	// will only return an array the length of the longest array. 
	// Pad that array with zeros if you want the traditional version.
	switch(A.m_type) {
	case INT:
	case DOUBLE:
	    return convolve(A.toArray(), B);
	case INDEX:
	case ARRAY:
	    if(A.isNumbers()) {
		switch(B.m_type) {
		case INT:
		case DOUBLE:
		    return convolve(A, B.toArray());
		case INDEX:
		case ARRAY:
		    if(B.isNumbers()) return convolve_arrays(A.getVector(), B.getVector());
		}
	    }
	}
	System.out.println("Variant::convolve types not supported: "+A.m_type+", "+B.m_type);
	return new Variant();
    }

    private static Variant convolve_arrays(Double[] f, Double[] g) {
	if(f.length > g.length) return convolve_arrays(g,f);
	double[] h = new double[g.length];
	for(int n = 0; n < h.length; n++) {
	    int MIN = Math.max(n - (f.length-1), 0);
	    h[n]    = 0.0;
	    for(int m = MIN; m <= n; m++) {
		h[n] += f[n-m]*g[m];
	       //System.out.println((h[n])+"= h["+n+"] = f["+(n-m)+"]*g["+m+"] = "+(f[n-m]*g[m])+" ... n = "+n+", m = "+m);
	    }
	}
	return Variant.createDoubleArray(h);
    }

    public static Variant choose(Variant N, Variant K) {
	switch(N.m_type) {
	case INT:
	    switch(K.m_type) {
	    case INT:
		return new Variant(MathUtils.binomialCoefficient((Integer)N.m_object, (Integer)K.m_object));
	    case ARRAY:
		int Ksize = K.size();
		Variant[] C = new Variant[Ksize];
		for(int i = 0; i < Ksize; i++) C[i] = choose(N, K.get(i));
		return new Variant(C);
	    }
	}
	System.out.println("Variant::choose types not supported: "+N.m_type+", "+K.m_type);
	return new Variant();
    }
    
    // Cash Flow Value under Depreciation. Length is max(CF length, D length).
    public static Variant CFV(Variant CF, Variant D) {
	Variant S = sub(new Variant(1), rshift(D.cumsum()));
	return convolve(CF, S);
    }

    public static Variant CFD(Variant CF, Variant D) {
	// Return Cash Flows corrected for Depreciation
	switch(CF.m_type) {
	case ARRAY:
	    if(CF.isNumbers()) {
		Double[] cf = CF.getVector();
		int N = cf.length;
		switch(D.m_type) {
		case ARRAY:
		    if(D.isNumbers()) {
			Double[] d = D.getVector();
			int      T = d.length;
			for(int n = 1; n < N; n++) {
			    double S = 0;
			    for(int i = Math.max(0,n-T); i < n; i++) S += cf[i] * d[n-1-i];
			    cf[n] += S;
			}
			return createDoubleArray(cf); 
		    }
		}
	    }
	}
	System.out.println("Variant::CFD types not supported: "+CF.m_type+", "+D.m_type);
	return new Variant();
    }

    public static Variant asset_service_steady_state(Variant CF, Variant D) {
	// Uncorrected cash flows, but some unstated assumptions. Needs looking into.
	Double   Vo = CF.sum().toDouble();
	Double[] d  = D.getVector();
	int      T  = d.length;
	double   S  = 0;
	for(int i = 0; i < T - 1; i++) S += (T - 1 - i) * d[i];  // triangular sum
	return new Variant( Vo / (T - S) );
    }

    // DEPRECATED
    public static Variant asset_service(Variant CF, Variant D, int N) {
	// We compute the schedule of payments (for N periods) required to service a sequence of cash flows CF
	// subject to depreciation schedule D
	switch(CF.m_type) {
	case INT:
	case DOUBLE:
	    switch(D.m_type) {
	    case ARRAY:
		if(D.isNumbers()) {
		    Double[] cd = D.cumsum().getVector();
		    int      T  = cd.length;
		    double[] V  = new double[N];
		    V[0]        = CF.toDouble();
		    for(int n = 1; n < N; n++) {
			double S = 0;
			for(int i = Math.max(0,n-T); i < n; i++) S += V[i] * (1 - cd[n-1-i]);
			V[n] = V[0] - S;
		    }
		    return createDoubleArray(V);
		}
	    }
	case ARRAY:
	    if(CF.isNumbers()) {
		Variant X = asset_service(new Variant(1), D, N);
		return convolve(CF, X);
	    }
	}
	System.out.println("Variant::asset_service types not supported: "+CF.m_type+", "+D.m_type);
	return new Variant();
    }

    /////////////////////////////////////////////////////////////////////////////
    /////////////////////////////// Helper Class /////////////////////////////////////

    static interface Function {
	public Variant copy();
	public Variant get(Variant key);
	public String  toString();
    }

    static class LinearFunction implements Function {
	private Variant m_map;  // INT_MAP to number
	private Variant m_intercept;
	private Variant m_slope;

	private LinearFunction() {}

	public LinearFunction(Variant map) {
	    m_map          = map;
	    Variant keys   = map.keys().toArray();
	    Variant values = map.get(keys);
	    Variant b      = Variant.OLS(keys, values);
	    m_intercept    = b.get(0);
	    m_slope        = b.get(1);
	}

	public Variant copy() {
	    LinearFunction L = new LinearFunction();
	    L.m_map = m_map.copy();
	    L.m_intercept = m_intercept.copy();
	    L.m_slope = m_slope.copy();
	    Variant result = new Variant();
	    result.m_type = VTYPE.FUNCTION;
	    result.m_object = L;
	    return result;
	}
	
	public Variant get(Variant key) {
	    Variant dflt = Variant.add(Variant.mult(m_slope, key), m_intercept); // linear
	    return m_map.get(key, dflt);
	}

	public String toString() {
	    StringBuffer buf = new StringBuffer();
	    buf.append(OPEN_BRACKET);
	    buf.append("LinearFunction");
	    buf.append(COLON_SPACE);
	    buf.append(OPEN_PAREN); 
	    buf.append("intercept");
	    buf.append(COLON_SPACE);
	    buf.append(m_intercept);
	    buf.append(COMMA_SPACE);
	    buf.append("slope");
	    buf.append(COLON_SPACE);
	    buf.append(m_slope);
	    buf.append(CLOSE_PAREN); 
	    buf.append(COMMA_SPACE);
	    buf.append(m_map);
	    buf.append(CLOSE_BRACKET);
	    return buf.toString();
	}

    }

   ///////////////////////////////////////////////////////////////////////////////
}