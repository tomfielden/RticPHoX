package com.tomfielden.random;

import java.util.*;
import java.io.*;

/*
The LazyRandomVariable is meant to be drop-in replacement for the numeric ../random/RandomVariable. 
Through lazy-execution we have the opportunity to capture and rewrite the parse tree 
so that productions such as X+X*Y become X(1+Y) which is the  and sum of independent random variables

Switchable immediate rewrite process. Ex: N(1,2) + N(3,4) = N(1+3, sqrt(2^2+4^2))
*/

public class LazyRandomVariable {

    public enum VTYPE {LEAF, ADD, SUB, MUL, DIV, NEG, INV, LOG, EXP, SQRT, POW};

    static private final String OPEN_PAREN  = "(";
    static private final String CLOSE_PAREN = ")";
    static private final String PLUS_STRING = " + ";
    static private final String MUL_STRING  = " * ";
    static private final String NEWLINE     = "\n\r";

    private static final double ZERO_TOL = 1E-15;
    private static final String INDENT   = "  ";

    ////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Members

    public VTYPE                          m_type;
    public ParametricRandomVariable       m_prv = null;
    public ArrayList<LazyRandomVariable>  m_parents = new ArrayList<LazyRandomVariable>();

    ////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Constructors

    public LazyRandomVariable(LazyRandomVariable X) {                                  // deep copy
	m_type     = X.m_type;
	if(m_type == VTYPE.LEAF) m_prv = new ParametricRandomVariable(X.m_prv);        // want new PID
	for(int i = 0; i < X.size(); i++) append(new LazyRandomVariable(X.get(i)));
    }

    public LazyRandomVariable copy() {  // deep copy of LRV only (i.e. no new PID for LEAF's)
	LazyRandomVariable X = new LazyRandomVariable(type());
	X.m_prv = m_prv;
	for(int i = 0; i < size(); i++) X.append(get(i).copy());
	return X;
    }

    private LazyRandomVariable(VTYPE type) {
	m_type = type;
    }

    private LazyRandomVariable(VTYPE type, ArrayList<LazyRandomVariable> parents) {
	m_type    = type;
	m_parents = parents;
    }

    private LazyRandomVariable(ParametricRandomVariable X) {
	m_type = VTYPE.LEAF;
	m_prv  = X;
    }

    private LazyRandomVariable(double d) {
	m_type = VTYPE.LEAF;
	m_prv  = new ParametricRandomVariable(d);
    }

    private void assumeIdentity(LazyRandomVariable A) {
	m_type    = A.m_type;
	m_prv     = A.m_prv;
	m_parents = A.m_parents;
    }

    // Discrete constructors //////////////////////////////////

    static public LazyRandomVariable Discrete(double[] x, double[] p) {
	return new LazyRandomVariable(ParametricRandomVariable.Discrete(x, p));
    }

    public static LazyRandomVariable Poisson(double lambda) {
	return new LazyRandomVariable(ParametricRandomVariable.Poisson(lambda));
    }
    
    // Continuous constructors ////////////////////////////////

    public static LazyRandomVariable Continuous(double[] x, double[] p) {
	return new LazyRandomVariable(ParametricRandomVariable.Continuous(x, p));
    }

    public static LazyRandomVariable Uniform(double a, double b) {
	return new LazyRandomVariable(ParametricRandomVariable.Uniform(a, b));
    }

    public static LazyRandomVariable Beta(double alpha, double beta) {
	return new LazyRandomVariable(ParametricRandomVariable.Beta(alpha, beta));
    }

    public static LazyRandomVariable Cauchy(double median, double scale) {
	return new LazyRandomVariable(ParametricRandomVariable.Cauchy(median, scale));
    }

    public static LazyRandomVariable ChiSquared(double df) {
	return new LazyRandomVariable(ParametricRandomVariable.ChiSquared(df));
    }

    public static LazyRandomVariable Exponential(double lambda) {
	return new LazyRandomVariable(ParametricRandomVariable.Exponential(lambda));
    }

    public static LazyRandomVariable F(double df1, double df2) {
        return new LazyRandomVariable(ParametricRandomVariable.F(df1, df2));
    }

    public static LazyRandomVariable Normal(double mu, double sigma) {
	return new LazyRandomVariable(ParametricRandomVariable.Normal(mu, sigma));
    }

    public static LazyRandomVariable LogNormal(double mu, double sigma) {
	return new LazyRandomVariable(ParametricRandomVariable.LogNormal(mu, sigma));
    }

    public static LazyRandomVariable Student(double median, double scale, double nu) {
	return new LazyRandomVariable(ParametricRandomVariable.Student(median, scale, nu));
    }

    public static LazyRandomVariable Triangular(double L, double M, double R) {
	return new LazyRandomVariable(ParametricRandomVariable.Triangular(L, M, R));
    }

    //////////////////////////////////////////////////////////////////////////////////
    // Access

    public VTYPE  type()  {return m_type;}

    public boolean isDouble() {
	if(type() != VTYPE.LEAF) return false;
	return get().isDouble();
    }

    public boolean isProperLeaf() {
	return type() == VTYPE.LEAF && ! isDouble();
    }

    public boolean isIndependent() {
	boolean hit = false;
	HashMap<Integer, Integer> PID_tally = new HashMap<Integer, Integer>();
	countPIDs(PID_tally);

	ArrayList<Integer> PIDs = new ArrayList<Integer>();
	for(Integer i: PID_tally.keySet()) PIDs.add(i);

	Collections.sort(PIDs, new Comparator(){
		public int compare(Object o1, Object o2) {
		    Integer a = (Integer) o1;
		    Integer b = (Integer) o2;
		    if(a < b) return -1; if(a == b) return 0; return 1;}
	    });

	for(int i = 0; i < PIDs.size(); i++) {
	    Integer PID = PIDs.get(i);
	    if(PID_tally.get(PID) > 1) {
		hit = true;
		System.out.println("WARNING: PID {"+PID+"} appears "+PID_tally.get(PID)+" times.");
	    }
	}
	return !hit;
    }

    private void countPIDs(HashMap<Integer, Integer> PID_tally) {
	for(int i = 0; i < size(); i++) get(i).countPIDs(PID_tally);
	if(isProperLeaf()) {
	    int pid = getPID();
	    if(PID_tally.containsKey(pid)) PID_tally.put(pid, PID_tally.get(pid) + 1);
	    else                           PID_tally.put(pid, 1);
	}
    }

    public double value() {return get().value();}

    public int    size()  {return m_parents.size();}

    private void                         set(ParametricRandomVariable prv)     { m_prv = prv; }
    public ParametricRandomVariable      get()                                 { if(m_prv == null) eval();return m_prv; }
    public           LazyRandomVariable  get(int i)                            { return m_parents.get(i); }
    public ArrayList<LazyRandomVariable> get(ArrayList<Integer> idx_list)      {
	ArrayList<LazyRandomVariable> sublist = new ArrayList<LazyRandomVariable>();
	for(int i = 0; i < idx_list.size(); i++) sublist.add(get(idx_list.get(i)));
	return sublist;
    }

    public HashSet<LazyRandomVariable> getByPID (Set<Integer> PID_set) {
	HashSet<LazyRandomVariable> set = new HashSet<LazyRandomVariable>();
	for(int i = 0; i < size(); i++) { if(PID_set.contains(getPID(i))) set.add(get(i)); }
	return set;
    }

    public int              getPID()      { if(isProperLeaf()) return get().id(); return 0; } 
    public int              getPID(int i) { return get(i).getPID(); }
    public HashSet<Integer> getPIDset()   {
	HashSet<Integer> set = new HashSet<Integer>();
	for(int i = 0; i < size(); i++) {int PID = getPID(i); if(PID > 0) set.add(PID);}
	return set;
    }
    public ArrayList<HashSet<Integer>> getPIDset(ArrayList<Integer> idx_list) {
	ArrayList<HashSet<Integer>> setlist = new ArrayList<HashSet<Integer>>();
	for(int i = 0; i < idx_list.size(); i++) setlist.add(get(idx_list.get(i)).getPIDset());
	return setlist;
    }

    private void prepend(       LazyRandomVariable p) {m_parents.add(0,p);}
    private void append (       LazyRandomVariable p) {m_parents.add(p);}
    private void replace(int i, LazyRandomVariable p) {m_parents.set(i,p);}

    private void remove (int i)                       {m_parents.remove(i);}
    private void remove (ArrayList<Integer> list)     { // Assume list is sorted ascending
	for(int i = list.size() - 1; i >= 0; i--) remove(list.get(i));
    }

    public void removeByPID (Set<Integer> set) {
	for(int i = size() - 1; i >= 0; i--) if(set.contains(getPID(i))) remove(i);
    }

    //public double[][] getData() {return get().getData();} 

    ////////////////////////////////////////////////////////////////////
    // Info

    public double P()           { return get().P();     }
    public double E()           { return get().E();     }
    public double median()      { return get().median();}
    public double stdev()       { return get().stdev(); }

    public double LT (double x) { return get().LT (x); }
    public double LTE(double x) { return get().LTE(x); }
    public double EQ (double x) { return get().EQ (x); }
    public double GT (double x) { return get().GT (x); }
    public double GTE(double x) { return get().GTE(x); }

    public double LT (LazyRandomVariable B) { return get().LT (B.get()); }
    public double LTE(LazyRandomVariable B) { return get().LTE(B.get()); }
    public double GT (LazyRandomVariable B) { return get().GT (B.get()); }
    public double GTE(LazyRandomVariable B) { return get().GTE(B.get()); }
    public double EQ (LazyRandomVariable B) { return get().EQ (B.get()); }

    ////////////////////////////////////////////////////////////////////
    // Arithmetic Operations

    private static LazyRandomVariable createOperation(VTYPE type, LazyRandomVariable parent) {
	LazyRandomVariable X = new LazyRandomVariable(type);
	X.append(parent);
	return X;
    }

    private static LazyRandomVariable createOperation(VTYPE type, LazyRandomVariable A, LazyRandomVariable B) {
	LazyRandomVariable X = new LazyRandomVariable(type);
	X.append(A);
	X.append(B);
	return X;
    }

    public LazyRandomVariable neg()                { return createOperation(VTYPE.NEG,  this); }
    public LazyRandomVariable reciprocal()         { return createOperation(VTYPE.INV,  this); }
    public LazyRandomVariable log       ()         { return createOperation(VTYPE.LOG,  this); }
    public LazyRandomVariable exp       ()         { return createOperation(VTYPE.EXP,  this); }
    public LazyRandomVariable sqrt      ()         { return createOperation(VTYPE.SQRT, this); }
    public LazyRandomVariable add       (double d) { return createOperation(VTYPE.ADD,  this, new LazyRandomVariable(d)); }
    public LazyRandomVariable radd      (double d) { return createOperation(VTYPE.ADD,  new LazyRandomVariable(d), this); }
    public LazyRandomVariable sub       (double d) { return createOperation(VTYPE.SUB,  this, new LazyRandomVariable(d)); }
    public LazyRandomVariable rsub      (double d) { return createOperation(VTYPE.SUB,  new LazyRandomVariable(d), this); }
    public LazyRandomVariable multiply  (double d) { return createOperation(VTYPE.MUL,  this, new LazyRandomVariable(d)); }
    public LazyRandomVariable rmultiply (double d) { return createOperation(VTYPE.MUL,  new LazyRandomVariable(d), this); }
    public LazyRandomVariable divide    (double d) { return createOperation(VTYPE.DIV,  this, new LazyRandomVariable(d)); }
    public LazyRandomVariable rdivide   (double d) { return createOperation(VTYPE.DIV,  new LazyRandomVariable(d), this); }
    public LazyRandomVariable pow       (double d) { return createOperation(VTYPE.POW,  this, new LazyRandomVariable(d)); }
    public LazyRandomVariable rpow      (double d) { return createOperation(VTYPE.POW,  new LazyRandomVariable(d), this); }

    ////////////////////////////////////////////////////////////////////
    // Binary Arithmetic Operations

    public LazyRandomVariable add     (LazyRandomVariable B) { return createOperation(VTYPE.ADD, this, B); }
    public LazyRandomVariable sub     (LazyRandomVariable B) { return createOperation(VTYPE.SUB, this, B); }
    public LazyRandomVariable multiply(LazyRandomVariable B) { return createOperation(VTYPE.MUL, this, B); }
    public LazyRandomVariable divide  (LazyRandomVariable B) { return createOperation(VTYPE.DIV, this, B); }
    public LazyRandomVariable pow     (LazyRandomVariable B) { return createOperation(VTYPE.POW, this, B); }

    //////////////////////////////// Jython Support ///////////////////////////////

    public LazyRandomVariable __neg__()                      {return neg();}

    public LazyRandomVariable __add__ (double d)                 {return add(d);}
    public LazyRandomVariable __radd__(double d)                 {return radd(d);}
    public LazyRandomVariable __add__ (LazyRandomVariable that)  {return add(that);}

    public LazyRandomVariable __sub__ (double d)                 {return sub(d);}
    public LazyRandomVariable __rsub__(double d)                 {return rsub(d);}
    public LazyRandomVariable __sub__ (LazyRandomVariable that)  {return sub(that);}

    public LazyRandomVariable __mul__ (double d)                 {return multiply(d);}
    public LazyRandomVariable __rmul__(double d)                 {return rmultiply(d);}
    public LazyRandomVariable __mul__ (LazyRandomVariable that)  {return multiply(that);}

    public LazyRandomVariable __div__ (double d)                 {return divide(d);}
    public LazyRandomVariable __rdiv__(double d)                 {return rdivide(d);}
    public LazyRandomVariable __div__ (LazyRandomVariable that)  {return divide(that);}

    public LazyRandomVariable __pow__ (double d)                 {return pow(d);}
    public LazyRandomVariable __rpow__(double d)                 {return rpow(d);}
    public LazyRandomVariable __pow__ (LazyRandomVariable that)  {return pow(that);}

    //////////////////////////////////////////////////////////////////////////////////
    // Display
    
    public void plot(String filename) throws IOException  {}

    public String toString() {
	switch(type()) {
	case LEAF: return m_prv.toString();

	case ADD:  return toStringParents(PLUS_STRING);
	case MUL:  return toStringParents(MUL_STRING);

	case SUB:  return String.format("%s - %s",  get(0).toString(),get(1).toString());
	case DIV:  return String.format("%s / %s",  get(0).toString(),get(1).toString());
	case NEG:  return String.format("-%s",      get(0).toString());
	case INV:  return String.format("1/%s",     get(0).toString());
	case LOG:  return String.format("ln(%s)",   get(0).toString());
	case EXP:  return String.format("exp(%s)",  get(0).toString());
	case SQRT: return String.format("sqrt(%s)", get(0).toString());
	case POW:  return String.format("%s^%s",    get(0).toString(),get(1).toString());
	}

	return "ERROR: LazyRandomVariable.toString() unknown type: "+type();
    }

    private String toStringParents(String OP) {
	StringBuffer buf = new StringBuffer();
	buf.append(OPEN_PAREN);
	for(int i = 0; i < size(); i++) {
	    buf.append(get(i).toString());
	    if(i < size() - 1) buf.append(OP);
	}
	buf.append(CLOSE_PAREN);
	return buf.toString();
    }

    public  void print() {print(0);}
    private void print(int indent) {
	StringBuffer buf = new StringBuffer();
	System.out.print(getIndent(indent));

	if(type() == VTYPE.LEAF) System.out.println(m_prv.toString());
	else {
	    switch(type()) {
	    case ADD:  System.out.println("ADD");  break;
	    case MUL:  System.out.println("MUL");  break;       
	    case SUB:  System.out.println("SUB");  break;
	    case DIV:  System.out.println("DIV");  break;
	    case NEG:  System.out.println("NEG");  break;
	    case INV:  System.out.println("INV");  break;
	    case LOG:  System.out.println("LOG");  break;
	    case EXP:  System.out.println("EXP");  break;
	    case SQRT: System.out.println("SQRT"); break;
	    case POW:  System.out.println("POW");  break;
	    }
	    buf.append(NEWLINE);
	    for(int i = 0; i < size(); i++) get(i).print(indent + 1);
	}
    }

    private String getIndent(int indent) {
	StringBuffer buf = new StringBuffer();
	for(int i = 0; i < indent; i++) buf.append(INDENT);
	return buf.toString();
    }

    private static void print(ArrayList<Integer> list) {
	if(list.size() == 0) return;
	System.out.print(list.get(0));
	for(int i = 1; i < list.size(); i++) System.out.print(","+list.get(i));
	System.out.println();
    }

    private static void print(Set<Integer> set) {
	for(Integer i: set) System.out.print(i+",");
	System.out.println();
    }

    private static void print(HashSet<LazyRandomVariable> set) {
	for(LazyRandomVariable L: set) L.print();
    }

    //////////////////////////////////////////////////////////////////////////////////
    // Evaluating the parse tree

    public void eval() { 

	LazyRandomVariable X = simplify();

	X.eval_simple();
	set(X.m_prv);
    }

    public LazyRandomVariable simplify() {
	LazyRandomVariable X = copy();

	X.replace_negs();
	X.collect();
	X.pummel();
	X.factor();

	return X;
    }

    //////////////////////////////////////////////////////////////////////////////////
    // Misc Operations

    public void replace_negs() {
	if(type() == VTYPE.NEG) {
	    m_type = VTYPE.MUL;
	    append(new LazyRandomVariable(-1));
	}
	for(int i = 0; i < size(); i++) get(i).replace_negs();
    }

    public boolean anneal() {
	boolean hit = false;
	for(int i = 0; i < size(); i++) { if(get(i).anneal()) hit = true; }  // depth-first

	if(size() == 1 && (type() == VTYPE.MUL || type() == VTYPE.ADD)) {assumeIdentity(get(0));hit = true;}

	return hit;
    }

    //////////////////////////////////////////////////////////////////////////////////
    // Collect terms

    public void pummel() {
	while(distribute()) while(collect()) ;
    }

    public boolean collect() {
	// Ensure that at each level we have {0,1} DOUBLES, {0,1} MULL, {0,1} ADD
	if(type() == VTYPE.LEAF) return false;

	boolean hit = false;

	for(int i = 0; i < size(); i++)                { if(get(i).collect())       hit = true;}  // depth-first
	if(type() == VTYPE.ADD || type() == VTYPE.MUL) { if(       collect(type())) hit = true;}

	if(remove_mul_1_add_0()) hit = true; 

	return hit;
    }

    private boolean collect(VTYPE type) {
	boolean hit = false;

	// collapse MUL's or ADD's
	ArrayList<LazyRandomVariable> new_parents = new ArrayList<LazyRandomVariable>();
	for(int i = 0; i < size(); i++) {
	    LazyRandomVariable parent = get(i);
	    if(parent.type() == type) {
		hit = true;
		for(int j = 0; j < parent.size(); j++) new_parents.add(parent.get(j));
	    } else 
		new_parents.add(parent);
	}
	if(hit) m_parents = new_parents;
	
	if(collect_doubles(type)) hit = true;

	return hit;
    }

    private boolean collect_doubles(VTYPE type) {
	ArrayList<ArrayList<Integer>> double_test = identify_doubles();
	ArrayList<Integer>            double_ids  = double_test.get(0);
	ArrayList<Integer>            other_ids   = double_test.get(1);

	if(double_ids.size() < 2) return false;

	double new_double = get(double_ids.get(0)).value();
	if(type == VTYPE.MUL) for(int i = 1; i < double_ids.size(); i++) new_double *= get(double_ids.get(i)).value();
	else                  for(int i = 1; i < double_ids.size(); i++) new_double += get(double_ids.get(i)).value();
	    
	ArrayList<LazyRandomVariable> new_parents = new ArrayList<LazyRandomVariable>();
	new_parents.add(new LazyRandomVariable(new_double));

	for(int i = 0; i < other_ids.size(); i++) new_parents.add(get(other_ids.get(i)));

	m_parents = new_parents;

	anneal();     // a singleton DOUBLE may be all we have left. Need to remove the OP.
	return true;
    }

    private ArrayList<ArrayList<Integer>> identify_doubles() {
	ArrayList<Integer> A    = new ArrayList<Integer>();
	ArrayList<Integer> notA = new ArrayList<Integer>();
	for(int i = 0; i < size(); i++) {
	    if(get(i).isDouble()) A.add(i);
	    else               notA.add(i);
	}
	ArrayList<ArrayList<Integer>> pair = new ArrayList<ArrayList<Integer>>();
	pair.add(A);
	pair.add(notA);
	return pair;
    }

    private boolean remove_mul_1_add_0() {
	// Assume number of DOUBLE's in MUL or ADD parents is {0,1}
	boolean hit = false;

	for(int i = 0; i < size(); i++) { if(get(i).remove_mul_1_add_0()) hit = true; }   // depth-first

	if(type() != VTYPE.ADD && type() != VTYPE.MUL) return hit;

	int idx_double = index_first_double();   // Assume {0,1} DOUBLE's in parents
	if(idx_double < 0) return hit;

	double d = get(idx_double).value();

	if(type() == VTYPE.ADD) if(isZero(d    )) {remove(idx_double); hit = true;}
	if(type() == VTYPE.MUL) if(isZero(d - 1)) {remove(idx_double); hit = true;}

	if(hit) anneal();
	return hit;
    }

    private int index_first_double() {
	for(int i = 0; i < size(); i++) if(get(i).isDouble()) return i;
	return -1;
    }
    //////////////////////////////////////////////////////////////////////////////////
    // Distribute DOUBLE's

    public boolean distribute() {
	// Assume collect() has been called. 
	boolean hit = false;
	if(type() == VTYPE.MUL) {
	    ArrayList<Integer> double_ids = identify_doubles().get(0);                      // ASSUME {0,1}
	    if(double_ids.size() == 1) {
		double d = get(double_ids.get(0)).value();  
		ArrayList<LazyRandomVariable> sums = get(identify_type(VTYPE.ADD).get(0));  // ASSUME {0,1}
		if(sums.size() == 1) {
		    hit = true;
		    LazyRandomVariable sum = sums.get(0);
		    sum.distribute(d);
		    if(size() > 2) remove(double_ids.get(0));
		    else           assumeIdentity(sum);
		}
	    }
	}
	for(int i = 0; i < size(); i++) {if(get(i).distribute()) hit = true;}
	return hit;
    }

    private void distribute(double d) {
	for(int i = 0; i < size(); i++) {
	    if(get(i).isDouble()) {
		replace(i,new LazyRandomVariable(d * get(i).value()));
	    } else {
		LazyRandomVariable p = new LazyRandomVariable(VTYPE.MUL);
		p.append(new LazyRandomVariable(d));
		p.append(get(i));
		replace(i,p);
	    }
	}
    }

    //// Factoring /////

    public void factor() {
	while(factor_deep()) ;
    }

    public boolean factor_deep() {
	boolean hit = false;
	for(int i = 0; i < size(); i++) { if(get(i).factor_deep()) hit = true; }  // depth-first (must be done)
	if(type() != VTYPE.ADD) return hit;
	factor_prep();
	ArrayList<Integer> product_ids = identify_type(VTYPE.MUL).get(0);  // mask
	if(product_ids.size() < 2)  {collect(); return hit;}
	Set<Integer> common_PIDs = intersect(getPIDset(product_ids));
	if(common_PIDs.size() == 0) {collect(); return hit;}
	hit = true;
	HashSet<LazyRandomVariable> common_set =    get(product_ids.get(0)).   getByPID(common_PIDs);
	for(int i = 0; i < product_ids.size(); i++) get(product_ids.get(i)).removeByPID(common_PIDs);
	LazyRandomVariable product = new LazyRandomVariable(VTYPE.MUL);
	for(LazyRandomVariable v: common_set) product.append(v);
	LazyRandomVariable sum = new LazyRandomVariable(VTYPE.ADD, get(product_ids));
	product.append(sum);
	remove(product_ids);
	append(product);
	anneal();
	while(collect()) ;  // was just collect()
	return hit;
    }

    public boolean factor_prep() {
	if(type() != VTYPE.ADD) return false;
	boolean hit = false;
	for(int i = 0; i < size(); i++) {
	    LazyRandomVariable parent = get(i);
	    if(parent.isProperLeaf()) {
		hit = true;
		LazyRandomVariable p = new LazyRandomVariable(VTYPE.MUL);
		p.append(parent);
		p.append(new LazyRandomVariable(1));
		replace(i,p);
	    }
	}
	return hit;
    }

    private static Set<Integer> intersect(ArrayList<HashSet<Integer>> set_list) {
	if(set_list.size() == 0) return new HashSet<Integer>();
	Set<Integer> C = set_list.get(0);
	for(int i = 1; i < set_list.size(); i++) {
	    C = intersect(C, set_list.get(i));
	    if(C.size() == 0) break;
	}
	return C;
    }

    private static HashSet<Integer> intersect(Set<Integer> A, Set<Integer> B) {
	if(A.size() < B.size()) return intersect(B, A);
	HashSet<Integer> C = new HashSet<Integer>();
	for(Integer i: A) if(B.contains(i)) C.add(i);
	return C;
    }

    private ArrayList<ArrayList<Integer>> identify_type(VTYPE type) {
	ArrayList<Integer> A    = new ArrayList<Integer>();
	ArrayList<Integer> notA = new ArrayList<Integer>();
	for(int i = 0; i < size(); i++) {
	    if(get(i).type() == type) A.add(i);
	    else                   notA.add(i);
	}
	ArrayList<ArrayList<Integer>> pair = new ArrayList<ArrayList<Integer>>();
	pair.add(A);
	pair.add(notA);
	return pair;
    }

    //// Actual Eval ////////////////////////////////////////////

    private void eval_simple() {
	if(type() == VTYPE.LEAF) return;

	//for(int i = 0; i < size(); i++) get(i).eval_simple(); // depth-first -- no longer needed.
	
	switch(type()) {
	case ADD: 
	    if     (size() == 0) set(new ParametricRandomVariable(0.0));
	    else if(size() == 1) set(get(0).get());
	    else {
		set(get(0).get().add(get(1).get()));
		for(int i = 2; i < size(); i++)
		    set(get().add(get(i).get()));
	    }
	    return;
	case SUB:
	    if(     size() == 0) set(new ParametricRandomVariable(0.0));
	    else if(size() == 1) set(get(0).get());
	    else if(size() == 2) set(get(0).get().sub(get(1).get()));
	    else System.out.println("ERROR LazyRandomVariable.eval() can't subtract " + size() + " operands.");
	    return;
	case MUL:
	    if     (size() == 0) set(new ParametricRandomVariable(1.0));
	    else if(size() == 1) set(get(0).get());
	    else {
		set(get(0).get().multiply(get(1).get()));
		for(int i = 2; i < size(); i++)
		    set(get().multiply(get(i).get()));
	    }
	    return;
	case NEG:
	    if(size() == 1) set(get(0).get().neg());
	    else System.out.println("ERROR LazyRandomVariable.eval() can't negate " + size() + " operands.");
	    return;
	case INV:
	    if     (size() == 0) set(new ParametricRandomVariable(1.0));
	    else if(size() == 1) set(get(0).get().reciprocal());
	    else System.out.println("ERROR LazyRandomVariable.eval() can't invert " + size() + " operands.");
	    return;
	case LOG:
	    if (size() == 1) set(get(0).get().log());
	    else System.out.println("ERROR LazyRandomVariable.eval() can't compute log of " + size() + " operands.");
	    return;
	case EXP:
	    if (size() == 1) set(get(0).get().exp());
	    else System.out.println("ERROR LazyRandomVariable.eval() can't compute exp of " + size() + " operands.");
	    return;
	case SQRT:
	    if (size() == 1) set(get(0).get().sqrt());
	    else System.out.println("ERROR LazyRandomVariable.eval() can't compute sqrt of " + size() + " operands.");
	    return;
	case POW:
	    if     (size() == 1) set(get(0).get());
	    else if(size() == 2) set(get(0).get().pow(get(1).get()));
	    else System.out.println("ERROR LazyRandomVariable.eval() can't compute pow of " + size() + " operands.");
	    return;
	}

	System.out.println("LazyRandomVariable.eval_simple type "+type()+" not implemented.");
    }

    // Util /////////////////////////////////////////////////////////////////////////////////

    private static String toStringDouble(double d) {
	if(isInt(d))          return String.format("%d", (int)d);
	else if(isInt(d*1E1)) return String.format("%.1f", d);
	else if(isInt(d*1E2)) return String.format("%.2f", d);
	else if(isInt(d*1E3)) return String.format("%.3f", d);
	else if(isInt(d*1E4)) return String.format("%.4f", d);
	else if(isInt(d*1E5)) return String.format("%.5f", d);
	else                  return String.format("%f", d);
    }

    private static boolean isInt(double d) {
	return Math.abs((double)(int)d - d) <= ZERO_TOL;
    }

    private static boolean isZero(double d) {
	return Math.abs(d) <= ZERO_TOL;
    }

} // end LazyRandomVariable
