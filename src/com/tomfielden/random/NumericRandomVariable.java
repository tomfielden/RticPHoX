package com.tomfielden.random;

import java.util.*;
import java.io.*;
import org.apache.commons.math.distribution.*;
import org.apache.commons.math.special.Erf;

////////////////////////////////////////////////////////////////////////////////////////////

/* A NumericRandomVariable is both Discrete and Continuous.
 * Discrete probability at infinity is really concentrated there, but
 * Continuous probability at infinity is merely between infinity and the adjacent value(s).
 *
 * CONVENTION: (-x)^y => x^y * cos(pi*y), i.e. Re() part: (-1)^y = cos(pi*y) + i*sin(pi*y)
 *
 *	    // A^B cases (45-12)
 *          ///////////////////////////////////////////////////
 *	    //                                               //    X := {0,1,I,-I}~{1/4,..,1/4}
 *	    //  (B)\(A)__-I__-2__-1__-.5___0__.5__1__2_+I_   //   -I := -Infinity
 *	    //  +I  |   +-I +-I +-1    0   0   0  1 +I +I    //   +I := +Infinity
 *	    //  +R  |   +-I   R   R    R   0   R  R  R +I    //  +-x := {-x,x}~{.5,.5}
 *	    //   0  |     1   1   1    1   X   1  1  1  1    //    R := real value
 *	    //  -R  |     0   R   R    R +-I   R  R  R  0    //   -2 := (-I, -1) 
 *	    //  -I  |     0   0 +-1  +-I +-I   I  1  0  0    //   -1 := -1
 *          //                                               //   .5 := (0,1)
 *          ///////////////////////////////////////////////////   
 */

public class NumericRandomVariable {
    static private final String COMMA         = ",";
    static private final String COMMA_SPACE   = ", ";
    static private final String NEWLINE       = "\n";
    static private final String OPEN_PAREN    = "(";
    static private final String CLOSE_PAREN   = ")";
    static private final String OPEN_BRACKET  = "[";
    static private final String CLOSE_BRACKET = "]";
    static private final String OPEN_ANGLE    = "<";
    static private final String CLOSE_ANGLE   = ">";

    static private int PARTITION_LIMIT = 500;
    static private final double DTOL = 1.E-15;

    //////////////////////////////////////////////////////////////////////////////////////
    // Members

    DiscreteRandomVariable   m_discrete;
    ContinuousRandomVariable m_continuous;

    //////////////////////////////////////////////////////////////////////////////////////
    // Special

    static public void force_partition(double[] partition) {
	ContinuousRandomVariable.JointRV.force_partition(partition);
    }

    //////////////////////////////////////////////////////////////////////////////////////
    // Constructors

    public NumericRandomVariable() {
	m_discrete   = new DiscreteRandomVariable();
	m_continuous = new ContinuousRandomVariable();
    }
    public NumericRandomVariable(DiscreteRandomVariable X) {
	m_discrete   = new DiscreteRandomVariable(X);
	m_continuous = new ContinuousRandomVariable();
    }
    public NumericRandomVariable(ContinuousRandomVariable X) {
	m_discrete   = new DiscreteRandomVariable();
	m_continuous = new ContinuousRandomVariable(X);
    }
    public NumericRandomVariable(DiscreteRandomVariable X, ContinuousRandomVariable Y) {
	m_discrete   = new DiscreteRandomVariable(X);
	m_continuous = new ContinuousRandomVariable(Y);
    }
    public NumericRandomVariable(NumericRandomVariable X) {
	m_discrete   = new DiscreteRandomVariable  (X.m_discrete);
	m_continuous = new ContinuousRandomVariable(X.m_continuous);
    }
    public NumericRandomVariable(NumericRandomVariable X, NumericRandomVariable Y) {
	m_discrete   = DiscreteRandomVariable.  mix(X.m_discrete,   Y.m_discrete);
	m_continuous = ContinuousRandomVariable.mix(X.m_continuous, Y.m_continuous);
    }

    static public NumericRandomVariable createNoCopy(DiscreteRandomVariable D) {
	NumericRandomVariable R = new NumericRandomVariable();
	R.setDiscrete(D);
	return R;
    }
    static public NumericRandomVariable createNoCopy(ContinuousRandomVariable C) {
	NumericRandomVariable R = new NumericRandomVariable();
	R.setContinuous(C);
	return R;
    }
    static public NumericRandomVariable createNoCopy(DiscreteRandomVariable   D, 
						     ContinuousRandomVariable C) {
	NumericRandomVariable R = new NumericRandomVariable();
	R.setDiscrete(D);
	R.setContinuous(C);
	return R;
    }

    // Discrete constructors //////////////////////////////////

    static public NumericRandomVariable Discrete(double d) {
	double[] x = new double[1];
	double[] p = new double[1];
	x[0] = d;
	p[0] = 1;
	return NumericRandomVariable.Discrete(x, p);
    }

    static public NumericRandomVariable Discrete(double[] x, double[] p) {
	return new NumericRandomVariable(new DiscreteRandomVariable(x,p));
    }
    static public NumericRandomVariable Discrete(Double[] X, Double[] P) {
	double[] x = new double[X.length];
	double[] p = new double[P.length];
	for(int i = 0; i < X.length; i++) {x[i] = X[i]; p[i] = P[i];}
	return new NumericRandomVariable(new DiscreteRandomVariable(x,p));
    }

    public static NumericRandomVariable Poisson(double lambda) {
	return new NumericRandomVariable(DiscreteRandomVariable.Poisson(lambda));
    }
    
    // Continuous constructors ////////////////////////////////

    public static NumericRandomVariable Continuous(double[] x, double[] p) {
	return new NumericRandomVariable(new ContinuousRandomVariable(x,p));
    }
    static public NumericRandomVariable Continuous(Double[] X, Double[] P) {
	double[] x = new double[X.length];
	double[] p = new double[P.length];
	for(int i = 0; i < X.length; i++) {x[i] = X[i]; p[i] = P[i];}
	return new NumericRandomVariable(new ContinuousRandomVariable(x,p));
    }

    public static NumericRandomVariable Uniform(double a, double b) {
	double[] x = new double[2];
	double[] p = new double[1];
	x[0] = a;
	x[1] = b;
	p[0] = 1;
	return Continuous(x,p);
    }

    public static NumericRandomVariable Beta(double alpha, double beta, int N) {
	try {
	    return new NumericRandomVariable(ContinuousRandomVariable.Beta(alpha, beta, N));
	} catch(org.apache.commons.math.MathException me) {
	    System.out.println(me);
	    return null;
	}
    }

    public static NumericRandomVariable Cauchy(double median, double scale, int N) {
	try {
	    return new NumericRandomVariable(ContinuousRandomVariable.Cauchy(median, scale, N));
	} catch(org.apache.commons.math.MathException me) {
	    System.out.println(me);
	    return null;
	}
    }

    public static NumericRandomVariable ChiSquared(double df, int N) {
	try {
	    return new NumericRandomVariable(ContinuousRandomVariable.ChiSquared(df, N));
	} catch(org.apache.commons.math.MathException me) {
	    System.out.println(me);
	    return null;
	}
    }

    public static NumericRandomVariable Exponential(double lambda, int N) {
	try {
	    return new NumericRandomVariable(ContinuousRandomVariable.Exponential(lambda, N));
	} catch(org.apache.commons.math.MathException me) {
	    System.out.println(me);
	    return null;
	}
    }

    public static NumericRandomVariable F(double df1, double df2, int N) {
	try {
	    return new NumericRandomVariable(ContinuousRandomVariable.F(df1, df2, N));
	} catch(org.apache.commons.math.MathException me) {
	    System.out.println(me);
	    return null;
	}
    }

    public static NumericRandomVariable Normal(double mu, double sigma, int N) {
	if(sigma == 0) return NumericRandomVariable.Discrete(mu);
	try {
	    return new NumericRandomVariable(ContinuousRandomVariable.Normal(mu, sigma, N));
	}  catch(org.apache.commons.math.MathException me) {
	    System.out.println(me);
	    return null;
	}
    }

    public static NumericRandomVariable Normal(NumericRandomVariable Mu, double sigma, int N) {
	double[][] data = Mu.ForceDiscrete().Discrete().xp();

	double[] param  = data[0];
	double[] weight = data[1];
	int      M      = param.length;

	if(M == 0) return new NumericRandomVariable();

	NumericRandomVariable[] pdf = new NumericRandomVariable[M];
	for(int i = 0; i < M; i++) pdf[i] = Normal(param[i], sigma, N).scale(weight[i]);
	
	return mixN(pdf);
    }

    public static NumericRandomVariable Normal(double mu, NumericRandomVariable Sigma, int N) {
	double[][] data = Sigma.ForceDiscrete().Discrete().xp();

	double[] param  = data[0];
	double[] weight = data[1];
	int M = param.length;
	if(M == 0) return new NumericRandomVariable();

	NumericRandomVariable[] pdf = new NumericRandomVariable[M];
	for(int i = 0; i < M; i++) pdf[i] = Normal(mu, param[i], N).scale(weight[i]);
	
	return mixN(pdf);
    }

    public static NumericRandomVariable LogNormal(double mu, double sigma, int N) {
	try {
	    return new NumericRandomVariable(ContinuousRandomVariable.LogNormal(mu, sigma, N));
	} catch (org.apache.commons.math.MathException me) {
	    System.out.println(me);
	    return null;
	}
    }

    public static NumericRandomVariable Student(double median, double scale, double nu, int N) {
	try {
	    return new NumericRandomVariable(ContinuousRandomVariable.Student(median, scale, nu,N));
	}  catch(org.apache.commons.math.MathException me) {
	    System.out.println(me);
	    return null;
	}
    }

    public static NumericRandomVariable LevyStable(double alpha, double beta,
						   double location, double scale, int N) {
	double[][] xp = LevyStableDistribution.Density(alpha, beta, location, scale, N);
	return Continuous(xp[0], xp[1]);
    }

    public static NumericRandomVariable Triangular(double L, double M, double R,int N) {
	return new NumericRandomVariable(ContinuousRandomVariable.Triangular(L,M,R, N));
    }

    //////////////////////////////////////////////////////////////////////////////////
    // Access

    public DiscreteRandomVariable   Discrete()   {return m_discrete;}
    public ContinuousRandomVariable Continuous() {return m_continuous;}

    public NumericRandomVariable ForceDiscrete() { 
	DiscreteRandomVariable D = DiscreteRandomVariable.mix(Continuous().Discrete(),Discrete());;
	return new NumericRandomVariable(D);
    }

    public void setDiscrete  (DiscreteRandomVariable   D) {m_discrete   = D;}
    public void setContinuous(ContinuousRandomVariable C) {m_continuous = C;}

    public double[][] getData() {
	double[][] d = Discrete().xp();
	double[][] c = Continuous().xh_array(); 
	return new double[][]{d[0],d[1],c[0],c[1]};
    }

    public double[][] getCumulativeData() {
	// TODO: include discrete distribution
	double[][] xp = Continuous().xp_array(); 
	double[] x = xp[0];
	double[] p = xp[1];

	int N = x.length;
	double[] c = new double[N];
	double cumulative = 0;
	c[0] = 0;
	for(int i = 1; i < N; i++) {
	    cumulative += p[i-1];
	    c[i] = cumulative;
	}
	return new double[][]{x,c};
    }

    public NumericRandomVariable[] split(double x) {
	DiscreteRandomVariable[]   D = Discrete().  split(x);  // will bisect at x.
	ContinuousRandomVariable[] C = Continuous().split(x);
	NumericRandomVariable A = new NumericRandomVariable(D[0], C[0]);
	NumericRandomVariable B = new NumericRandomVariable(D[1], C[1]);
	return new NumericRandomVariable[]{A,B};
    }

    public NumericRandomVariable blunt(double x) {
	NumericRandomVariable[] LR = split(x);
	NumericRandomVariable   L  = LR[0];
	NumericRandomVariable   R  = LR[1];

	double p = L.P();

	DiscreteRandomVariable A = new DiscreteRandomVariable(new double[]{x}, new double[]{p});
	R.setDiscrete(DiscreteRandomVariable.mix(A, R.Discrete()));
	return R;
    }

    ////////////////////////////////////////////////////////////////////
    // Info
	
    public int size() {return Discrete().size() + Continuous().size();}

    public double P() {return Discrete().P() + Continuous().P();}

    public double LT(double x) {
	NumericRandomVariable X = split(x)[0];
	double px = EQ(x);
	if(px > 0) return X.P() - px/2;
	else       return X.P();
    }

    public double LTE(double x) {
	NumericRandomVariable X = split(x)[0];
	double px = EQ(x);
	if(px > 0) return X.P() + px/2;
	else       return X.P();
    }

    public double EQ(double x) {return Discrete().P(x);}

    public double GT(double x) {
	NumericRandomVariable X = split(x)[1];
	double px = EQ(x);
	if(px > 0) return X.P() - px/2;
	else       return X.P();
    }

    public double GTE(double x) {
	NumericRandomVariable X = split(x)[1];
	double px = EQ(x);
	if(px > 0) return X.P() + px/2;
	else       return X.P();
    }

    private static NumericRandomVariable CMP(NumericRandomVariable A, NumericRandomVariable B) {
	double[] partition = new double[]{Double.NEGATIVE_INFINITY, 0, Double.POSITIVE_INFINITY};
	NumericRandomVariable.force_partition(partition);
	NumericRandomVariable C = A.sub(B);
	NumericRandomVariable.force_partition(null);
	return C;
    }

    public double LT (NumericRandomVariable B) {return CMP(this,B).LT (0);}
    public double LTE(NumericRandomVariable B) {return CMP(this,B).LTE(0);}
    public double GT (NumericRandomVariable B) {return CMP(B,this).LT (0);}
    public double GTE(NumericRandomVariable B) {return CMP(B,this).LTE(0);}
    public double EQ (NumericRandomVariable B) {
	System.out.println("RandomVariable.EQ(A,B) not yet implemented");
	return 0;
    }

    // NumericRandomVariable
    public double E() {
	double d = Discrete().E();
	double c = Continuous().E();
	if(Double.isNaN(d) || Double.isNaN(c)) {
	    if(!Double.isNaN(d)) return d;
	    if(!Double.isNaN(c)) return c;
	    return Double.NaN;
	}
	return d+c;  // each side is understood to return fractional expectation
    }

    // NumericRandomVariable
    public double variance() {
	double mean = E();
	double d = Discrete().  variance(mean);
	double c = Continuous().variance(mean);
	if(Double.isNaN(d) || Double.isNaN(c)) {
	    if(!Double.isNaN(d)) return d;
	    if(!Double.isNaN(c)) return c;
	    return Double.NaN;
	}
	return d+c;  // each side is understood to return fractional variance
    }

    public double stdev()  {
	return Math.sqrt(variance());
    }

    public double median() {return percentile(0.5);} 

    public double percentile(double p) { 
	if(Discrete().size() == 0) return Continuous().Percentile(p);
	return Double.NEGATIVE_INFINITY;
    }

    ////////////////////////////////////////////////////////////////////
    // Operations

    public NumericRandomVariable scale(double d) {
	return createNoCopy(Discrete().scale(d), Continuous().scale(d));
    }
    static public NumericRandomVariable mix(NumericRandomVariable A, NumericRandomVariable B) {
	DiscreteRandomVariable   D = DiscreteRandomVariable.  mix(A.Discrete(),  B.Discrete());
	ContinuousRandomVariable C = ContinuousRandomVariable.mix(A.Continuous(),B.Continuous());
	return createNoCopy(D,C);
    }
    static public NumericRandomVariable mix2(NumericRandomVariable A, NumericRandomVariable B) {
	// Assume that A and B are will double the partition count which must be reduced.
	DiscreteRandomVariable   D = DiscreteRandomVariable.  mix (A.Discrete(),  B.Discrete());
	ContinuousRandomVariable C = ContinuousRandomVariable.mix2(A.Continuous(),B.Continuous());
	return createNoCopy(D,C);
    }
    static public NumericRandomVariable mixN(NumericRandomVariable[] A) {
	// The density at each point is effectively summed across each component. 
	// Assume that any weighting has already been applied.
	// The pairing polarity (left-to-right v. right-to-left) must alternate:
	//     1  2 3  4  5
	//      12   34   5
	//      12     345
	//        12345
	if(A.length == 0) return new NumericRandomVariable();
	return mixN(A, true);
    } 
    static private NumericRandomVariable mixN(NumericRandomVariable[] A, boolean LeftRight) {
	int N = A.length;
	if(N == 1) return A[0];
	//System.out.println("mixN("+N+")"); // DEBUG
	int half   = N / 2;
	int offset = N % 2;
	int M      = half + offset;
	NumericRandomVariable[] B = new NumericRandomVariable[M];
	if(LeftRight) {
	    for(int i = 0; i < half; i++) B[i]        = mix2(A[2*i], A[2*i+1]);
	    if(offset == 1)               B[M-1]      = A[N-1];
	} else {
	    for(int i = 0; i < half; i++) B[i+offset] = mix2(A[2*i+offset], A[2*i+1+offset]);
	    if(offset == 1)               B[0]        = A[0];
	}
	return mixN(B,!LeftRight);
    }

    ////////////////////////////////////////////////////////////////////
    // Arithmetic Operations

    public NumericRandomVariable neg() {
	return createNoCopy(Discrete().neg(), Continuous().neg());
    }
    public NumericRandomVariable reciprocal() {
	return createNoCopy(Discrete().reciprocal(), Continuous().reciprocal());
    }
    public NumericRandomVariable log() {
	return createNoCopy(Discrete().log(), Continuous().log());
    }
    public NumericRandomVariable exp() {
	return createNoCopy(Discrete().exp(), Continuous().exp());
    }
    public NumericRandomVariable sqrt() {
	return createNoCopy(Discrete().sqrt(), Continuous().sqrt());
    }
    public NumericRandomVariable pow(double d) {
	// ASSUME: d != +-I
	// Handle: d == 0, N/2
	DiscreteRandomVariable D = Discrete().pow(d);
	if(strict_equal(d,0)) {
	    double p = Continuous().P();
	    if(p == 0) return createNoCopy(D);
	    DiscreteRandomVariable E = new DiscreteRandomVariable();
	    E.head().insert(new Node(1,p));
	    return createNoCopy(DiscreteRandomVariable.mix(D, E));
	}

	if(half_fraction(d)) {
	    ContinuousRandomVariable[] neg_pos = Continuous().split(0);
	    ContinuousRandomVariable neg = neg_pos[0];
	    ContinuousRandomVariable pos = neg_pos[1];
	    ContinuousRandomVariable C   = pos.pow(d);
	    double p = neg.P();
	    if(p == 0) createNoCopy(D, C); // nothing to do.
	    DiscreteRandomVariable E = new DiscreteRandomVariable();
	    E.head().insert(new Node(0,p));
	    return createNoCopy(DiscreteRandomVariable.mix(D, E), C);
	}

	ContinuousRandomVariable C = Continuous().pow(d);

	return createNoCopy(D, C);
    }
    public NumericRandomVariable rpow(double d) {
	// ASSUME d != +-I
	// Handle d == 0, 1
	// RECALL: 0^+C = 0, 0^-C = +-I, 1^C = 1
	DiscreteRandomVariable D = Discrete().rpow(d);
	if(strict_equal(d,0)) {
	    ContinuousRandomVariable[] neg_pos = Continuous().split(0);
	    double neg = neg_pos[0].P();
	    double pos = neg_pos[1].P();
	    if(strict_equal(neg+pos,0)) return createNoCopy(D);
	    DiscreteRandomVariable E = new DiscreteRandomVariable();
	    if(!strict_equal(pos, 0)) E.head().insert(new Node(0,pos));
	    E.setNegInf(neg/2);
	    E.setPosInf(neg/2);
	    return createNoCopy(DiscreteRandomVariable.mix(D, E));
	}
	if(strict_equal(d,1)) {
	    double p = Continuous().P();
	    if(strict_equal(p,0)) return createNoCopy(D);
	    DiscreteRandomVariable E = new DiscreteRandomVariable();
	    E.head().insert(new Node(1, p));
	    return createNoCopy(DiscreteRandomVariable.mix(D, E));
	}

	// Handle the explicit INFINITY cases
	ContinuousRandomVariable C = Continuous().rpow(d);
	DiscreteRandomVariable   E = new DiscreteRandomVariable();
	double neg = Continuous().getNegInf();
	double pos = Continuous().getPosInf();
	if(d < -1) {
	    if(neg > 0) E.head().insert(new Node(0,neg));
	    if(pos > 0) {C.addNegInf(pos/2); C.addPosInf(pos/2);}
	} else if(d < 0) {
	    if(neg > 0) {C.addNegInf(neg/2); C.addPosInf(neg/2);}
	    if(pos > 0) E.head().insert(new Node(0,pos));
	} else if(d < 1) {
	    if(neg > 0) C.addPosInf(neg);
	    if(pos > 0) E.head().insert(new Node(0,pos));
	} else {
	    if(neg > 0) E.head().insert(new Node(0,neg));
	    if(pos > 0) C.addPosInf(pos);
	}

	return createNoCopy(DiscreteRandomVariable.mix(D, E), C);
    }
    public NumericRandomVariable add(double d) {
	return createNoCopy(Discrete().add(d), Continuous().add(d));
    }
    public NumericRandomVariable sub(double d) {
	return createNoCopy(Discrete().sub(d), Continuous().sub(d));
    }
    public NumericRandomVariable multiply(double d) {
	if(strict_equal(d,0)) {
	    DiscreteRandomVariable A = Discrete().multiply(d);
	    DiscreteRandomVariable B = new DiscreteRandomVariable();
	    B.head().insert(new Node(0,Continuous().coreP()));
	    B.setNegInf(Continuous().getNegInf());
	    B.setPosInf(Continuous().getPosInf());
	    return createNoCopy(DiscreteRandomVariable.mix(A,B));
	}
	return createNoCopy(Discrete().multiply(d), Continuous().multiply(d));
    }
    public NumericRandomVariable divide(double d) {
	return createNoCopy(Discrete().divide(d), Continuous().divide(d));
    }

    ////////////////////////////////////////////////////////////////////
    // Binary Arithmetic Operations

    public NumericRandomVariable add(NumericRandomVariable X2) {
	DiscreteRandomVariable   D = 
	    DiscreteRandomVariable.  add(Discrete(),   X2.Discrete());
	ContinuousRandomVariable CC = 
	    ContinuousRandomVariable.add(Continuous(), X2.Continuous());
	ContinuousRandomVariable CD = 
	    ContinuousRandomVariable.add(Continuous(), X2.Discrete());
	ContinuousRandomVariable DC = 
	    ContinuousRandomVariable.add(X2.Continuous(), Discrete());

	ContinuousRandomVariable CDC = ContinuousRandomVariable.mix(CD, DC);
	ContinuousRandomVariable C   = ContinuousRandomVariable.mix(CDC, CC);

	return createNoCopy(D, C);
    }

    public NumericRandomVariable sub(NumericRandomVariable X2) {
	DiscreteRandomVariable   D = 
	    DiscreteRandomVariable.  sub(Discrete(),   X2.Discrete());
	ContinuousRandomVariable CC = 
	    ContinuousRandomVariable.sub(Continuous(), X2.Continuous());
	ContinuousRandomVariable CD = 
	    ContinuousRandomVariable.sub(Continuous(), X2.Discrete());
	ContinuousRandomVariable DC = 
	    ContinuousRandomVariable.sub(X2.Continuous(), Discrete());

	ContinuousRandomVariable CDC = ContinuousRandomVariable.mix(CD, DC);
	ContinuousRandomVariable C   = ContinuousRandomVariable.mix(CDC, CC);

	return createNoCopy(D, C);
    }

    public NumericRandomVariable multiply(NumericRandomVariable B) {
	DiscreteRandomVariable   DD = 
	    DiscreteRandomVariable.  multiply(Discrete(),   B.Discrete());
	ContinuousRandomVariable CC = 
	    ContinuousRandomVariable.multiply(Continuous(), B.Continuous());

	NumericRandomVariable CD = multiply(Continuous(), B.Discrete());
	NumericRandomVariable DC = multiply(Discrete(),   B.Continuous());

	NumericRandomVariable pure  = createNoCopy(DD, CC);
	NumericRandomVariable mixed = mix(CD, DC);

	return mix(pure, mixed);
    }

    public static NumericRandomVariable multiply(ContinuousRandomVariable A, DiscreteRandomVariable B) {
	// A: -I..---..I
	// B: -I..-0-..I
	    
	DiscreteRandomVariable[] nzp = B.negative_zero_positive();
	
	DiscreteRandomVariable D = new DiscreteRandomVariable();

	double p = nzp[1].coreP();
	if(p > 0) {double q = A.P(); if(q > 0) D.head().insert(new Node(0, p*q));}

	ContinuousRandomVariable.OP_CASE op_case = ContinuousRandomVariable.OP_CASE.MULTIPLY;

	// The following do not handle INFINTY
	ContinuousRandomVariable negC = 
	    ContinuousRandomVariable.operate(A, nzp[0].neg(), op_case).neg();
	ContinuousRandomVariable posC = 
	    ContinuousRandomVariable.operate(A, nzp[2],       op_case);

	ContinuousRandomVariable C = ContinuousRandomVariable.mix(negC, posC);

	// Handle A probability at INFINITY
	double neg_factor = nzp[0].coreP();
	double pos_factor = nzp[2].coreP();
	C.setNegInf(A.getPosInf() * neg_factor + A.getNegInf() * pos_factor);
	C.setPosInf(A.getPosInf() * pos_factor + A.getNegInf() * neg_factor);
	
	// Handle B probability at INFINITY
	if(B.hasNegInf() || B.hasPosInf()) {
	    double[] np = A.neg_pos_P();
	    double   neg = B.getNegInf();
	    double   pos = B.getPosInf();
	    D.setNegInf(np[0] * pos + np[1] * neg);
	    D.setPosInf(np[0] * neg + np[1] * pos);
	}

	return createNoCopy(D, C);
    }

    public static NumericRandomVariable multiply(DiscreteRandomVariable A, ContinuousRandomVariable B) {
	return multiply(B, A);
    }
    
    public NumericRandomVariable divide(NumericRandomVariable B) {
	return multiply(B.reciprocal());
    }

    public NumericRandomVariable pow(NumericRandomVariable B) {

	DiscreteRandomVariable   DD = DiscreteRandomVariable.  pow(Discrete(),  B.Discrete());
	ContinuousRandomVariable CC = ContinuousRandomVariable.pow(Continuous(),B.Continuous());
	NumericRandomVariable    CD =                          pow(Continuous(),B.Discrete());
	NumericRandomVariable    DC =                          pow(Discrete(),  B.Continuous());

	NumericRandomVariable pure  = createNoCopy(DD, CC);
	NumericRandomVariable mixed = mix(CD, DC);

	return mix(pure, mixed);
    }

    public static NumericRandomVariable pow(ContinuousRandomVariable A, DiscreteRandomVariable B) {
	// Special cases: A^{0, +-I}. 
	// Notice: (-I,I)^0                   = 1
	//         (-I,-1)(-1,1)(1,I)^I       = (+-Ic)(  0 )(I)
	//         (-I,-1)(-1,0)(0,1)(1,I)^-I = (  0 )(+-Ic)(I)(0)

	DiscreteRandomVariable[] Bnzp = B.negative_zero_positive();
	DiscreteRandomVariable   Bnp  = DiscreteRandomVariable.mix(Bnzp[0], Bnzp[2]);

	ContinuousRandomVariable.OP_CASE op_pow  = ContinuousRandomVariable.OP_CASE. POW;

	ContinuousRandomVariable C = ContinuousRandomVariable.operate(A, Bnp, op_pow);
	DiscreteRandomVariable   D = new DiscreteRandomVariable();

	ContinuousRandomVariable[] A_minus_1    = A.split(-1);           // [-I,-1)(-1,I]
	ContinuousRandomVariable   A_lt_minus_1 = A_minus_1[0];          // [-I,-1)
	ContinuousRandomVariable   A_gt_minus_1 = A_minus_1[1];          // (-1, I]
	ContinuousRandomVariable[] A_plus_1     = A_gt_minus_1.split(1); // (-1, 1)(1,I]
	ContinuousRandomVariable   A_pm_1       = A_plus_1[0];           // (-1, 1)
	ContinuousRandomVariable   A_gt_1       = A_plus_1[1];           // ( 1, I]
	ContinuousRandomVariable[] A_m1_p1      = A_pm_1.split(0);       // (-1, 0)(0,1)
	ContinuousRandomVariable   A_m1         = A_m1_p1[0];            // (-1, 0)
	ContinuousRandomVariable   A_p1         = A_m1_p1[1];            // ( 0, 1)

	double one = Bnzp[1].P() * A.P();
	if(one > 0) D.head().insert(new Node(1, one));

	double zero = 0;
	// There may be a discrete zero encoded on last node of C from half-powers. Get it.
	if(C.last() != null && C.last().p() > 0) {
	    zero += C.last().p();
	    C.last().p(0);
	}
	zero += B.getPosInf() * A_pm_1.P();
	zero += B.getNegInf() * A_lt_minus_1.P();
	zero += B.getNegInf() * A_gt_1.P();
	if(zero > 0) D.head().insert(new Node(0, zero));
	
	double pos_inf = 0;
	pos_inf += A_gt_1.P()       * B.getPosInf();
	pos_inf += A_p1.P()         * B.getNegInf();

	double left_pos_inf   = A_lt_minus_1.P() * B.getPosInf();
	double middle_neg_inf = A_m1.P()         * B.getNegInf();
	C.addNegInf(left_pos_inf / 2 + middle_neg_inf / 2);
	C.addPosInf(left_pos_inf / 2 + middle_neg_inf / 2);

	return createNoCopy(D, C);
    }

    public static NumericRandomVariable pow(DiscreteRandomVariable A, ContinuousRandomVariable B) {
	// Special cases: {0, +-I}^B. 
	// Notice: 0^(-I,0)(0,I) = (+-I)( 0)
	//         I^(-I,0)(0,I) = (  0)( I)
	//        -I^(-I,0)(0,I) = (  0)(-I)

	DiscreteRandomVariable[] Anzp = A.negative_zero_positive();
	DiscreteRandomVariable   Anp  = DiscreteRandomVariable.mix(Anzp[0], Anzp[2]);

	ContinuousRandomVariable.OP_CASE op_rpow = ContinuousRandomVariable.OP_CASE.RPOW;

	ContinuousRandomVariable C = ContinuousRandomVariable.operate(B, Anp, op_rpow);
	DiscreteRandomVariable   D = new DiscreteRandomVariable();

	ContinuousRandomVariable[] B_neg_pos = B.split(0);   // [-I,0)(0,I]
	ContinuousRandomVariable   B_neg     = B_neg_pos[0]; // [-I,0)
	ContinuousRandomVariable   B_pos     = B_neg_pos[1]; // ( 0,I]

	double neg = B_neg.P();
	double pos = B_pos.P();

	double zero = 0;
	zero += Anzp[1].P()   * pos;
	zero += A.getPosInf() * neg;
	zero += A.getNegInf() * neg;
	if(zero > 0) D.head().insert(new Node(0, zero));

	double neg_inf = 0;
	neg_inf += Anzp[1].P()   * neg / 2;
	neg_inf += A.getNegInf() * pos;
	D.setNegInf(neg_inf);

	double pos_inf = 0;
	pos_inf += Anzp[1].P()   * neg / 2;
	pos_inf += A.getPosInf() * neg;
	D.setPosInf(pos_inf);

	return createNoCopy(D, C);
    }

    //////////////////////////////// Jython Support ///////////////////////////////

    public NumericRandomVariable __neg__()                             {return neg();}

    public NumericRandomVariable __add__ (double d)                    {return add(d);}
    public NumericRandomVariable __radd__(double d)                    {return add(d);}
    public NumericRandomVariable __add__ (NumericRandomVariable that)  {return add(that);}

    public NumericRandomVariable __sub__ (double d)                    {return sub(d);}
    public NumericRandomVariable __rsub__(double d)                    {return neg().add(d);}
    public NumericRandomVariable __sub__ (NumericRandomVariable that)  {return sub(that);}

    public NumericRandomVariable __mul__ (double d)                    {return multiply(d);}
    public NumericRandomVariable __rmul__(double d)                    {return multiply(d);}
    public NumericRandomVariable __mul__ (NumericRandomVariable that)  {return multiply(that);}

    public NumericRandomVariable __div__ (double d)                    {return divide(d);}
    public NumericRandomVariable __rdiv__(double d)                    {return reciprocal().multiply(d);}
    public NumericRandomVariable __div__ (NumericRandomVariable that)  {return divide(that);}

    public NumericRandomVariable __pow__ (double d)                    {return pow(d);}
    public NumericRandomVariable __rpow__(double d)                    {return rpow(d);}
    public NumericRandomVariable __pow__ (NumericRandomVariable that)  {return pow(that);}

    //////////////////////////////////////////////////////////////////////////////////
    // Display
    
    public String toString() {
	StringBuffer line = new StringBuffer();

	line.append(Discrete().toString());
	line.append(NEWLINE);
	line.append(Continuous().toString());

	return line.toString();
    }

    private static void write_line(Writer writer, double x, double y) throws IOException {
	writer.write(str(x) + "\t" + str(y) + "\n" );
    }

    private static void write_spike(Writer writer, Node a, double h) throws IOException {
	write_line(writer, a.x(), h);
	write_line(writer, a.x(), a.p());
	write_line(writer, a.x(), h);
    }

    private static double write_step(Writer writer, Node b, double h) throws IOException {
	if(!Double.isInfinite(b.next().x())) { 
	    //if(b.next().p() > 0) {
	    write_line(writer, b.x(), h);
	    double new_h = b.p() / (b.next().x() - b.x());
	    write_line(writer, b.x(), new_h);
	    return new_h;
	}
	write_line(writer, b.x(), h);
	write_line(writer, b.x(), 0);
	return 0;
    }

    public void plot(String filename) throws IOException  {
	Writer writer = new OutputStreamWriter(new FileOutputStream(filename));
	try {
	    Node a = Discrete  ().head().next();
	    Node b = Continuous().head().next();
	    
	    double h = 0;

	    while(a != Discrete().tail() || b != Continuous().tail()) {
		if(strict_equal(a.x(), b.x())) {
		    write_spike(writer, a, h);
		    h = write_step(writer, b, h);
		    a = a.next();
		    b = b.next();
		} else if(a.x() < b.x()) {
		    write_spike(writer, a, h);
		    a = a.next();
		} else {
		    h = write_step(writer, b, h);
		    b = b.next();
		}
	    }
	    
	} finally {writer.close();}
    }

    //////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////

    static public class DiscreteRandomVariable extends SimpleRandomVariable {

	/* Assume probability at each node is > 0 (else what's the point?)
	 */

	public static enum OP_CASE {ADD, MULTIPLY, POW};

	private static double PROBABILITY_TOL = 1.E-16;  // least probability node

	////////////////////////////////////////////////////////////////////
	// Constructors

	public DiscreteRandomVariable()                       {super();}
	public DiscreteRandomVariable(double[] x, double p[]) {super(x,p);}
	public DiscreteRandomVariable(DiscreteRandomVariable A) {super(A);}

	public static DiscreteRandomVariable Poisson(double lambda) {
	    // f(i; lambda) = lambda**i * exp(-lambda) / k!
	    // Notice: f(i) = lambda/i * f(i-1), f(0) = exp(-lambda)

	    DiscreteRandomVariable A = new DiscreteRandomVariable();
	    Node a = A.head();
	    
	    a = a.insert(new Node(0, Math.exp(-lambda)));  // f(0)
	    double p = 0;
	    int i = 1;
	    while(p < 1-DTOL && i < PARTITION_LIMIT) {
		p += a.p();
		a = a.insert(new Node(i, a.p() * lambda/i));
		i++;
	    }
	    
	    return A;
	}

	/////////////////////////////////////////////////////////////////////
	// Access

	public double P(double x) {
	    Node node = first();
	    if(node == null) return 0;
	    while(node != tail()) {
		if(node.x() == x) return node.p();
		node = node.next();
	    }
	    return 0;
	}

	// DiscreteRandomVariable
	public double E() {
	    // See notes: 3/1/2012 p.2
	    if(hasPosInf() || hasNegInf()) {
		if(getPosInf() > getNegInf()) return tail().x();
		else                          return head().x();
	    }

	    Node node = first();
	    if(node == null) return Double.NaN;    // will be caught at mixed level.

	    double E = 0;
	    while(node != tail()) {
		E += node.x() * node.p();
		node = node.next();
	    }

	    return E;
	}

	// DiscreteRandomVariable
	public double variance(double mean) {
	    // See notes: 3/1/2012 p.2
	    if(hasPosInf() || hasNegInf()) {
		if(getPosInf() > getNegInf()) return tail().x();
		else                          return head().x();
	    }

	    Node node = first();
	    if(node == null) return Double.NaN;    // will be caught at mixed level.

	    double V = 0;
	    while(node != tail()) {
		double cx = node.x() - mean;
		V += (cx*cx) * node.p();
		node = node.next();
	    }

	    return V;
	}

	/////////////////////////////////////////////////////////////////////
	// Operations

	public DiscreteRandomVariable scale(double f) {
	    // multiply each probability value by f
	    DiscreteRandomVariable B = new DiscreteRandomVariable();
	    B.setPosInf(getPosInf() * f);
	    B.setNegInf(getNegInf() * f);
	    Node a = first();
	    if(a == null) return B;
	    Node b = B.head();
	    while(a != tail()) {
		b = b.insert(new Node(a.x(), a.p() * f)); // action
		a = a.next();
	    } 
	    return B;
	}

	public DiscreteRandomVariable[] split(double x) {
	    // Will bisect the value that equals x giving each side a portion.
	    DiscreteRandomVariable A = new DiscreteRandomVariable();
	    DiscreteRandomVariable B = new DiscreteRandomVariable();
	    DiscreteRandomVariable[] result = new DiscreteRandomVariable[] {A, B};
	    
	    A.setNegInf(getNegInf());
	    B.setPosInf(getPosInf());

	    Node a    = A.head();
	    Node b    = B.head();
	    Node node = first();  // local place-holder
	    if(node == null) return result;

	    // Fill in A with known-good nodes
	    while(node != tail() && node.x() <= x) {
		a = a.insert(new Node(node));
		node = node.next();
	    }
	    
	    if(strict_equal(a.x(), x)) {
		b = b.insert(new Node(a));
		b.p(b.p()/2);
		a.p(a.p()/2);
	    }
	    
	    // Finish B with remaining nodes
	    while(node != tail()) {
		b = b.insert(new Node(node));
		node = node.next();
	    }
	    
	    return result;
	}

	public DiscreteRandomVariable[] negative_zero_positive() {
	    DiscreteRandomVariable[] neg_pos = split(0);
	    DiscreteRandomVariable negative = neg_pos[0];
	    DiscreteRandomVariable positive = neg_pos[1];
	    Node nb = negative.last();
	    Node pb = positive.first();

	    double p_zero = 0;
	    if(nb != null && nb.x() == 0) {p_zero += nb.p(); nb.remove();}
	    if(pb != null && pb.x() == 0) {p_zero += pb.p(); pb.remove();}

	    DiscreteRandomVariable zero = new DiscreteRandomVariable();
	    if(p_zero > 0) zero.head().insert(new Node(0, p_zero));
	    
	    return new DiscreteRandomVariable[]{negative, zero, positive};
	}

	public DiscreteRandomVariable[] small_one_big() {
	    // ASSUME non-negative
	    DiscreteRandomVariable[] small_big = split(1);
	    DiscreteRandomVariable small = small_big[0];
	    DiscreteRandomVariable big   = small_big[1];
	    Node last  = small.last();
	    Node first = big.first();

	    double p_one = 0;
	    if(last  != null && last.x()  == 1) {p_one += last.p();  last.remove();}
	    if(first != null && first.x() == 1) {p_one += first.p(); first.remove();}

	    DiscreteRandomVariable one = new DiscreteRandomVariable();
	    if(p_one > 0) one.head().insert(new Node(1, p_one));
	    
	    return new DiscreteRandomVariable[]{small, one, big};
	}

	/*
	public DiscreteRandomVariable[] odd_fractional_even() {
	    // Infinity is considered fractional.

	    DiscreteRandomVariable odd        = new DiscreteRandomVariable();
	    DiscreteRandomVariable fractional = new DiscreteRandomVariable();
	    DiscreteRandomVariable even       = new DiscreteRandomVariable();

	    fractional.setNegInf(getNegInf());
	    fractional.setPosInf(getPosInf());

	    Node n = first();
	    if(n == null) return new DiscreteRandomVariable[]{odd, fractional, even};

	    Node o = odd.head();
	    Node f = fractional.head();
	    Node e = even.head();

	    while(n != tail()) {
		int x = (int)n.x();
		if(strict_equal(x, n.x())) {
		    if(x % 2 == 0) e = e.insert(new Node(n));
		    else           o = o.insert(new Node(n));
		} else             f = f.insert(new Node(n));
		n = n.next();
	    }

	    return new DiscreteRandomVariable[]{odd, fractional, even};
	}
	*/

	public static DiscreteRandomVariable mix(DiscreteRandomVariable A, 
						  DiscreteRandomVariable B) {
	    // Notice both lists start and stop in sync (at +/-infinity)
	    // 

	    Node a = A.first();
	    Node b = B.first();

	    if(a == null) {
		DiscreteRandomVariable C = new DiscreteRandomVariable(B);
		C.addNegInf(A.getNegInf());
		C.addPosInf(A.getPosInf());
		return C;
	    }
	    if(b == null) {
		DiscreteRandomVariable C = new DiscreteRandomVariable(A);
		C.addNegInf(B.getNegInf());
		C.addPosInf(B.getPosInf());
		return C;
	    }
	    
	    DiscreteRandomVariable C = new DiscreteRandomVariable();
	    C.setNegInf(A.getNegInf() + B.getNegInf());
	    C.setPosInf(A.getPosInf() + B.getPosInf());
	    
	    Node c = C.head();
	    
	    // A:  -I......|---|-|---|----|......I
	    // B:  -I..|-|---|-----|-|-|-|--|-|..I
	    
	    while(!(a == A.tail() && b == B.tail())) {
		int cmp = compare(a,b);
		Node d = new Node();
		
		// decide where new node will live (x)
		if(cmp == -1) d.x(a.x());
		else          d.x(b.x());
		
		// apply (discrete) probability to fresh node (d)
		if(cmp == -1)     d.p(a.p());
		else if(cmp == 0) d.p(a.p() + b.p());
		else              d.p(b.p());
		
		// move our pointers along
		if(cmp == -1)      a = a.next();
		else if(cmp == 0) {a = a.next(); b = b.next();}
		else               b = b.next();
		
		// attach our fresh node.
		c = c.insert(d);
	    }

	    return C;
	}

	public DiscreteRandomVariable neg() {
	    DiscreteRandomVariable B = new DiscreteRandomVariable();
	    
	    B.setNegInf(getPosInf());
	    B.setPosInf(getNegInf());
	    
	    Node b = B.head();
	    Node a = first();
	    
	    if(a == null) return B;
	    
	    while(a != tail()) {
		b.insert(new Node(-a.x(), a.p()));  // action
		a = a.next();
	    }
	    
	    return B;
	}

	private DiscreteRandomVariable pos_reciprocal() {
	    // Assume all numbers are non-negative. 
	    DiscreteRandomVariable B = new DiscreteRandomVariable();
	    Node b = B.head();
	    Node a = first();
	    
	    if(hasPosInf() || hasNegInf()) b = b.insert(new Node(0, getPosInf()+getNegInf()));

	    if(a == null) return B;
	    
	    if(strict_equal(a.x(), 0)) {
		B.setPosInf(a.p());  
		a = a.next();
	    }
	    
	    while(a != tail()) {
		b.insert(new Node(1/a.x(), a.p())); // action. No zero's allowed
		a = a.next();
	    }
	    
	    return B;
	}

	private DiscreteRandomVariable neg_reciprocal() {
	    // Assume all numbers are non-positive. 
	    DiscreteRandomVariable B = new DiscreteRandomVariable();
	    Node b = B.head();
	    Node a = last();
	    
	    if(hasPosInf() || hasNegInf()) b.insert(new Node(0, getPosInf()+getNegInf()));
	    
	    if(a == null) return B;
	    
	    if(strict_equal(a.x(), 0)) {
		B.setNegInf(a.p());  
		a = a.prev();
	    }
	    
	    while(a != head()) {
		b = b.insert(new Node(1/a.x(), a.p())); // action. No zero's allowed
		a = a.prev();
	    }
	    
	    return B;
	}

	public DiscreteRandomVariable reciprocal() {
	    DiscreteRandomVariable[] neg_pos = split(0);
	    return mix(neg_pos[0].neg_reciprocal(), neg_pos[1].pos_reciprocal());
	}

	public DiscreteRandomVariable log() {
	    DiscreteRandomVariable B = new DiscreteRandomVariable();
	    B.setPosInf(getPosInf());

	    Node a = first();
	    if(a == null) return B;
	    if(a.x() < 0) return split(0)[1].log();
	    if(strict_equal(a.x(),0)) {
		B.setNegInf(a.p());
		a = a.next();
	    }
	    
	    Node b = B.head();
	    while(a != tail()) {
		b = b.insert(new Node(Math.log(a.x()), a.p())); // action
		a = a.next();
	    }
	    
	    return B;
	}

	public DiscreteRandomVariable exp() {
	    DiscreteRandomVariable B = new DiscreteRandomVariable();
	    B.setPosInf(getPosInf());
	    
	    Node b = B.head();
	    Node a = first();
	    
	    if(a == null) return B;
	    if(hasNegInf()) {
		b = b.insert(new Node());
		b.p(getNegInf());
	    }
	    
	    while(a != tail()) {
		b = b.insert(new Node(Math.exp(a.x()), a.p())); // action
		a = a.next();
	    }
	    
	    return B;
	}

	public DiscreteRandomVariable sqrt() {
	    return pow(.5);
	}

	private DiscreteRandomVariable pos_sq() {
	    // ASSUME List is only non-negative elements
	    DiscreteRandomVariable B = new DiscreteRandomVariable();
	    B.setPosInf(getPosInf());
	    
	    Node a = first();
	    if(a == null) return B;
	    
	    Node b = B.head();
	    while(a != tail()) {
		b = b.insert(new Node(a.x()*a.x() ,a.p())); // action
		a = a.next();
	    }
	    
	    return B;
	}

	public DiscreteRandomVariable sq() {
	    DiscreteRandomVariable[] pos_neg = split(0);
	    return mix(pos_neg[0].neg().pos_sq(),pos_neg[1].pos_sq());
	}

	public DiscreteRandomVariable pow(double d) {
	    if(d == 1)  return new DiscreteRandomVariable(this);
	    // We call up to the higher-level function rather than the usual call-down.
	    DiscreteRandomVariable B = new DiscreteRandomVariable();
	    B.head().insert(new Node(d, 1));
	    return pow(this, B);
	}

	public DiscreteRandomVariable rpow(double d) {
	    // We call up to the higher-level function rather than the usual call-down.
	    DiscreteRandomVariable A = new DiscreteRandomVariable();
	    A.head().insert(new Node(d,1));
	    return pow(A, this);
	}

	public DiscreteRandomVariable add(double d) {
	    DiscreteRandomVariable B = new DiscreteRandomVariable();
	    B.setPosInf(getPosInf());
	    B.setNegInf(getNegInf());
	
	    Node a = first();
	    if(a == null) return B;
	    
	    Node b = B.head();
	    while(a != tail()) {
		b = b.insert(new Node(a.x()+d, a.p())); // Action
		a = a.next();
	    }
	    
	    return B;
	}

	public DiscreteRandomVariable sub(double d) {
	    DiscreteRandomVariable B = new DiscreteRandomVariable();
	    B.setPosInf(getPosInf());
	    B.setNegInf(getNegInf());
	
	    Node a = first();
	    if(a == null) return B;
	    
	    Node b = B.head();
	    while(a != tail()) {
		b = b.insert(new Node(a.x()-d, a.p())); // Action
		a = a.next();
	    }
	    
	    return B;
	}

	public DiscreteRandomVariable pos_multiply(double d) {
	    // ASSUME d > 0
	    DiscreteRandomVariable B = new DiscreteRandomVariable();
	    B.setPosInf(getPosInf());
	    B.setNegInf(getNegInf());
	
	    Node a = first();
	    if(a == null) return B;
	    
	    Node b = B.head();
	    while(a != tail()) {
		b = b.insert(new Node(a.x()*d, a.p())); // Action
		a = a.next();
	    }
	    
	    return B;
	}

	public DiscreteRandomVariable zero_multiply() {
	    DiscreteRandomVariable B = new DiscreteRandomVariable();
	    B.setPosInf(getPosInf());
	    B.setNegInf(getNegInf());
	
	    Node a = first();
	    if(a == null) return B;
	    
	    double p = 0;
	    while(a != tail()) {
		p += a.p();
		a = a.next();
	    }

	    B.head().insert(new Node(0, p));
	    
	    return B;
	}

	public DiscreteRandomVariable neg_multiply(double d) {
	    // ASSUME: d < 0
	    DiscreteRandomVariable B = new DiscreteRandomVariable();
	    B.setPosInf(getNegInf());
	    B.setNegInf(getPosInf());
	
	    Node a = last();
	    if(a == null) return B;
	    
	    Node b = B.head();
	    while(a != head()) {
		b = b.insert(new Node(a.x()*d, a.p())); // Action
		a = a.prev();
	    }
	    
	    return B;
	}

	
	public DiscreteRandomVariable multiply(double d) {
	    if(strict_equal(d,0)) return zero_multiply();
	    if(d < 0) return neg_multiply(d);
	    return pos_multiply(d);
	}

	public DiscreteRandomVariable zero_divide() {
	    DiscreteRandomVariable B = new DiscreteRandomVariable();
	    B.setPosInf(getPosInf());
	    B.setNegInf(getNegInf());
	
	    Node a = first();
	    if(a == null) return B;
	    
	    double p = 0;
	    while(a != tail() && a.x() < 0) {
		p += a.p();
		a = a.next();
	    }

	    B.addNegInf(p);

	    if(strict_equal(a.x(),0)) {
		B.addNegInf(a.p()/2);
		B.addPosInf(a.p()/2);
	    }

	    p = 0;
	    while(a != tail()) {
		p += a.p();
		a = a.next();
	    }

	    B.addPosInf(p);
	    
	    return B;
	}

	public DiscreteRandomVariable divide(double d) {
	    if(strict_equal(d,0)) return zero_divide();
	    if(d < 0) return neg_multiply(1/d);
	    return pos_multiply(1/d);
	}

	///////////////////////////////////////////////////////////////////////
	// Arithmetic

	public static DiscreteRandomVariable add(DiscreteRandomVariable A,
						 DiscreteRandomVariable B) {

	    DiscreteRandomVariable C = operate(A, B, OP_CASE.ADD);

	    // TODO: Handle INFINITY

	    C.setPosInf(A.getPosInf() + B.getPosInf());
	    C.setNegInf(A.getNegInf() + B.getNegInf());

	    return C;
	}

	public static DiscreteRandomVariable sub(DiscreteRandomVariable A,
						 DiscreteRandomVariable B) {
	    return add(A, B.neg());
	}

	public static DiscreteRandomVariable multiply(DiscreteRandomVariable A,
						      DiscreteRandomVariable B) {

	    DiscreteRandomVariable C = operate(A, B, OP_CASE.MULTIPLY);

	    // TODO: Handle INFINITY

	    C.setPosInf(A.getPosInf() * B.getPosInf() + A.getNegInf() * B.getNegInf());
	    C.setNegInf(A.getPosInf() * B.getNegInf() + A.getNegInf() * B.getPosInf());

	    return C;
	}

	public static DiscreteRandomVariable divide(DiscreteRandomVariable A,
						    DiscreteRandomVariable B) {
	    return multiply(A, B.reciprocal());
	}

	public static DiscreteRandomVariable pow(DiscreteRandomVariable A,
						 DiscreteRandomVariable B) {

	    DiscreteRandomVariable[] Anzp = A.negative_zero_positive();
	    DiscreteRandomVariable[] Bnzp = B.negative_zero_positive();

	    DiscreteRandomVariable   An   = Anzp[0];
	    DiscreteRandomVariable   Az   = Anzp[1];  // special case
	    DiscreteRandomVariable   Ap   = Anzp[2];

	    DiscreteRandomVariable   Bn   = Bnzp[0];
	    DiscreteRandomVariable   Bz   = Bnzp[1];  // special case
	    DiscreteRandomVariable   Bp   = Bnzp[2];

	    DiscreteRandomVariable   Anp  = mix(An, Ap);
	    DiscreteRandomVariable   Bnp  = mix(Bn, Bp);
	    
	    DiscreteRandomVariable   Q    = operate(Anp, Bnp, OP_CASE.POW);

	    // Handle extreme cases (copy: see top of file)
	    ///////////////////////////////////////////////////
	    //                                               //    X := {0,1,I,-I}~{1/4,..,1/4}
	    //  (B) ___-I__-2__-1__-.5___0__.5__1__2_+I_ (A) //   -I := -Infinity
	    //  +I  | +-I +-I +-1    0   0   0  1 +I +I      //   +I := +Infinity
	    //  +R  | +-I   R   R    R   0   R  R  R +I      //  +-x := {-x,x}~{.5,.5}
	    //   0  |   1   1   1    1   X   1  1  1  1      //    R := real value
	    //  -R  |   0   R   R    R +-I   R  R  R  0      //   -2 := (-I, -1) 
	    //  -I  |   0   0 +-1  +-I +-I   I  1  0  0      //   -1 := -1
            //                                               //   .5 := (0,1)
            ///////////////////////////////////////////////////   

	    DiscreteRandomVariable[] Ans1b  = An.neg().small_one_big();
	    DiscreteRandomVariable[] Aps1b  = Ap.small_one_big();

	    double a_neg_inf = A.getNegInf();
	    double a_neg_big = Ans1b[2].coreP();
	    double a_neg_one = Ans1b[1].coreP();
	    double a_neg_sml = Ans1b[0].coreP();
	    double a_zero    = Az.P();
	    double a_pos_sml = Aps1b[0].coreP();
	    double a_pos_one = Aps1b[1].coreP();
	    double a_pos_big = Aps1b[2].coreP();
	    double a_pos_inf = A.getPosInf();

	    double b_pos_inf = B.getPosInf();
	    double b_pos     = Bp.coreP();
	    double b_zero    = Bz.P();
	    double b_neg     = Bn.coreP();
	    double b_neg_inf = B.getNegInf();

	    double neg_inf = 0;
	    double pos_inf = 0;
	    double zero    = 0;
	    double pos_one = 0;
	    double neg_one = 0;

	    neg_inf += b_pos_inf * a_neg_inf / 2;
	    neg_inf += b_pos_inf * a_neg_big / 2;
	    neg_inf += b_pos     * a_neg_inf / 2;
	    neg_inf += b_neg     * a_zero / 2;
	    neg_inf += b_neg_inf * a_zero / 2;
	    neg_inf += b_neg_inf * a_neg_sml / 2;
	    neg_inf += b_zero    * a_zero / 4;

	    neg_one += b_pos_inf * a_neg_one / 2;
	    neg_one += b_neg_inf * a_neg_one / 2;

	    zero    += b_pos_inf * (a_neg_sml + a_zero + a_pos_sml);
	    zero    += b_pos     * a_zero;
	    zero    += b_zero    * a_zero / 4;
	    zero    += b_neg     * (a_neg_inf + a_pos_inf);
	    zero    += b_neg_inf * (a_neg_inf + a_neg_big + a_pos_big + a_pos_inf);

	    pos_one += b_pos_inf * a_neg_one / 2;
	    pos_one += b_pos_inf * a_pos_one;
	    pos_one += b_zero    * (a_neg_inf + a_neg_big + a_neg_one + a_neg_sml);
	    pos_one += b_zero    * (a_pos_sml + a_pos_one + a_pos_big + a_pos_inf);
	    pos_one += b_zero    * a_zero / 4;
	    pos_one += b_neg_inf * a_pos_one;
	    pos_one += b_neg_inf * a_neg_one / 2;

	    pos_inf += b_pos_inf * a_neg_inf / 2;
	    pos_inf += b_pos_inf * a_neg_big / 2;
	    pos_inf += b_pos_inf * a_pos_inf;
	    pos_inf += b_pos_inf * a_pos_big;
	    pos_inf += b_pos     * a_neg_inf / 2;
	    pos_inf += b_pos     * a_pos_inf;
	    pos_inf += b_zero    * a_zero / 4;
	    pos_inf += b_neg     * a_zero / 2;
	    pos_inf += b_neg_inf * a_zero / 2;
	    pos_inf += b_neg_inf * a_pos_sml;
	    pos_inf += b_neg_inf * a_neg_sml / 2;

	    DiscreteRandomVariable Z = new DiscreteRandomVariable();
	    if(pos_one  > 0) Z.head().insert(new Node( 1, pos_one));
	    if(zero     > 0) Z.head().insert(new Node( 0, zero));
	    if(neg_one  > 0) Z.head().insert(new Node(-1, neg_one));
	    Z.setPosInf(pos_inf);
	    Z.setNegInf(neg_inf);

	    // DEBUG
	    //System.out.println("Q: "+Q);
	    //System.out.println("Z: "+Z);

	    return mix(Q, Z);
	}

	////////////////////////////////////////////////////////////////////////////

	private static DiscreteRandomVariable operate(DiscreteRandomVariable A,
						      DiscreteRandomVariable B, 
						      OP_CASE op_case) {
	    // Does NOT process INFINITY
	    
	    int N = A.size();
	    int M = B.size();

	    if(N == 0 || M == 0) return new DiscreteRandomVariable();

	    int K = N*M;
	    double[] Y = new double[K];
	    double[] P = new double[K];

	    DiscreteRandomVariable C = new DiscreteRandomVariable();
	    
	    Node a = A.first();
	    for(int i = 0; i < N; i++) {
		Node b = B.first();
		for(int j = 0; j < M; j++) {
		    switch(op_case) {
		    case ADD:
			Y[i*M + j] = a.x() + b.x();
			break;
		    case MULTIPLY:
			Y[i*M + j] = a.x() * b.x();
			break;
		    case POW:
			if(a.x() < 0) {
			    if(half_fraction(b.x()))           // capture -1/2 powers
				 Y[i*M + j] = 0;               // cos(-pi/2) not quite zero.
			    else Y[i*M + j] = Math.pow(-a.x(), b.x()) * Math.cos(Math.PI*b.x()); 
			} else   Y[i*M + j] = Math.pow( a.x(), b.x()); 
			break;
		    }

		    P[i*M + j] = a.p() * b.p();  // always
		    b = b.next();
		}
		a = a.next();
	    }

	    // This is not the most efficient way to do this, but it works.
	    int[] index = isort(Y);
	    Y = sort(Y, index);
	    P = sort(P, index);

	    // Has built-in threshold compression.
	    // Accepts duplicates
	    double p = 0;   // probability captured below threshold.
	    Node c = C.head();
	    c = c.insert(new Node(Y[0], P[0]));
	    for(int k = 1; k < K; k++) {
		if(P[k] < PROBABILITY_TOL) {p += P[k]; continue;}

		if(Double.isNaN(Y[k])) {
		    C.addNegInf(P[k]/2);
		    C.addPosInf(P[k]/2);
		} else if(Double.isInfinite(Y[k])) {
		    if(Y[k] < 0) C.addNegInf(P[k]);
		    else         C.addPosInf(P[k]);
		} else if(strict_equal(Y[k], Y[k-1])) {
		    c.p(c.p()+P[k]);
		} else {
		    c = c.insert(new Node(Y[k], P[k]));
		}

		// allocate accumulated probability fragments
		c.prev().p(c.prev().p()+p/2);
		c.p(c.p()+p/2);
		p = 0;
	    }

	    // allocate left-over probability frags.
	    if(p > PROBABILITY_TOL) 
		if(Double.isNaN(Y[K-1])) {
		    C.addNegInf(p/2);
		    C.addPosInf(p/2);
		} else if(Double.isInfinite(Y[K-1]))
		    C.addPosInf(p);
		else
		    c = c.insert(new Node(Y[K-1], p)); 

	    return C;
	}

	///////////////////////////////////////////////////////////////////////
	// Display

	private static void write_line(Writer writer, double x, double y) throws IOException {
	    writer.write(str(x) + "\t" + str(y) + "\n" );
	}

	public void plot(String filename) throws IOException  {
	    // see: 11/30/10p1
	    Writer writer = new OutputStreamWriter(new FileOutputStream(filename));
	    try {
		Node node = first();
		if(node == null) 
		    writer.write("0\t0\n");
		else {
		    while(node != tail()) {
			write_line(writer, node.x(), 0);
			write_line(writer, node.x(), node.p());
			write_line(writer, node.x(), 0);
			node = node.next();
		    }
		}
	    }
	    finally {
		writer.close();
	    }
	}
	
    } // end DiscreteRandomVariable

    //////////////////////////////////////////////////////////////////////////////////

    static public class ContinuousRandomVariable extends SimpleRandomVariable {

	/*  The last() node must have probability ZERO.
	 *  Probability of core (non-endpoint) nodes is between adjacent nodes.
	 */

	private static double PROBABILITY_TOL = 1.E-16;  // least probability node
	private static int    PARTITION_MIN   = 100;

	public static enum OP_CASE {ADD, MULTIPLY, POW, RPOW};

	/////////////////////////////////////////////////////////////////////////
	// Constructors

	public ContinuousRandomVariable()                           {super();}
	public ContinuousRandomVariable(double[] x, double p[])     {super(x,p);}
	public ContinuousRandomVariable(ContinuousRandomVariable A) {super(A);}
	
	// Distribution are rescaled rather than encoding probability at Infinity

	public static ContinuousRandomVariable Beta(double alpha, double beta, int N) 
	    throws org.apache.commons.math.MathException {
	    if(N <= 0) N = PARTITION_MIN; // default
	    
	    BetaDistributionImpl Beta = new BetaDistributionImpl(alpha, beta);

	    ContinuousRandomVariable A = new ContinuousRandomVariable();
	    Node a = A.head();
	    
	    // The probability is applied to last node, position to new node.
	    double n = N-1; // n is number of cells, N is number of endpoints.
	    a = a.insert(new Node(0));
	    double p = 0;
	    for(int i = 1; i < N; i++) {
		double z  = i/n;
		double p1 = Beta.cumulativeProbability(z);
		a.p(p1 - p);
		a = a.insert(new Node(z)); 
		p = p1;
	    }
	    
	    return A;
	}

	// OBSOLETE - using Student_z instead
	private static double cauchy_z(int i, int N, int L, double lambda, int M) {
	    // implements double-exponential sampling
	    // See: 12/15/2010p7

	    double n = i*2.0/(N-1)-1;

	    // linear sampling region
	    if(-lambda <= n && n <= lambda) return L/lambda * n;

	    double alpha = (M-Math.log(L))/(1-lambda);
	    
	    double z;
	    if(n > 0) z =  L * Math.exp(( n-lambda)*alpha);
	    else      z = -L * Math.exp((-n-lambda)*alpha);
	    
	    //System.out.println("i,n,z "+i+", "+n+", "+z); // DEBUG

	    return z;
	}

	// NB: Cauchy is Student(nu = 1)
	public static ContinuousRandomVariable Cauchy(double median, double scale, int N) 
	    throws org.apache.commons.math.MathException {
	    if(N <= 0) N = PARTITION_MIN; // default

	    // See: 12/15/2010 p.7 && 1/23/2012 p.5

	    double   knee  = 0.8;          // percentage points within quasi-linear sampling range
	    double   W     = 12;           // number of standard deviations of 1/2-width
	    if(N >= 500) W = 16;           // attempt to avoid too much probability in the tails.
	    double       L = 1E10;

	    CauchyDistributionImpl   F = new CauchyDistributionImpl();
	    ContinuousRandomVariable A = new ContinuousRandomVariable();

	    Node a = A.head();
	    
	    double min_point = Student_z(0, N, knee, W, L);
	    double max_point = Student_z(N, N, knee, W, L);
	    double c         = F.cumulativeProbability(min_point); // first cumulative probability
	    double d         = F.cumulativeProbability(max_point); // last cumulative probability
	    double MP        = (c + 1-d);                          // "missing" probability

	    a = a.insert(new Node(min_point*scale + median));
	    
	    // The probability is applied to last node, position to new node.
	    for(int i = 1; i < N; i++) {
		double z  = Student_z(i, N, knee, W, L);
		double c1 = F.cumulativeProbability(z);
		a.p((c1 - c)/(1-MP));          // scaled to include endpoint fragments
		a = a.insert(new Node(z*scale + median)); 
		c = c1;
	    }
	    
	    return A;
	}

	private static double chi_sqaured_z(int i, int N, double R, double df) {
	    // See: 12/15/2010p9
	    double b  = df*1.5/Math.log(df+1)+R;
	    double mu = Math.max(df-2,0);
	    double a  = Math.min(mu,b);
	    return i*(a+b)/(N-1)+mu-a;
	}

	public static ContinuousRandomVariable ChiSquared(double df, int N) 
	    throws org.apache.commons.math.MathException {
	    if(N <= 0) N = PARTITION_MIN; // default
	    
	    double       R = 16;          // standard width of sampling interval
	    if(N >= 500) R = 20;          // attempt to avoid too much probability in the tails.
	    
	    ChiSquaredDistributionImpl F = new ChiSquaredDistributionImpl(df); 
	    
	    ContinuousRandomVariable A = new ContinuousRandomVariable();
	    Node a = A.head();
	    
	    double z_max = chi_sqaured_z(N-1, N, R, df);
	    double f = 1/F.cumulativeProbability(z_max);
	    
	    a = a.insert(new Node(0));
	    double p = 0;
	    // The probability is applied to last node, position to new node.
	    for(int i = 1; i < N; i++) {
		double z  = chi_sqaured_z(i, N, R, df);
		double p1 = F.cumulativeProbability(z);
		a.p((p1 - p)*f);          // scaled to include endpoint fragments
		a = a.insert(new Node(z)); 
		p = p1;
	    }
	    
	    return A;
	}

	private static double exponential_z(int i, int N, double lambda) {
	    // See: 12/15/10p11
	    double last_z = Math.min(2*lambda*(lambda+2)+1,8*lambda);
	    return i * last_z / (N-1);
	}

	public static ContinuousRandomVariable Exponential(double lambda, int N) 
	    throws org.apache.commons.math.MathException {
	    if(N <= 0) N = PARTITION_MIN; // default
	    
	    double       R = 16;          // standard width of sampling interval
	    if(N >= 500) R = 20;          // attempt to avoid too much probability in the tails.
	    
	    ExponentialDistributionImpl F = new ExponentialDistributionImpl(lambda); 
	    
	    ContinuousRandomVariable A = new ContinuousRandomVariable();
	    Node a = A.head();
	    
	    double z_max = exponential_z(N-1, N, lambda);
	    double f = 1/F.cumulativeProbability(z_max);
	    
	    a = a.insert(new Node(0));
	    double p = 0;
	    // The probability is applied to last node, position to new node.
	    for(int i = 1; i < N; i++) {
		double z  = exponential_z(i, N, lambda);
		double p1 = F.cumulativeProbability(z);
		a.p((p1 - p)*f);          // scaled to include endpoint fragments
		a = a.insert(new Node(z)); 
		p = p1;
	    }
	    
	    return A;
	}

	private static double F_z(int i, int N, double df1, double df2) {
	    // See: 12/16/10p1..4
	    // Studied the necessary limits based on df1, df2 to find empirical formulas,
	    double Rmax   = 3.7;
	    double Rmin   = Math.max(1.2, 2-.8/6*Math.log(df1));
	    double R      = Math.max(Rmin, Math.min(Rmax, 2-2.0/3*(Math.log(df2)-3)));
	    double offset = -2;                          // start second exp early
	    double X      = Math.exp(i*R/(N-1))-1;
	    return Math.exp(X+offset)-Math.exp(offset);  // double-exponentially sampled
	}

	public static ContinuousRandomVariable F(double df1, double df2, int N) 
	    throws org.apache.commons.math.MathException {
	    if(N <= 0) N = PARTITION_MIN; // default
	    
	    FDistributionImpl F = new FDistributionImpl(df1, df2); 
	    
	    ContinuousRandomVariable A = new ContinuousRandomVariable();
	    Node a = A.head();
	    
	    double z_max = F_z(N-1,N,df1,df2);
	    double f = 1/F.cumulativeProbability(z_max);

	    // System.out.println("missed: "+(1-F.cumulativeProbability(z_max))); // DEBUG
	    
	    a = a.insert(new Node(0));
	    double p = 0;
	    // The probability is applied to last node, position to new node.
	    for(int i = 1; i < N; i++) {
		double z  = F_z(i,N,df1,df2);
		double p1 = F.cumulativeProbability(z);
		a.p((p1 - p)*f);          // scaled to include endpoint fragments
		a = a.insert(new Node(z)); 
		p = p1;
	    }

	    return A;
	}

	public static ContinuousRandomVariable Normal(double mu, double sigma, int N) 
	    throws org.apache.commons.math.MathException {
	    // We extend W standard deviations on either side of mean.
	    
	    if(N <= 0) N = PARTITION_MIN; // default
	    
	    double   W = 6;               // number of standard deviations of 1/2-width
	    if(N >= 500) W = 8;           // attempt to avoid too much probability in the tails.
	    double   v = sigma*sigma;     // variance
	    
	    NormalDistributionImpl Normal = new NormalDistributionImpl(); // Apache Math-Commons
	    
	    ContinuousRandomVariable A = new ContinuousRandomVariable();
	    Node a = A.head();
	    
	    double f =1/(1-(Normal.cumulativeProbability(-W)+1-Normal.cumulativeProbability(W)));
	    double c = Normal.cumulativeProbability(-W); // first cumulative probability
	    
	    // -I...|p0--|p1--|p2--..--|0...I
	    //      x0   x1   x2       xN-1
	    
	    a = a.insert(new Node((-W)*sigma + mu));
	    
	    // The probability is applied to last node, position to new node.
	    for(int i = 1; i < N; i++) {
		double z  = W*(2.0*i/(N-1)-1.0);
		double c1 = Normal.cumulativeProbability(z);
		a.p((c1 - c)*f);          // scaled to include endpoint fragments
		a = a.insert(new Node(z*sigma + mu)); 
		c = c1;
	    }
	    
	    return A;
	}

	public static double LogNormalCumulativeProbability(double x, double mu, double sigma2) 
	    throws org.apache.commons.math.MathException {
	    return 0.5*(1 + org.apache.commons.math.special.Erf.erf((Math.log(x) - mu)/(2*sigma2)));
	}

	public static ContinuousRandomVariable LogNormal(double mu, double sigma, int N) 
	    throws org.apache.commons.math.MathException {
	    // Cheap Version.
	    ContinuousRandomVariable normal = Normal(mu, sigma, N);
	    return normal.exp();
	}

	private static double Student_z_value(double z, double knee, double W, double L) {
	    // See notes: 1/23/2012 p.3
	    // knee is percent breakpoint for quasi-linear region
	    // w is the target value for end of quasi-linear region. Will actually be 2w. (why?)

	    double ki   = 1/knee;
	    double kiz  = ki*z;
	    double LWki = L/W - ki*ki;

	    if(LWki < 1) return 2*W*Math.pow(kiz,2);  // restrict to quasi-linear region only

	    double y = Math.log(LWki) / Math.log(ki);

	    return W*Math.pow(kiz,2) + W*Math.pow(kiz,y);
	}

	private static double Student_z(int i, int N, double knee, double W, double L) {
	    double z = i/(N/2.0) - 1;

	    if(z < 0) return -Student_z_value(-z, knee, W, L);
	    else      return  Student_z_value( z, knee, W, L);
	}

	public static ContinuousRandomVariable Student(double median, double scale, double nu, int N) 
	    throws org.apache.commons.math.MathException {
	    // We extend W standard deviations on either side of mean.
	    
	    if(N <= 0) N = PARTITION_MIN;  // default
	    
	    double   knee  = 0.8;          // percentage points within quasi-linear sampling range
	    double   W     = 12;           // number of standard deviations of 1/2-width
	    if(N >= 500) W = 16;           // attempt to avoid too much probability in the tails.

	    double            L = 10;
	    if(nu <= 1)       L = 1E10;
	    else if(nu <= 2)  L = 1E6;
	    else if(nu <= 3)  L = 1E4;
	    else if(nu <= 5)  L = 250;
	    else if(nu <= 6)  L = 100;
	    else if(nu <= 10) L = 25;
	    
	    TDistributionImpl T = new TDistributionImpl(nu); // Apache Math-Commons
	    
	    ContinuousRandomVariable A = new ContinuousRandomVariable();
	    Node a = A.head();
	    
	    double min_point = Student_z(0, N, knee, W, L);
	    double max_point = Student_z(N, N, knee, W, L);
	    double c         = T.cumulativeProbability(min_point); // first cumulative probability
	    double d         = T.cumulativeProbability(max_point); // last cumulative probability
	    double MP        = (c + 1-d);                          // "missing" probability

	    // System.out.println("Student-t missing probability: "+MP);  // DEBUG
	    
	    a = a.insert(new Node(min_point * scale + median));
	    
	    // The probability is applied to last node, position to new node.
	    for(int i = 1; i < N; i++) {
		double z  = Student_z(i, N, knee, W, L);
		double c1 = T.cumulativeProbability(z);
		a.p((c1 - c)/(1-MP));          // new probability scaled to include endpoint fragments
		a = a.insert(new Node(z * scale + median)); 
		c = c1;
	    }
	    
	    return A;
	}

	public static ContinuousRandomVariable Triangular(double L, double M, double R, int N) {
	    // Assuming N is even. and triangle(L,M,R)
	    int    N2  = N/2;
	    double dxL = (M-L)/N2;
	    double dxR = (R-M)/N2;
	    
	    ContinuousRandomVariable A = new ContinuousRandomVariable();
	    Node a = A.head();
	    for(int i = 0; i < N; i++) {
		if(i < N/2) {
		    double x = L+i*dxL;
		    double h = 2/(R-L)/(M-L)*(x+dxL/2-L); // midpoint
		    a = a.insert(new Node(x,h*dxL));
		} else {
		    double x = M+(i-N2)*dxR;
		    double h = 2/(R-L)/(R-M)*(R - (x+dxR/2)); // midpoint
		    a = a.insert(new Node(x,h*dxR));
		}
	    }
	    a.insert(new Node(R));
	    
	    return A;
	}

	/////////////////////////////////////////////////////////////////////
	// Access

	public DiscreteRandomVariable Discrete() {
	    // Concentrate interval probability at midpoint to form discrete distribution.
	    DiscreteRandomVariable D = new DiscreteRandomVariable();
	    D.setNegInf(getNegInf());
	    D.setPosInf(getPosInf());

	    if(head().next() == tail()) return D;
	    Node node = first();
	    Node d = D.head();
	    
	    while(node.next() != tail()) {
		d = d.insert(new Node((node.x() + node.next().x())/2,node.p()));
		node = node.next();
	    }

	    return D;
	}

	public double Median() {return Percentile(0.5);}

	public double Percentile(double p) {
	    Node node = head();
	    double P = P();
	    double Q = 0;
	    p = p * P;        // normalize. p is a fraction of the total probability.

	    while(Q+node.p() <= p && node != tail()) {
		Q += node.p();
		node = node.next();
	    }

	    // Handle end-cases
	    if(node == head()) return head().x();
	    if(node == tail()) {
		if(p <= Q) return node.prev().x();
		return tail().x();
	    }

	    //System.out.println("Q, (x-,p-,x+): "+Q+" ("+node.x()+", "+
	    //		       node.next().p()+", "+node.next().x()+")"); // DEBUG

	    double dq = node.p();
	    double dp = p - Q;
	    return ((1-dp/dq) * node.x() + dp/dq * node.next().x());
	}

	// ContinuousRandomVariable
	public double E() {
	    // See notes: 3/1/12 p.2
	    if(hasPosInf() || hasNegInf()) {
		if(getPosInf() > getNegInf()) return tail().x();
		else                          return head().x();
	    }
	    Node node = first();
	    if(node == null) return Double.NaN;

	    double E = 0;
	    double xL = node.x();
	    while(node.next() != tail()) {
		double xR = node.next().x();
		E += node.p() * (xL + xR) / 2;   // Assuming PU
		xL = xR;
		node = node.next();
	    }
	    
	    return E;
	}

	// ContinuousRandomVariable
	public double variance(double mean) {
	    if(hasPosInf() || hasNegInf()) {
		if(getPosInf() > getNegInf()) return tail().x();
		else                          return head().x();
	    }
	    Node node = first();
	    if(node == null) return Double.NaN;

	    double u = mean;
	    double V = 0;
	    double a = node.x();
	    while(node.next() != tail()) {
		double b = node.next().x();
		V += node.p() * ((a*a+b*b+a*b)/3 - u*(a+b-u));   // Assuming PU
		a = b;
		node = node.next();
	    }
	    
	    return V;
	}

	public double[][] xp_array() {
	    int N = size();
	    if(N <= 1) return null;     // No singletons!
	    double[] x = new double[N];
	    double[] p = new double[N];
	    Node node = first();
	    for(int i = 0; i < N-1; i++) {
		x[i] = node.x(); 
		p[i] = node.p();
		node = node.next();
	    }
	    x[N-1] = node.x();
	    p[N-1] = 0;
	    return new double[][]{x,p};
	}
	
	public double[] neg_pos_P() {
	    double p = getNegInf();
	    double q = getPosInf();

	    Node a = first();
	    if(a == null) return new double[]{p,q};

	    // handle singlton case
	    if(a == last()) {
		if(a.x() < 0) p += a.p();
		else          q += a.p();
		return new double[]{p,q};
	    }

	    while(a.next().x() < 0) {p += a.p(); a = a.next();}

	    if(a == last()) return new double[]{p,q};

	    if(a.p() > 0 && a.x() < 0) {
		double dx = a.next().x() - a.x();
		p += a.p() * (-a.x())     / dx;
		q += a.p() * a.next().x() / dx;
		a = a.next();
	    }
		
	    while(a != tail()) {q += a.p(); a = a.next();}

	    return new double[]{p,q};
	}
    
	/////////////////////////////////////////////////////////////////////
	// Operations

	public ContinuousRandomVariable scale(double f) {
	    // multiply each probability value by f
	    ContinuousRandomVariable B = new ContinuousRandomVariable();
	    B.setPosInf(getPosInf() * f);
	    B.setNegInf(getNegInf() * f);
	    Node a = first();
	    if(a == null) return B;
	    Node b = B.head();
	    while(a != tail()) {
		b = b.insert(new Node(a.x(), a.p() * f)); // action
		a = a.next();
	    } 
	    return B;
	}
	
	public ContinuousRandomVariable reduce(int n) {
	    // Skip n-1 nodes. Need to look ahead to figure out how to split d's
	    // see: 12/1/10p2
	    ContinuousRandomVariable B = new ContinuousRandomVariable();
	    B.setPosInf(getPosInf());
	    B.setNegInf(getNegInf());
	    
	    Node a = first();
	    if(a == null || a == last() && n < 2) return B;
	    
	    Node b = B.head().insert(new Node(a));
	    a = a.next();
	    
	    while(a != tail()) {
		// Push "end" out as far as you can (<= n and tail)
		double p   = 0;
		Node start = a;
		Node end   = start;
		for(int i = 0; i < n-1; i++) {
		    if(end == last()) break;
		    p += end.p();
		    end = end.next();
		}
		
		a = end;
		b = b.insert(new Node(a));
		a = a.next();
		b.prev().add_p(p);  // back-annotate the probability we accumulated
	    }
	    
	    return B;
	}

	public void reduceMiddle(int M) {
	    // See: notes 5/30/12 p.4 and 5/31/12 p.1 for special case
	    // given a target of N partition elements, remove middle partition elements
	    int N = size();
	    double factor = (double)M / (double)N;
	    Node a = first();

	    for(int i = 1; i < N-1; i++) {
		if(i > 1 && (int)((i-1)*factor) == (int)((i-2)*factor)) {
		    double x1  = a.prev().x();       double p1 = a.prev().p();
		    double x2  = a.x();              double p2 = a.p();
		    double x3  = a.next().x();       double p3 = a.next().p();
		    double x4  = a.next().next().x();
		    double h1  = p1 / (x2-x1); 
		    double h2  = p2 / (x3-x2); 
		    double h3  = p3 / (x4-x3);
		    double x12 = x1*x1; double x22 = x2*x2; 
		    double x32 = x3*x3; double x42 = x4*x4;
		    double P   = p1+p2+p3;
		    double E   = 0.5 * (h1*(x22-x12)+h2*(x32-x22)+h3*(x42-x32));
		    double y   = (x2+x3)/2.0;
		    double q1  = (P*(y+x4)-2*E)/(x4-x1);
		    double q2  = P - q1;
		    if(q1 < 0 || q2 < 0) { 
			//System.out.println("reduceMiddle("+M+") i: "+i+" q1:"+q1+" q2:"+q2);
			//System.out.println("x:"+x1+","+x2+","+x3+","+x4+"; y:"+y);
			//System.out.println("p:"+p1+","+p2+","+p3);
			double E_P = E/P;
			if(x1 < E_P && E_P < x4) {
			    if(E_P < y) {
				q1 = P;
				q2 = 0;
				y  = 2*E_P - x1;
			    } else {
				q1 = 0;
				q2 = P;
				y  = 2*E_P - x4;
			    }
			} else {
			    // TODO: Prove that this case can't happen.
			    System.out.println("reduceMiddle("+M+") punt."); // DEBUG
			    a.prev().p(P);
			    Node next = a.next().next();
			    a.next().remove();
			    a.remove();
			    a = next;
			    continue;
			}
		    }
		    a.prev().p(q1);
		    a.next().p(q2);
		    a.next().x(y);
		    Node next = a.next();
		    a.remove();
		    a = next;
		} else
		    a = a.next();
	    }
	}

	public static ContinuousRandomVariable mix(ContinuousRandomVariable A, 
						   ContinuousRandomVariable B) {
	    // Notice both lists start and stop in sync (at +/-infinity)

	    Node a = A.first();
	    Node b = B.first();

	    if(a == null) {
		ContinuousRandomVariable C = new ContinuousRandomVariable(B);
		C.addNegInf(A.getNegInf());
		C.addPosInf(A.getPosInf());
		return C;
	    }
	    if(b == null) {
		ContinuousRandomVariable C = new ContinuousRandomVariable(A);
		C.addNegInf(B.getNegInf());
		C.addPosInf(B.getPosInf());
		return C;
	    }
	    
	    ContinuousRandomVariable C = new ContinuousRandomVariable();
	    C.setNegInf(A.getNegInf() + B.getNegInf());
	    C.setPosInf(A.getPosInf() + B.getPosInf());

	    Node c = C.head();

	    // A:  -I......|---|-|---|----|......I
	    // B:  -I..|-|---|-----|-|-|-|--|-|..I
	    
	    // Need to get the first node in place. Fill in probability later.
	    int cmp0 = compare(a,b);
	    if(cmp0 == -1)     {c = c.insert(new Node(a.x())); a = a.next();}
	    else if(cmp0 == 0) {c = c.insert(new Node(a.x())); a = a.next(); b = b.next();}
	    else               {c = c.insert(new Node(b.x())); b = b.next();}
	    
	    while(!(a == A.tail() && b == B.tail())) {
		double da = a.x() - a.prev().x();
		double db = b.x() - b.prev().x();
		double dc;
		int cmp = compare(a,b);
		Node d = new Node();
		
		// decide where new node will live (x)
		if(cmp == -1) {dc = a.x() - c.x(); d.x(a.x());}
		else          {dc = b.x() - c.x(); d.x(b.x());}
		
		// now apply probability to recent node (c)
		c.p(0);
		if(a.prev().p() > 0) c.p(c.p() + a.prev().p() * dc/da);
		if(b.prev().p() > 0) c.p(c.p() + b.prev().p() * dc/db);
		
		// move our pointers along
		if(cmp == -1)      a = a.next();
		else if(cmp == 0) {a = a.next(); b = b.next();}
		else               b = b.next();
		
		// attach our fresh node.
		c = c.insert(d);
	    }

	    return C;
	}

	public static ContinuousRandomVariable mix2(ContinuousRandomVariable A, 
						    ContinuousRandomVariable B) {
	    // Assume A and B will double the partition count which must be reduced by half
	    ContinuousRandomVariable C = mix(A,B);
	    C = C.removeSlivers();
	    int N = Math.max(A.size(),B.size());
	    C.reduceMiddle(N);
	    return C;
	}

	public ContinuousRandomVariable removeSlivers() {
	    ContinuousRandomVariable C = new ContinuousRandomVariable();
	    C.setNegInf(getNegInf());
	    C.setPosInf(getPosInf());

	    Node a = first();
	    Node c = C.head();
	    
	    double x = a.x();
	    double p = a.p();

	    while(a.next() != tail()) {
		a = a.next();
		boolean WIDE     = a.x() - x > DTOL;
		boolean PROBABLE = p > DTOL;
		if(WIDE && PROBABLE) {
		    c = c.insert(new Node(x,p));
		    x = a.x();
		    p = a.p();
		} else {
		    p += a.p();
		}
	    }

	    c.insert(new Node(a)); // last node
	    return C;
	}


	public ContinuousRandomVariable[] split(double x) {
	    ContinuousRandomVariable A = new ContinuousRandomVariable();
	    ContinuousRandomVariable B = new ContinuousRandomVariable();
	    ContinuousRandomVariable[] result = new ContinuousRandomVariable[] {A, B};
	    
	    A.setNegInf(getNegInf());
	    B.setPosInf(getPosInf());

	    Node a    = A.head();
	    Node b    = B.head();
	    Node node = first();  // local place-holder
	    if(node == null) return result;

	    if(x <= first().x()) {
		if(A.hasNegInf()) a.insert(new Node(x));
		while(node != tail()) {
		    b = b.insert(new Node(node));
		    node = node.next();
		}
		return result;
	    }

	    if(x >= last().x()) {
		while(node != tail()) {
		    a = a.insert(new Node(node));
		    node = node.next();
		}
		if(B.hasPosInf()) b.insert(new Node(x));
		return result;
	    }

	    // ASSUME: x in (first.x(), last.x())
	    // Fill in A with known-good nodes
	    while(node != tail() && node.next().x() <= x) {
		a = a.insert(new Node(node));
		node = node.next();
	    }

	    if(strict_equal(node.x(),x)) {
		a = a.insert(new Node(x));               // cap-off
	    } else  {
		if(node.p() > 0) {
		    double dx = node.next().x() - node.x();
		    a = a.insert(new Node(node.x(), node.p() * (x - node.x())/dx));
		    b = b.insert(new Node(x, node.p() * (node.next().x() - x)/dx));
		    a = a.insert(new Node(x));           // cap-off
		} else {
		    a = a.insert(new Node(node.x()));    // cap-off
		}
		node = node.next();
	    }

	    //System.out.println("node: "+node.prev()+" A: "+A+" B: "+B); // DEBUG

	    // Finish B with remaining nodes
	    while(node != tail()) {
		b = b.insert(new Node(node));
		node = node.next();
	    }
	    
	    return result;
	}

	private ContinuousRandomVariable pos_reciprocal() {
	    // Assume all numbers are non-negative. 
	    // A: (x0,p0)->(x1,p1)->(x2,p2)->.....->(xn-1,pn-1)->(xn,0)->I
	    // B: (xn,pn-1)->(xn-1,pn-2)->...->(x2,p1)->(x1,p0)->(x0,0)->I

	    ContinuousRandomVariable B = new ContinuousRandomVariable();
	    Node b = B.head();
	    Node a = first();
	    
	    if(a == null) return B;
	    
	    if(hasPosInf()) b = b.insert(new Node(0, getPosInf()));

	    if(strict_equal(a.x(), 0)) {
		B.setPosInf(a.p());
		a = a.next();
	    }

	    double p = 0;                      // delayed action
	    while(a != tail()) {
		b.insert(new Node(1/a.x(),p)); // reciprocal -- no zero's allowed
		p = a.p();                     // delay
		a = a.next();
	    }
	    
	    return B;
	}

	private ContinuousRandomVariable neg_reciprocal() {
	    // Assume all numbers are non-positive. 
	    // -I..[)[)[)..0  -or- -I..[)[)[0

	    ContinuousRandomVariable B = new ContinuousRandomVariable();
	    Node b = B.head();
	    Node a = first();

	    if(a == null) return B;

	    if(strict_equal(last().x(), 0)) {
		B.setNegInf(last().prev().p());
		last().prev().p(0);
		last().remove();
	    }

	    double p = 0;             // delayed action
	    if(hasNegInf()) {
		b.insert(new Node());
		p = getNegInf();
	    }
	    
	    while(a != tail()) {
		b.insert(new Node(1/a.x(),p)); // reciprocal -- no zero's allowed
		p = a.p();                     // delay
		a = a.next();
	    }
	    
	    return B;
	}

	public ContinuousRandomVariable reciprocal() {
	    ContinuousRandomVariable[] neg_pos = split(0);
	    return mix(neg_pos[0].neg_reciprocal(), neg_pos[1].pos_reciprocal());
	}

	/////////////////////////////////////////////////////////////////////////
	// Arithmetic

	public ContinuousRandomVariable neg() {
	    ContinuousRandomVariable B = new ContinuousRandomVariable();
	    
	    B.setNegInf(getPosInf());
	    B.setPosInf(getNegInf());
	    
	    Node b = B.head();
	    Node a = first();
	    
	    if(a == null) return B;
	    
	    double p = 0;                      // delayed action
	    while(a != tail()) {
		b.insert(new Node(-a.x(), p)); // action
		p = a.p();                     // delay
		a = a.next();
	    }

	    return B;
	}

	public ContinuousRandomVariable log() {
	    ContinuousRandomVariable B = new ContinuousRandomVariable();
	    B.setPosInf(getPosInf());

	    Node a = first();
	    if(a == null) return B;
	    if(a.x() < 0) return split(0)[1].log();
	    if(strict_equal(a.x(),0)) {
		B.setNegInf(a.p());
		a = a.next();
	    }
	    
	    Node b = B.head();
	    while(a != tail()) {
		b = b.insert(new Node(Math.log(a.x()), a.p())); // action
		a = a.next();
	    }
	    
	    return B;
	}

	public ContinuousRandomVariable exp() {
	    ContinuousRandomVariable B = new ContinuousRandomVariable();
	    B.setPosInf(getPosInf());
	    
	    Node b = B.head();
	    Node a = first();
	    
	    if(a == null) return B;
	    if(hasNegInf()) {
		b = b.insert(new Node());
		b.p(getNegInf());
	    }
	    
	    while(a != tail()) {
		b = b.insert(new Node(Math.exp(a.x()), a.p())); // action
		a = a.next();
	    }
	    
	    return B;
	}

	public ContinuousRandomVariable add(double d) {
	    ContinuousRandomVariable B = new ContinuousRandomVariable();
	    B.setPosInf(getPosInf());
	    B.setNegInf(getNegInf());
	
	    Node a = first();
	    if(a == null) return B;
	
	    Node b = B.head();
	    while(a != tail()) {
		b = b.insert(new Node(a.x() + d, a.p())); // action
		a = a.next();
	    }
	    
	    return B;
	}

	public ContinuousRandomVariable sub(double d) {
	    ContinuousRandomVariable B = new ContinuousRandomVariable();
	    B.setPosInf(getPosInf());
	    B.setNegInf(getNegInf());
	    
	    Node a = first();
	    if(a == null) return B;
	    
	    Node b = B.head();
	    while(a != tail()) {
		b = b.insert(new Node(a.x() - d, a.p())); // action
		a = a.next();
	    }
	    
	    return B;
	}

	public ContinuousRandomVariable multiply(double d) {
	    // ASSUME d != 0
	    if(d < 0) return neg().multiply(-d);
	    
	    ContinuousRandomVariable B = new ContinuousRandomVariable();
	    B.setPosInf(getPosInf());
	    B.setNegInf(getNegInf());
	    
	    Node a = first();
	    if(a == null) return B;
	    
	    Node b = B.head();
	    while(a != tail()) {
		b = b.insert(new Node(a.x() * d, a.p())); // action
		a = a.next();
	    }
	    
	    return B;
	}

	public ContinuousRandomVariable zero_divide() {
	    ContinuousRandomVariable A = new ContinuousRandomVariable();
	    A.setNegInf(getNegInf());
	    A.setPosInf(getPosInf());
	    double[] np = neg_pos_P();
	    A.addNegInf(np[0]);
	    A.addPosInf(np[1]);
	    return A;
	}

	public ContinuousRandomVariable divide(double d) {
	    if(strict_equal(d,0)) return zero_divide();
	    return multiply(1/d);
	}

	public ContinuousRandomVariable sqrt() {return pow(0.5);}
	public ContinuousRandomVariable sq()   {return pow(2);}

	public ContinuousRandomVariable pow(double d) {
	    // ASSUME: d != 0, d != +/-INFINITY, 
	    // IF      d == N/2 THEN only non-negative portion survives, other is lost.
	    if(d < 0)   return reciprocal().pow(-d);
	    if(d ==  1) return new ContinuousRandomVariable(this);

	    // ASSUME d > 0
	    ContinuousRandomVariable[] neg_pos = split(0);
	    ContinuousRandomVariable   neg     = neg_pos[0];
	    ContinuousRandomVariable   pos     = neg_pos[1];

	    //System.out.println("neg: "+neg+" pos: "+pos); // DEBUG

	    Node p = pos.first();
	    if(p != null) 
		while(p != pos.tail()) {
		    p.x(Math.pow(p.x(),d));
		    p = p.next();
		}
	    
	    if(half_fraction(d)) { // capture N/2 powers and assume neg == 0
		pos.setNegInf(neg.getNegInf());  // won't be returning neg side at all
		return pos;
	    } 
		
	    double r = Math.cos(Math.PI * d); // Real Value: Re(x) = r*x
	    // ASSUME: r != 0 (no N/2 powers)
	    if(r > 0) neg = neg.neg(); // fold over (to +x axis)
	    Node n = neg.first();
	    if(n != null) 
		while(n != neg.tail()) {
		    if(r > 0) n.x(r * Math.pow( n.x(),d));  // already negated.
		    else      n.x(r * Math.pow(-n.x(),d));
		    n = n.next();
		}
	
	    return mix(neg, pos);
	}

	public ContinuousRandomVariable rpow(double d) {
	    // ASSUME: d != {0, +1, +-I}
	    // Does NOT handle INFINITY
	    if(d < 0) return rpow_neg(d);
	    if(d < 1) return rpow(1/d).reciprocal();

	    ContinuousRandomVariable B = new ContinuousRandomVariable();
	    Node b = B.head();

	    Node a = first();
	    if(a != null) {
		while(a != tail()) {
		    b = b.insert(new Node(Math.pow(d, a.x()), a.p()));
		    a = a.next();
		}
	    }

	    return B;
	}

	public ContinuousRandomVariable rpow_neg(double d) {
	    // Infinities and Special Cases NOT HANDLED except those internally generated
	    // ASSUME d < 0 && d != -I
	    // See: 12/13/10p2
	    // The problem here is the projected oscillation from cos(pi*c)
	    // This causes the entire set of partition cells to be chopped up and reassembled.
	    // First we create 2-endpoint CRV's based new cell widths

	    Node a = first();
	    if(a == null) return new ContinuousRandomVariable();

	    Node T = new Node(); // Hold all the rows by the tail
	    Node t = T;

	    double offset = rpow_neg_limits_offset(d);

	    double inf = 0;  // gather any exploded-cell probability here.

	    // build Tail-based list of 2-element CRV's
	    while(a != last()) {
		ContinuousRandomVariable row = new ContinuousRandomVariable();
		double[] limits = rpow_neg_limits(d, a.x(), a.next().x(), offset);
		if(Double.isInfinite(limits[0]) || Double.isInfinite(limits[1]))
		    inf = a.p();
		else {
		    row.head().insert(new Node(limits[1]));
		    row.head().insert(new Node(limits[0], a.p()));
		    t = t.prev(row.head());
		}
		a = a.next();
	    }

	    ContinuousRandomVariable C = mixer(T);
	    
	    if(inf > 0) {
		C.setPosInf(inf/2);
		C.setNegInf(inf/2);
	    }

	    int N = size();
	    int K = C.size();
	    if(K/N == 0) return C;
	    return C.reduce(K/N);
	}

	// Use this for external testing. Make sure it's above the "other" so jython finds it.
	public static double[] rpow_neg_limits(double d, double c1, double c2) {
	    double offset = rpow_neg_limits_offset(d);
	    return rpow_neg_limits(d, c1, c2, offset);
	}

	public static double rpow_neg_limits_offset(double d) {
	    return Math.atan(Math.log(-d)/Math.PI)/Math.PI; // in (0, 0.5)
	}

	public static double[] rpow_neg_limits(double d, double c1, double c2, double offset) {
	    // See: 12/13/10p3
	    // find [min_max(d^c1, d^c2)] where d < 0 && d != -I
	    // A typical cell will be oscillated by Re(.) and have new limits
	    // Find nearest lower integer
	    
	    if(d > -1) return rpow_neg_limits(1/d, -c2, -c1);  // d in (-1,0)

	    double r1 = NumericRandomVariable.pow(d, c1);
	    double r2 = NumericRandomVariable.pow(d, c2);

	    double x1 = c1 - offset;
	    double x2 = c2 - offset;

	    int n1 = right_integer(x1);
	    int n2 = right_integer(x2);

	    //System.out.println("n1,n2, r1,r2, offset: "+n1+","+n2+", "+r1+","+r2+", "+offset); 

	    double max = Math.max(r1, r2);
	    double min = Math.min(r1, r2);

	    if(n2 - n1 == 0) return new double[]{min, max};

	    if(n2 - n1 == 1) 
		if(NumericRandomVariable.isEven(n2))
		    return new double[]{NumericRandomVariable.pow(d,n2-1+offset), max};
		else 
		    return new double[]{min, NumericRandomVariable.pow(d, n2-1+offset)};
	    
	    if(NumericRandomVariable.isEven(n2)) {
		max = Math.max(r2, NumericRandomVariable.pow(d,n2-2+offset));
		min = Math.min(r2, NumericRandomVariable.pow(d,n2-1+offset));
	    } else {
		max = Math.max(r2, NumericRandomVariable.pow(d,n2-1+offset));
		min = Math.min(r2, NumericRandomVariable.pow(d,n2-2+offset));
	    }

	    return new double[]{min, max};
	}

	///////////////////////////////////////////////////////////////////////
	// Binary Arithmentic

	public static ContinuousRandomVariable add(ContinuousRandomVariable A,
						   DiscreteRandomVariable B) {
	    return operate(A, B, OP_CASE.ADD);
	}

	public static ContinuousRandomVariable sub(ContinuousRandomVariable A,
						   DiscreteRandomVariable B) {
	    return operate(A, B.neg(), OP_CASE.ADD);
	}

	public static ContinuousRandomVariable sub(DiscreteRandomVariable A,
						   ContinuousRandomVariable B) {
	    return operate(B.neg(), A, OP_CASE.ADD); 
	}

	////////////////////////////////////////////////////////////////////////////

	private static ContinuousRandomVariable operate(ContinuousRandomVariable A,
							DiscreteRandomVariable B, 
							OP_CASE op_case) {

	    // DOES NOT handle INFINITY or discrete effects such as "times zero"
	    // Create #B scaled and altered versions of A. Mix them together.
	    // ASSUME order preserved within rows.
	    // ASSUME no duplicates within rows.

	    if(A.first() == null || B.first() == null) 
		return new ContinuousRandomVariable();

	    Node T = new Node(); // Hold all the rows by the tail
	    Node t = T;

	    Node b = B.first();

	    double neg_inf = 0;
	    double pos_inf = 0;
	    double zero    = 0;  // Discrete zero probability. Encoded onto last node!

	    while(b != B.tail()) {
		switch(op_case) {
		case ADD:      
		    ContinuousRandomVariable add_row = A.add(b.x()).scale(b.p());  
		    t = t.prev(add_row.head());
		    break;
		case MULTIPLY: 
		    ContinuousRandomVariable mult_row = A.multiply(b.x()).scale(b.p());
		    t = t.prev(mult_row.head());
		    break;
		case POW:
		    // must handle half-powers. 
		    if(half_fraction(b.x())) {
			ContinuousRandomVariable[] neg_pos = A.split(0);
			ContinuousRandomVariable   neg     = neg_pos[0];
			ContinuousRandomVariable   pos     = neg_pos[1];
			ContinuousRandomVariable   pow_row = pos.pow(b.x()).scale(b.p());
			zero += neg.P() * b.p();   // This is the bugaboo.
			pos_inf += pow_row.getPosInf();
			t = t.prev(pow_row.head());
		    } else {
			ContinuousRandomVariable pow_row = A.pow(b.x()).scale(b.p());
			neg_inf += pow_row.getNegInf();
			pos_inf += pow_row.getPosInf();
			t = t.prev(pow_row.head());
		    }
		    break;
		case RPOW:
		    ContinuousRandomVariable rpow_row = A.rpow(b.x()).scale(b.p());
		    neg_inf += rpow_row.getNegInf();
		    pos_inf += rpow_row.getPosInf();
		    t = t.prev(rpow_row.head());
		    break;
		}
		b = b.next();
	    }

	    ContinuousRandomVariable C = mixer(T);
	    C.setNegInf(neg_inf);
	    C.setPosInf(pos_inf);

	    int N = A.size();
	    int K = C.size();
	    if(K/N != 0) C = C.reduce(K/N);

	    if(zero > 0 && C.last() != null) C.last().p(zero);  // encode discrete zero, if any

	    return C;
	}

       	private static ContinuousRandomVariable mixer(Node T) {
	    // ASSUME: T is tail-based list of CRV's to be mixed.

	    ContinuousRandomVariable C = new ContinuousRandomVariable();
	    Node c = C.head();
	    Node t = T;
	    while(true) {
		// (1) scan for next y-value
		double ymin = Double.POSITIVE_INFINITY;
		t = T.prev();
		while(t != null) {
		    if(t.next().x() < ymin) ymin = t.next().x();
		    t = t.prev();
		}
		if(Double.isInfinite(ymin) && ymin > 0) break;  // EXIT POINT
		c = c.insert(new Node(ymin));

		// (2) assign probability to c.prev() node and trim node.x()==ymin nodes
		if(c != C.head()) {
		    t = T.prev();
		    double dy = c.x() - c.prev().x();
		    while(t != null) {
			if(t.p() > 0) 
			    c.prev().add_p( t.p()*(dy/(t.next().x()-t.x())) );
			if(t.next().x() <= ymin) {
			    t.x(t.next().x());
			    t.p(t.next().p());
			    t.next().remove();
			}
			t = t.prev(); 
		    }
		    
		    // sweep away the flitty stuff
		    if(c.prev().p() > 0 && c.prev().p() < PROBABILITY_TOL) c.prev().remove();  
		}
	    }

	    return C;
	}

	///////////////////////////////////////////////////////////////////
	// Continuous * Continuous Arithmetic

    private static ContinuousRandomVariable operate(ContinuousRandomVariable X1, 
						    ContinuousRandomVariable X2, 
						    JointRV.OP_CASE op_case) {
	JointRV joint = new JointRV(X1, X2, op_case);
	return joint.getY();
    }

    public static ContinuousRandomVariable add(ContinuousRandomVariable X1, 
					       ContinuousRandomVariable X2) {
	ContinuousRandomVariable Y = operate(X1, X2, JointRV.OP_CASE.ADD);

	// Rule: near-infinite addition.
	// Y++ = (X1++)(X2 - X2--/2) + (X2++)(X1 - X1--/2)
	// Y-- = (X1--)(X2 - X2++/2) + (X2--)(X1 - X1++/2)

	double PX1   = X1.P();
	double PX2   = X2.P();
	double PX1pp = X1.getPosInf();
	double PX1nn = X1.getNegInf();
	double PX2pp = X2.getPosInf();
	double PX2nn = X2.getNegInf();

	Y.setPosInf(PX1pp * (PX2 - PX2nn/2) + PX2pp * (PX1 - PX1nn/2));
	Y.setNegInf(PX1nn * (PX2 - PX2pp/2) + PX2nn * (PX1 - PX1pp/2));

	return Y;
    }

    public static ContinuousRandomVariable sub(ContinuousRandomVariable X1, 
					       ContinuousRandomVariable X2) {
	return add(X1, X2.neg());
    }

    public static ContinuousRandomVariable multiply(ContinuousRandomVariable X1, 
						    ContinuousRandomVariable X2) {
	ContinuousRandomVariable[] X1pair = X1.split(0);
	ContinuousRandomVariable[] X2pair = X2.split(0);

	ContinuousRandomVariable X1n = X1pair[0].neg();
	ContinuousRandomVariable X1p = X1pair[1];
	ContinuousRandomVariable X2n = X2pair[0].neg();
	ContinuousRandomVariable X2p = X2pair[1];

	ContinuousRandomVariable Ypp = operate(X1p, X2p, JointRV.OP_CASE.MULTIPLY);
	ContinuousRandomVariable Ynp = operate(X1n, X2p, JointRV.OP_CASE.MULTIPLY);
	ContinuousRandomVariable Ypn = operate(X1p, X2n, JointRV.OP_CASE.MULTIPLY);
	ContinuousRandomVariable Ynn = operate(X1n, X2n, JointRV.OP_CASE.MULTIPLY);

	ContinuousRandomVariable Yp = mix(Ypp, Ynn);
	ContinuousRandomVariable Yn = mix(Ynp, Ypn).neg();

	ContinuousRandomVariable Y = mix(Yn, Yp);

	// Rule: handle near-infinite multiplication
	// Y++ = (X1>0)(X2++) + (X1++)(X2+) + (X1<0)(X2--) + (X1--)(X2-)
	// Y-- = (X1>0)(X2--) + (X1++)(X2-) + (X1<0)(X2++) + (X1--)(X2+)

	double PX1p  = X1p.P();
	double PX1n  = X1n.P();
	double PX1pp = X1.getPosInf();
	double PX1nn = X1.getNegInf();
	double PX2p  = X2p.P();
	double PX2n  = X2n.P();
	double PX2pp = X2.getPosInf();
	double PX2nn = X2.getNegInf();

	Y.setPosInf(PX1p * PX2pp + PX1pp * PX2p + PX1n * PX2nn + PX1nn * PX2n);
	Y.setNegInf(PX1p * PX2nn + PX1pp * PX2n + PX1n * PX2pp + PX1nn * PX2p);

	return Y;
    }

    public static ContinuousRandomVariable divide(ContinuousRandomVariable X1, 
						  ContinuousRandomVariable X2) {
	return multiply(X1, X2.reciprocal());
    }

    public static ContinuousRandomVariable pow(ContinuousRandomVariable X1, 
					       ContinuousRandomVariable X2) {
	
	// NOTE: We only handle (+C)^C and simply ignore -C.

	ContinuousRandomVariable[] X1pair = X1.split(0)[1].split(1);
	ContinuousRandomVariable[] X2pair = X2.split(0);

	ContinuousRandomVariable X10 = X1pair[0].reciprocal(); // [0,1]
	ContinuousRandomVariable X11 = X1pair[1];              // [1,inf]
	ContinuousRandomVariable X2n = X2pair[0].neg();        // [-inf,0]
	ContinuousRandomVariable X2p = X2pair[1];              // [0,inf]

	// Notice: neg exponent -> reciprocal of result.
	ContinuousRandomVariable Y0p = operate(X10, X2p, JointRV.OP_CASE.POW).reciprocal();
	ContinuousRandomVariable Y0n = operate(X10, X2n, JointRV.OP_CASE.POW);// dbl reciprocal.
	ContinuousRandomVariable Y1p = operate(X11, X2p, JointRV.OP_CASE.POW);
	ContinuousRandomVariable Y1n = operate(X11, X2n, JointRV.OP_CASE.POW).reciprocal(); 

	ContinuousRandomVariable Y = mix(mix(Y0p,Y1p), mix(Y0n, Y1n));

	// TODO: apply infinity rules
	return Y;
    }

    static public class Ei {
	//Euler-Masceroni Constant:
	static private final double gamma = 0.57721566490153286060651209008240243104215933593992;

	static public double getEi(double x) {
	    int n = getNumTerms(x);
	    return getEi(x, n);
	}

	static public double getEiDiff(double a, double b) {
	    int n = Math.max(getNumTerms(a), getNumTerms(b));
	    return Math.log(b/a) + getSeriesSum(b,n) - getSeriesSum(a,n);
	}

	static public double getEi(double x, int n) {
	    return gamma + Math.log(x) + getSeriesSum(x,n);
	}

	static int getNumTerms(double x) {
	    // empirically determined envelope curve for 1 to 100
	    if(x < 1.E-8) return 1;     // we can get away with 1 (but never 0) terms
	    int N = (int)(1.5625*(1+1./(5*x+5))*(x-99)+189);
	    return Math.min(500, N);   // Don't get carried away.
	}

	static double getSeriesSum(double x, int n) {
	    double v = 1./n;
	    for(int k = n; k > 1; k--) 
		v = (x * v) / k +  1./(k-1);  // careful with calculation order
	    return x * v;
	}

    } // end Ei


	public static class JointRV {
	    public  static enum OP_CASE    {ADD, MULTIPLY, POW};
	    private static enum BLOCK_CASE {BELOW, ABOVE};
	    public  static double[]  FORCE_PARTITION = null; // force a particular partition
	    
	    OP_CASE   m_op_case;       // operation case
	    
	    double   m_pos_inf_x1;
	    double   m_neg_inf_x1;
	    double   m_pos_inf_x2;
	    double   m_neg_inf_x2;
	    
	    double[] m_x1;            // working support of X1 (N)
	    double[] m_x2;            // support of X2         (M)
	    
	    double[][] m_y;           // y-values back-annotated to lattice (N x M)
	    double[][] m_p = null;    // Continuous * Continuous probability density (N-1 x M-1)
	    
	    // if this remains null, we return empty except for INFINITY
	    double[] m_partition = null;  // the working support of our Y = X1*X2 rv 
	    
	    int[]    m_left;          // cell search reduction. Reusable. (M-1)
	    int[]    m_right;         // cell search reduction. Reusable. (M-1)

	    /////////////////////////////////////////////////////////////////////

	    static public void force_partition(double[] partition) {FORCE_PARTITION = partition;}
	    
	    public JointRV(ContinuousRandomVariable X1, 
			   ContinuousRandomVariable X2, 
			   OP_CASE op_case) {

		m_op_case = op_case;
		
		m_pos_inf_x1 = X1.getPosInf();
		m_neg_inf_x1 = X1.getNegInf();
		m_pos_inf_x2 = X2.getPosInf();
		m_neg_inf_x2 = X2.getNegInf();
		
		double[][] X1xp = X1.xp_array();
		double[][] X2xp = X2.xp_array();
		
		if(X1xp == null || X2xp == null) return;
		
		double[] x1 = X1xp[0];
		double[] x2 = X2xp[0];
		double[] p1 = X1xp[1];
		double[] p2 = X2xp[1];
		
		int N = x1.length;
		int M = x2.length;
		
		if(op_case == OP_CASE.POW) m_x1 = log(x1);
		else                       m_x1 = x1;
		m_x2 = x2;
		
		m_y = create_y(m_x1, m_x2, op_case);
		m_p = create_p(p1,p2,x1,x2);
		
		m_left  = new int[M-1];   // recyclable
		m_right = new int[M-1];   // recyclable
		
		int L = 2 * Math.max(N, M);                                // 2x "oversample"
		L = Math.max(PARTITION_MIN, Math.min(L, PARTITION_LIMIT)); // the govenor.
		m_partition = create_partition(x1, x2, L, op_case);
		if(FORCE_PARTITION != null) m_partition = FORCE_PARTITION; // override
	    }
	    
	    //////////////////////////////////////////////////////////////////////
	    
	    public ContinuousRandomVariable getY() {
		if(m_partition == null) return new ContinuousRandomVariable();

		int      L = m_partition.length;
		double[] y = m_partition;
		double[] p = new double[L-1];
		
		zero(p);
		
		// fill the probabilities into our new y-array.
		for(int k = 0; k < L-1; k++) {
		    // See: 11/18/2010p1 and 12/2/2010p1
		    
		    double yL = y[k];
		    double yR = y[k+1];
		    
		    setSearchLimits(m_y, yL, yR, m_left, m_right);
		    
		    int      M = m_x2.length;
		    double x10; double x11; double x20; double x21; // reusable.
		    
		    // walk the involved cells. Will unconditionally check each row
		    for(int j = 0; j < M-1; j++) 
			for(int i = m_left[j]; i < m_right[j]; i++) {
			    
			    x10 = m_x1[i]; x11 = m_x1[i+1];
			    x20 = m_x2[j]; x21 = m_x2[j+1];
			    
			    if(m_p[i][j] > 0) 
				p[k] += m_p[i][j] * cell_area(x10,x11,x20,x21,yL,yR,m_op_case);
			}
		}
		
		if(m_op_case == OP_CASE.POW) y = exp(y); // now put it back.
		
		ContinuousRandomVariable Y = new ContinuousRandomVariable(y,p);
		setInfiniteProbability(Y);
		
		return Y;
	    }
	    
	    private void setInfiniteProbability(ContinuousRandomVariable Y) {
		// PosInf
		double p1 = m_pos_inf_x1;
		double p2 = m_pos_inf_x2;
		Y.setPosInf(p1 + p2 - p1*p2);
		// NegInf
		double n1 = m_neg_inf_x1;
		double n2 = m_neg_inf_x2;
		Y.setNegInf(n1 + n2 - n1*n2);
	    }
	    
	    //////////////////////////////////////////////////////////////////////
	    // The Probability Support
	    
	    static public double cell_area(double x10, double x11, double x20, double x21,
					   double yL,  double yR, OP_CASE op_case) {
		double area = 0;
		double a; double b; double c; double d;
		switch(op_case) {
		case ADD:
		    a = yR-x21;
		    b = yR-x20;
		    if(a < x10) a = x10; else if(a > x11) a = x11;
		    if(b < x10) break;   else if(b > x11) b = x11;
		    area = (a-x10)*(x21-x20) + (yR-x20)*(b-a) - (b*b - a*a)/2;
		    c = yL-x21;
		    d = yL-x20;
		    if(c < x10) c = x10; else if(c > x11) break;   // should never happen
		    if(d < x10) break;   else if(d > x11) d = x11;
		    area -= (c-x10)*(x21-x20) + (yL-x20)*(d-c) - (d*d - c*c)/2;
		    break;
		case MULTIPLY:
		    a = yR/x21;
		    b = x11; if(x20 > 0) b = yR/x20;
		    if(a < x10) a = x10; else if(a > x11) a = x11;
		    if(b < x10) break;   else if(b > x11) b = x11;
		    area = (a-x10)*(x21-x20) + yR*Math.log(b/a)-x20*(b-a);
		    if(yL <= 0) break;  // nothing to subtract. bottom cell.
		    c = yL/x21;
		    d = x11; if(x20 > 0) d = yL/x20;
		    if(c < x10) c = x10; else if(c > x11) break;   // should never happen
		    if(d < x10) break;   else if(d > x11) d = x11;
		    area -= (c-x10)*(x21-x20) + yL*Math.log(d/c)-x20*(d-c);
		    break;
		case POW:
		    // convention: yR, yL, a, b are log version of actual.
		    a = yR/x21; if(a < x10) a = x10;
		    if(x20 <= 0) b = x11; else b = yR/x20; if(b > x11) b = x11;
		    if(a >= x11) area = (x21 - x20) * (Math.exp(x11) - Math.exp(x10));
		    else {
			if(a > x10) area = (x21 - x20) * (Math.exp(a) - Math.exp(x10));
			area += yR*Ei.getEiDiff(a,b) - x20 * (Math.exp(b) - Math.exp(a));
		    }
		    //System.out.println("+(a,b): ("+Math.exp(a)+", "+Math.exp(b)+") "+area);
		    if(yL <= 0) break; // was k==0, i.e. 0,0 cell
		    c = yL/x21; if(c < x10) c = x10;
		    if(x20 <= 0) d = x11; else d = yL/x20; if(d > x11) d = x11;
		    if(x10 < d && c < x11) {
			if(c > x10) area -= (x21 - x20) *(Math.exp(c) - Math.exp(x10)); 
			area -= yL*Ei.getEiDiff(c,d) - x20 * (Math.exp(d) - Math.exp(c));
		    }
		    //System.out.println("-(c,d): ("+Math.exp(c)+", "+Math.exp(d)+") "+area);
		    break;
		}
		
		return area;
	    }
	    
	    //////////////////////////////////////////////////////////////////////
	    // Cell Search Limits & Support
	    
	    static public void setSearchLimits(double[][] y, double yL, double yR, 
					       int[] left, int[] right) {
		int N = y.length;
		int M = y[0].length;
		int[] below_point = getBelowCornerPoint(y, yL);
		int[] above_point = getAboveCornerPoint(y, yR);
		for(int i = 0; i < M-1; i++) {left[i] = 0; right[i] = N-1;} // defaults
		setHorizontalLimits(y, below_point, left,  yL, BLOCK_CASE.BELOW);
		setHorizontalLimits(y, above_point, right, yR, BLOCK_CASE.ABOVE);
	    }
	    
	    static public int[] getBelowCornerPoint(double[][] y, double yL) {
		int N = y.length;
		int M = y[0].length;
		if(y[0]  [0]   <= yL) return new int[]{0,   0  };
		if(y[0]  [M-1] <= yL) return new int[]{0,   M-1};
		if(y[N-1][0]   <= yL) return new int[]{N-1, 0  };
		return                       new int[]{N-1, M-1};
	    }
	    
	    static public int[] getAboveCornerPoint(double[][] y, double yR) {
		int N = y.length;
		int M = y[0].length;
		if(y[N-1][M-1] >= yR) return new int[]{N-1, M-1};
		if(y[N-1][0]   >= yR) return new int[]{N-1, 0  };
		if(y[0]  [M-1] >= yR) return new int[]{0,   M-1};
		return                       new int[]{0,   0  };
	    }
	    
	    static public void setHorizontalLimits(double[][] y, int[] start, int[] limits, 
						   double yBound, BLOCK_CASE block_case) {
		//System.out.println("setHLimits(start="+start[0]+","+start[1]+",y="+yBound+")");
		
		int N = y.length;
		int M = y[0].length;
		int L = N * M;                                    // set a search limit
		
		int[] d = getInitialDelta(start);                 // delta
		int[] p = new int[]{start[0], start[1]};          // current point
		int[] q = p;
		while(L-- > 0) {
		    //System.out.println("On step: d="+d+", p="+p+")"); // DEBUG
		    
		    // find a valid direction
		    for(int c = 0; c < 4; c++) {
			q = getSearchPoint(p, d, c);
			int i = q[0]; int j = q[1];
			
			//System.out.println("try i,j="+i+","+j);  // DEBUG
			if(i < 0 || j < 0 || i >= N || j >= M) ;
			else if(block_case == BLOCK_CASE.BELOW && y[i][j] <= yBound) break;
			else if(block_case == BLOCK_CASE.ABOVE && y[i][j] >= yBound) break;
			
			if(c == 3) {
			    //System.out.println("search terminated. no more moves."); // DEBUG
			    return;
			}
		    }
		    
		    //System.out.println("q="+q);  // DEBUG
		    
		    // record the new horizontal limit. It's what we're here for.
		    if(q[0] != start[0] && p[0] == q[0]) {
			int i = Math.min(p[1], q[1]);
			limits[i] = q[0];
			//System.out.println("add search limit: limit["+i+"]="+limits[i]); // DBG
		    }
		    d = new int[]{q[0]-p[0], q[1]-p[1]};
		    p = q;
		    if(p[0] == start[0] && p[1] == start[1]) return;
		}
		
		System.out.println("setHorizontalLimits ran out of steps."); // ERROR
		System.out.println("start: "+start[0]+","+start[1]+
				   " yBound: "+yBound+" case: "+block_case);
	    }
	    
	    static public int[] getInitialDelta(int[] start) {
		// Each corner is different. Must get this initial orientation correct.
		if(start[0] == 0) {
		    if(start[1] == 0) return new int[]{0,1};
		    return new int[]{1,0};
		}
		if(start[1] == 0) return new int[]{-1,0};
		return new int[]{0,-1};
	    }
	    
	    static public int[] getSearchPoint(int[] p, int[] d, int dir) {
		switch(dir) {
		case 0: 
		    if(d[0] == 0) return new int[]{p[0]+d[1],p[1]}; 
		    else          return new int[]{p[0],     p[1]-d[0]};
		case 1:           return new int[]{p[0]+d[0],p[1]+d[1]};
		case 2: 
		    if(d[0] == 0) return new int[]{p[0]-d[1],p[1]};
		    else          return new int[]{p[0],     p[1]+d[0]};
		case 3:           return new int[]{p[0]-d[0],p[1]-d[1]};
		}
		return new int[]{-1,-1}; // intentional error.
	    }
	    
	    //////////////////////////////////////////////////////////////////////
	    // Construction Helpers
	    
	    static private double[] create_partition(double[] x1, double[] x2, 
						     int L, OP_CASE op_case) {
		double[] y  = new double[L];
		double[] y1 = partition(x1,L);
		double[] y2 = partition(x2,L);
		
		switch(op_case) {
		case ADD:      
		    for(int i = 0; i < L; i++) y[i] = y1[i] + y2[i];
		    return y;
		case MULTIPLY: 
		    for(int i = 0; i < L; i++) y[i] = y1[i] * y2[i];
		    return y;
		case POW:
		    double[] logy1 = log(y1);
		    for(int i = 0; i < L; i++) y[i] = y2[i] * logy1[i];    // builds log(y)
		    return y;
		}
		
		return null;
	    }
	    
	    static private double[][] create_y(double[] x1, double[] x2, OP_CASE op_case) {
		int N = x1.length;
		int M = x2.length;
		
		double[][] y = new double[N][M];
		
		switch(op_case) {
		case ADD:
		    for(int i = 0; i < N; i++) for(int j = 0; j < M; j++) y[i][j] = x1[i]+x2[j];
		    return y;
		case MULTIPLY:
		    for(int i = 0; i < N; i++) for(int j = 0; j < M; j++) y[i][j] = x1[i]*x2[j];
		    return y;
		case POW:
		    // ASSUME x1 is log(X1)
		    for(int i = 0; i < N; i++) for(int j = 0; j < M; j++) y[i][j] = x1[i]*x2[j];
		    return y;
		}
		
		return null;
	    }
	    
	    static private double[][] create_p(double[] p1, double[] p2,
					       double[] x1, double[] x2) {
		int N = x1.length;
		int M = x2.length;
		double[][] p = new double[N-1][M-1];
		for(int i = 0; i < N-1; i++) 
		    for(int j = 0; j < M-1; j++)
			p[i][j] = p1[i]/(x1[i+1]-x1[i])*p2[j]/(x2[j+1]-x2[j]);
		return p;
	    }
	    
	    static public double[] partition(double[] x, int L) {
		// ASSUME: Finite points only.
		// Will choose first and last. N is length of x, L is length of result.

		int      N = x.length;
		double[] y = new double[L];
		
		y[0] = x[0];
		for(int j = 1; j < L-1; j++) {
		    double f  = (double)j * (N-1) / (L-1);     // standardized (fractional) index
		    int    i  = (int)f;
		    y[j] = x[i] + (x[i+1] - x[i]) * (f - (double)i);
		}
		y[L-1] = x[N-1];
		
		return y;
	    }
	    
	    static private double sum(double[] x) {
		if(x == null) return 0;
		int N = x.length;
		double v = 0;
		for(int i = 0; i < N; i++) v += x[i];
		return v;
	    }
	    
	    static private double[] log(double[] x) {
		// ASSUME: x > 0
		if(x == null) return null;
		int N = x.length;
		double[] logx = new double[N];
		for(int i = 0; i < N; i++) logx[i] = Math.log(x[i]);
		return logx;
	    }
	    
	    static private double[] exp(double[] x) {
		if(x == null) return null;
		int N = x.length;
		double[] expx = new double[N];
		for(int i = 0; i < N; i++) expx[i] = Math.exp(x[i]);
		return expx;
	    }
	    
	    static private void zero(double[] x) {
		if(x == null) return;
		int N = x.length;
		for(int i = 0; i < N; i++) x[i] = 0;
	    }
	    
	    /////////////////////////////////////////////////////////////
	    // Display
	    
	    static String toString(double[] y) {
		int N = y.length;
		StringBuffer line = new StringBuffer();
		line.append(OPEN_PAREN);
		for(int i = 0; i < N; i++) {
		    line.append(str(y[i]));
		    if(i < N-1) line.append(COMMA_SPACE);
		}
		line.append(CLOSE_PAREN);
		return line.toString();
	    }
	    
	    static String toString(double[][] y) {
		int N = y.length;
		int M = y[0].length;
		StringBuffer line = new StringBuffer();
		line.append(OPEN_PAREN);
		for(int j = M-1; j >= 0; j--) {
		    line.append(OPEN_PAREN);
		    for(int i = 0; i < N; i++) {
			line.append(str(y[i][j]));
			if(i < N-1) line.append(COMMA_SPACE);
		    }
		    line.append(CLOSE_PAREN);
		    if(j > 0) line.append(NEWLINE);
		}
		return line.toString();
	    }
	    
	} // end JointRV
	
	///////////////////////////////////////////////////////////////////
	// Density Interpolation

	// See notes: 3/2/12 p.4
	public static double[] newton_gaussian_step(double[] x, 
						    double a, double b, double c,
						    double p, double q) {
	    final double s2   = Math.sqrt(2);
	    final double is2p = 1/Math.sqrt(2*Math.PI);

	    double u = x[0]; 
	    double s = x[1];

	    try {
		double at = (a-u)/s/s2;
		double bt = (b-u)/s/s2;
		double ct = (c-u)/s/s2;

		double fa = Erf.erf(at);
		double fb = Erf.erf(bt);
		double fc = Erf.erf(ct);

		double g1 = 0.5 * (fb - fa) - p;
		double g2 = 0.5 * (fc - fb) - q;

		System.out.println("g1 = "+g1+" g2 = "+g2); // DEBUG

		double ea = is2p * Math.exp(-at*at) / s; 
		double eb = is2p * Math.exp(-bt*bt) / s; 
		double ec = is2p * Math.exp(-ct*ct) / s; 

		double A = ea - eb;
		double B = s2 * (at*ea - bt*eb);
		double C = eb - ec;
		double D = s2 * (bt*eb - ct*ec);
		double E = A*D - B*C;
		
		x[0] = u - (D*g1 - B*g2) / E;
		x[1] = s - (A*g2 - C*g1) / E;

	    } catch (org.apache.commons.math.MathException e) {
		System.out.println("ContinuousRandomVariable::newtown_step_g error: "+e);
		return null;
	    }

	    return x;
	}

	// Experimental
	public double variance2(double mean) {
	    if(hasPosInf() || hasNegInf()) {
		if(getPosInf() > getNegInf()) return tail().x();
		else                          return head().x();
	    }
	    Node node = first();
	    if(node == null) return Double.NaN;

	    double[] B = get_exponential_factors();
	    double[] A = get_exponential_scalars(B);

	    System.out.print("B = ");
	    for(int i = 0; i < B.length; i++) {
		System.out.print(String.format("%e, ", B[i]));
	    }
	    System.out.println();
	    System.out.print("A = ");
	    for(int i = 0; i < B.length; i++) {
		System.out.print(String.format("%e, ", A[i]));
	    }
	    System.out.println();
	    
	    return 0;
	}

	public double[] get_exponential_factors() {
	    Node node = first();
	    if(node == null) return new double[0];

	    int N = size();
	    if(N <= 2) return new double[]{0};

	    double[] BL = new double[N-2];

	    // ...|p0|p1|p2|p3|p4|0... N = 6

	    // Gather the left-stored exponential factors
	    for(int i = 0; i < N-2; i++) {
		double a = node.x();
		double b = node.next().x();
		double c = node.next().next().x();
		double p = node.p();
		double q = node.next().p();
		BL[i] = newton_exp(a, b, c, p, q);
		node = node.next();
	    }

	    if(N == 3) return new double[]{BL[0],BL[0]};  // no middle cells
	    
	    // mix left- and right- exponential factors for all middle cells
	    double[] B = new double[N-1];
	    node = first();
	    for(int i = 1; i < N-2; i++) {
		double p = node.p();
		double q = node.next().p();
		double r = node.next().next().p();
		double P = p + 2*q + r;
		B[i] = BL[i-1] * (p + q)/P + BL[i] * (q + r)/P;
		node = node.next();
	    }

	    B[0]   = BL[0];
	    B[N-2] = BL[N-3];
	    return B;
	}

	public double[] get_exponential_scalars(double[] B) {
	    Node node = first();
	    double[] A = new double[B.length];

	    for(int i = 0; i < B.length; i++) {

		double a  = node.x();
		double b  = node.next().x();
		double p  = node.p();

		double ea = Math.exp(B[i]*a);
		double eb = Math.exp(B[i]*b);

		A[i] = B[i] * p / (eb - ea);

		node = node.next();
	    }

	    return A;
	}

	// See notes: 3/2/12 p.7
	public static double newton_exp(double a, double b, double c, double p, double q) {
	    double B = p/(b-a) < q/(c-b) ? +1 : -1;  // select polarity
	    
	    int count = 0; int count_limit = 100;

	    while(count++ < count_limit) {
		double B_last = B;
		B = newton_exp_step(B, a, b, c, p, q);
		if(Math.abs(B - B_last) < 1E-15) break;
	    }

	    // DEBUG
	    double ea = Math.exp(B*a);
	    double eb = Math.exp(B*b);
	    double ec = Math.exp(B*c);

	    double A = B * p / (eb - ea);

	    double g = p * (ec - eb) - q * (eb - ea);
      	    System.out.println("steps: "+count+" g = "+g+" p ~ "+(A/B*(eb - ea))+" q ~ "+(A/B*(ec-eb)));
	    // \DEBUG

	    return B;
	}

	public static double newton_exp_step(double B, 
					       double a, double b, double c, 
					       double p, double q) {
	    double ea = Math.exp(B*a);
	    double eb = Math.exp(B*b);
	    double ec = Math.exp(B*c);

	    return B - (p*(ec - eb) - q*(eb - ea))/(p*(c*ec - b*eb) - q*(b*eb - a*ea));
	}
	
	///////////////////////////////////////////////////////////////////
	// Display

	public static String str(double x) {return (new Double(x)).toString();}

	private static void write_line(Writer writer, double x, double y) throws IOException {
	    writer.write(str(x) + "\t" + str(y) + "\n" );
	}

	private static double write_step(Writer writer, Node a, Node b, double last_h) 
	    throws IOException {
	    // returns the height where it left off.
	    // Will accept b == tail()
	    double h = 0;
	    if(a.p() > 0) h = a.p() / (b.x() - a.x());
	    write_line(writer, a.x(), last_h);
	    write_line(writer, a.x(), h);
	    if(!Double.isInfinite(b.x())) write_line(writer, b.x(), h);
	    return h;
	}
	
	public void plot(String filename) throws IOException  {
	    // see: 11/30/10p1
	    Writer writer = new OutputStreamWriter(new FileOutputStream(filename));
	    try {
		Node node = first();
		if(node == null) writer.write("0\t0\n");
		else {
		    double last_h = 0;
		    while(node != tail()) {
			last_h = write_step(writer, node, node.next(), last_h);
			node = node.next();
		    }
		}
	    }
	    finally {
		writer.close();
	    }
	}

	private static double[][] emptyDoublePair() {
	    double[][] d = new double[2][];
	    d[0] = new double[0];
	    d[1] = new double[0];
	    return d;
	}

	public double[][] xh_array() {
	    // i: 0  1  2   j: 0  1  2  3  4  5
	    //   x0 x1 x2 ->  x0 x0 x1 x1 x2 x2
	    //   p0 p1 0       0 h0 h0 h1 h1  0
	    double[][] xp = xp_array();
	    if(xp == null) return emptyDoublePair();
	    double[]   x  = xp[0];
	    double[]   p  = xp[1];
	    int        N  = x.length * 2;
	    if(x.length < 2) return emptyDoublePair();
	    double[]   x2 = new double[N];
	    double[]   h2 = new double[N];
	    double last_h = 0;
	    for(int i = 0; i < x.length-1; i++) {
		int i2   = i*2;
		x2[i2]   = x[i]; 
		x2[i2+1] = x[i];
		h2[i2]   = last_h;
		last_h   = p[i] / (x[i+1] - x[i]);
		h2[i2+1] = last_h;
	    }
	    x2[N-2] = x[x.length-1];
	    x2[N-1] = x[x.length-1];
	    h2[N-2] = last_h;
	    h2[N-1] = 0;
	    return new double[][]{x2,h2};
	}
	
    } // end ContinuousRandomVariable

    //////////////////////////////////////////////////////////////////////////////////

    static public class SimpleRandomVariable {
	Node m_head;
	Node m_tail;

	SimpleRandomVariable() {
	    m_head = new Node(Double.NEGATIVE_INFINITY);
	    m_tail = new Node(Double.POSITIVE_INFINITY);
	    m_head.append(m_tail);
	}

	public SimpleRandomVariable(double[] x, double p[]) {
	    // Allows p.length == x.length - 1 (continuous) == x.length (discrete)
	    m_head = new Node(Double.NEGATIVE_INFINITY);
	    m_tail = new Node(Double.POSITIVE_INFINITY);
	    Node node = m_head;
	    int N = x.length;
	    for(int i = 0; i < N-1; i++) node = node.insert(new Node(x[i],   p[i]));
	    if(p.length < x.length)      node = node.insert(new Node(x[N-1], 0));
	    else                         node = node.insert(new Node(x[N-1], p[N-1]));
	    node.append(tail());
	}

	public SimpleRandomVariable(SimpleRandomVariable A) {
	    m_head = new Node(Double.NEGATIVE_INFINITY);
	    m_tail = new Node(Double.POSITIVE_INFINITY);
	    m_head.append(m_tail);
	    setNegInf(A.getNegInf());
	    setPosInf(A.getPosInf());
	    Node a = A.first();
	    if(a == null) return;
	    Node b = head();
	    while(a != A.tail()) {b = b.insert(new Node(a)); a = a.next();}
	}

	//////////////////////////////////////////////////////////////////////////
	// Access

	public Node head() {return m_head;}
	public Node tail() {return m_tail;}

	public Node first() {if(head().next() != tail()) return head().next(); return null;}
	public Node last()  {if(tail().prev() != head()) return tail().prev(); return null;}

	public boolean hasPosInf()         {return tail().p() > 0;}
	public boolean hasNegInf()         {return head().p() > 0;}
	public double  getPosInf()         {return tail().p();}
	public double  getNegInf()         {return head().p();}
	public void    setPosInf(double v) {tail().p(v);}
	public void    setNegInf(double v) {head().p(v);}
	public void    addPosInf(double v) {tail().p(tail().p() + v);}
	public void    addNegInf(double v) {head().p(head().p() + v);}
	
	public double[][] xp() {
	    int N = size();
	    if(N == 0) return new double[2][0];
	    double[] x = new double[N];
	    double[] p = new double[N];
	    Node node = first();
	    for(int i = 0; i < N; i++) {
		x[i] = node.x(); 
		p[i] = node.p();
		node = node.next();
	    }
	    return new double[][]{x,p};
	}

	////////////////////////////////////////////////////////////////////
	// Info
	
	public int size() {
	    // SLOW. Does not include head & tail
	    int N = 0;
	    Node node = first();
	    if(node == null) return 0;
	    while(node != tail()) {N++; node = node.next();} 
	    return N;
	}

	public double P() {
	    // Good example of how to traverse list
	    double p = 0;
	    Node node = head();
	    while(node != null) {p += node.p(); node = node.next();} 
	    return p;
	}

	public double coreP() {return P() - getNegInf() - getPosInf();}
	
	/////////////////////////////////////////////////////////////////////
	// Operations

	/////////////////////////////////////////////////////////////////////
	// Display

	public String toString() {
	    StringBuffer line = new StringBuffer();
	    line.append(OPEN_PAREN);
	    if(head().p() > 0) {
		line.append(head().toString());
		line.append(COMMA);
	    }
	    Node node = first();
	    if(node == null) {
		if(head().p() > 0 && tail().p() > 0) line.append(COMMA);
	    } else 
		while(node != tail()) {
		    line.append(node.toString());
		    node = node.next();
		}
	    if(tail().p() > 0) line.append(tail().toString());
	    line.append(CLOSE_PAREN);
	    return line.toString();
	}


	
    } // end SimpleRandomVariable

    //////////////////////////////////////////////////////////////////////////////////
    public static class Node {
	/*
	 * Doubly linked list. 
	 * m_x is x-value unless infinite (0)
	 * m_p is associated probability possbibly 
	 *           discrete (at the node) 
	 *        or continusous (between this node and the next
	 *           in which case, last before tail must be zero
      	 */

	private Node m_next = null;
	private Node m_prev = null;
	private double  m_x = 0;
	private double  m_p = 0;
	
	public Node() {}
	public Node(Node node)          {m_x = node.m_x; m_p = node.m_p;}
	public Node(double x)           {m_x = x;}
	public Node(double x, double p) {m_x = x; m_p = p;}

	public Node prev() {return m_prev;}
	public Node next() {return m_next;}

	public double   x()         {return m_x;}
	public double   p()         {return m_p;}
	public void     x(double x) {m_x  = x;}
	public void     p(double p) {m_p  = p;}
	public void add_x(double x) {m_x += x;}
	public void add_p(double p) {m_p += p;}

	public Node append (Node node) {m_next = node; node.m_prev = this; return node;}
	public Node prepend(Node node) {m_prev = node; node.m_next = this; return node;}
	public Node prev   (Node node) {m_prev = node; return node;}
	public Node next   (Node node) {m_next = node; return node;}

	public Node insert(Node node) {
	    Node that = this.next();
	    if(that == null) {append(node); return node;}
	    this.append(node);
	    node.append(that);
	    return node;
	}

	public Node remove() {
	    Node p = prev();
	    Node n = next();
	    if(p != null) {p. append(n); return p;}
	    if(n != null) {n.prepend(p); return n;}
	    return null;
	}

	//////////////////////////////////////////////////////////
	// Display
	
	public String toString() {
	    StringBuffer line = new StringBuffer();

	    line.append(str(x()));

	    if(p() > 0) {
		line.append(OPEN_PAREN);
		line.append(str(p()));
		line.append(CLOSE_PAREN);
	    } else {
		if(prev() != null && next() != null)
		    line.append(COMMA);
	    }

	    return line.toString();
	}

    } // end Node

    //////////////////////////////////////////////////////////////////////////////////
    // Statics

    public static Boolean is_zero(double x) {
	return Math.abs(x) <= DTOL;
    }

    public static Boolean strict_equal(double x, double y) {
	return Math.abs(y-x) <= DTOL;
    }

    private static int compare(Node A, Node B) {
	if(strict_equal(A.x(), B.x()))  return 0;
	if(A.x() < B.x()) return -1;
	return 1;
    }

    public static int[] isort(double[] arr) {
        int N = arr.length;
        if(N == 0) return new int[0];
        int     [] index = new int[N];
        SortPair[] array = new SortPair[N];
        for(int i = 0; i < N; i++) array[i] = new SortPair(arr[i], i);
        Arrays.sort(array);
        for(int i = 0; i < N; i++) index[i] = array[i].getIndex();
        return index;
    }

    public static double[] sort(double[] arr, int[] index) {
        int N = arr.length;
        double[] x = new double[N];
        for(int i = 0; i < N; i++) x[i] = arr[index[i]];
	return x;
    }

    static class SortPair implements Comparable<SortPair> {
	private double m_value;
	private int    m_index;

	public SortPair() {}
	public SortPair(double value, int index) {m_value = value; m_index = index;}
	public void set(double value, int index) {m_value = value; m_index = index;}

	public double getValue() {return m_value;}
	public int    getIndex() {return m_index;}

	public int compareTo( SortPair other ) { 
	    // may need some work on what to do with null's
	    // System.out.println("in compareTo()");
	    if ( this == null && other == null ) return 0;
	    if(m_value == other.m_value) return 0;
	    if(m_value <  other.m_value) return -1;
	    return 1;
	}
    }

    public static String str(double x) {return (new Double(x)).toString();}

    public static String toString(double[] y) {
	int N = y.length;
	StringBuffer line = new StringBuffer();
	line.append(OPEN_PAREN);
	for(int i = 0; i < N; i++) {
	    line.append(str(y[i]));
	    if(i < N-1) line.append(COMMA_SPACE);
	}
	line.append(CLOSE_PAREN);
	return line.toString();
    }

    public static boolean half_fraction(double d) {
	if(d < 0) return strict_equal(2*d % 2, -1);  // d == -N/2 for N = {0,1,2,...}
	else      return strict_equal(2*d % 2, +1);  // d ==  N/2 for N = {0,1,2,...}
    }

    public static double pow(double a, double b) {
	// NOT INTENDED FOR SPECIAL CASES
	// ASSUME: a != 0, b != {0,+-1} 
	if(a > 0) return Math.pow(a,b);
	if(half_fraction(b)) return 0;
	return Math.cos(Math.PI * b) * Math.pow(-a,b);
    }

    public static int left_integer(double a) {
	if(a > 0) return (int)a;
	if(a < 0) return -((int)-a)-1;
	return 0;
    }

    public static int right_integer(double a) {
	if(a > 0) return (int)a+1;
	if(a < 0) return -((int)-a);
	return 0;
    }

    public static boolean isEven(int n) {
	if(n < 0) return (-n % 2 == 0);
	else      return   n % 2 == 0;
    }

    public static int[] allocate(int N, double[] weight) {
	int i_max_weight = -1;
	double total_weight = 0;
	double max_weight = -1;
	for(int i = 0; i < weight.length; i++) {
	    if(weight[i] > max_weight) {max_weight = weight[i]; i_max_weight = i;}
	    total_weight += weight[i];
	}
	int N_allocated = 0;
	int[] count = new int[weight.length];
	for(int i = 0; i < weight.length; i++) {
	    count[i] = (int)(N * weight[i] / (double)total_weight);
	    N_allocated += count[i];
	}
	count[i_max_weight] = N - N_allocated;  // fix off-by-one error, if present.
	return count;
    }

} // end NumericRandomVariable


