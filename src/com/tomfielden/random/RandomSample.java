package com.tomfielden.random;

import java.util.*;
import java.io.*;
import org.apache.commons.math.MathException;
import org.apache.commons.math.distribution.*;
import org.apache.commons.math.random.*;
import org.apache.commons.math.stat.descriptive.SummaryStatistics;

public class RandomSample {
    static private final String COMMA         = ",";
    static private final String COMMA_SPACE   = ", ";
    static private final String NEWLINE       = "\n";
    static private final String OPEN_PAREN    = "(";
    static private final String CLOSE_PAREN   = ")";

    static int  DEFAULT_BINS = 500;
    static long m_seed = 17399225432l;

    RandomData randomData = new RandomDataImpl(); 

    public static void SetSeed(long seed) {m_seed = seed;}
    public static void SetSeed() {

    }

    double[] m_sample;

    // Constructors

    public RandomSample()               { m_sample = null;}
    public RandomSample(Double[] a)     { create_new(a); }
    public RandomSample(double[] a)     { create_new(a); }
    public RandomSample(RandomSample A) { create_new(A.sample()); }

    private void create_new(double[] a) {
	m_sample = new double[a.length];
	for(int i = 0; i < a.length; i++) m_sample[i] = a[i];
    }
    private void create_new(Double[] a) {
	m_sample = new double[a.length];
	for(int i = 0; i < a.length; i++) m_sample[i] = a[i];
    }

    // Static Constructors

    public static RandomSample Cauchy(double median, double scale, int N) {
	try {
	    CauchyDistributionImpl distribution = new CauchyDistributionImpl(median, scale);
	    distribution.reseedRandomGenerator(m_seed);
	    return new RandomSample(distribution.sample(N));
	} catch (org.apache.commons.math.MathException me) {
	    System.out.println(me);
	    return null;
	}
    }

    public static RandomSample ChiSquared(double df, int N) {
	try {
	    ChiSquaredDistributionImpl distribution = new ChiSquaredDistributionImpl(df);
	    //distribution.reseedRandomGenerator(m_seed);
	    return new RandomSample(distribution.sample(N));
	} catch (org.apache.commons.math.MathException me) {
	    System.out.println(me);
	    return null;
	}
    }

    public static RandomSample Exponential(double lambda, int N) {
	try {
	    ExponentialDistributionImpl distribution = new ExponentialDistributionImpl(lambda);
	    //distribution.reseedRandomGenerator(m_seed);
	    return new RandomSample(distribution.sample(N));
	} catch (org.apache.commons.math.MathException me) {
	    System.out.println(me);
	    return null;
	}
    }

    public static RandomSample Normal(double mu, double sigma, int N) {
	try {
	    NormalDistributionImpl distribution = new NormalDistributionImpl(mu, sigma);
	    //distribution.reseedRandomGenerator(m_seed);
	    return new RandomSample(distribution.sample(N));
	} catch (org.apache.commons.math.MathException me) {
	    System.out.println(me);
	    return null;
	}
    }

    ////////////////////////////////////////////////////////////////////////////////
    // Access

    public int      size()   {return m_sample.length;}
    public double[] sample() {return m_sample;}
    public double   sum()    {
	double s = 0; 
	for(int i = 0; i < size(); i++) s += m_sample[i]; 
	return s;
    }

    public double       get(int i) {return sample()[getIndex(i, size())];}
    public RandomSample get(int M, int start, int step) {
	double[] a = sample();
	double[] c = new double[M];
	for(int i = 0; i < M; i++) c[i] = get(start + i*step);
	return new RandomSample(c);
    }
    public RandomSample get(Boolean[] mask) {
	double[] a = sample();
	int N = size();
	int M = 0;
	for(int i = 0; i < mask.length; i++) if(mask[i]) M++;
	double[] c = new double[M];
	int j = 0;
	for(int i = 0; i < N; i++) if(mask[i]) c[j++] = a[i];
	return new RandomSample(c);
    }

    public List binStats(int N) {
	EmpiricalDistributionImpl empirical = new EmpiricalDistributionImpl(N);
	empirical.load(m_sample);
	return empirical.getBinStats();
    }

    public NumericRandomVariable getDiscreteRandomVariable(int bins) {
	double N = (double) size();

	double[] x = new double[bins];
	double[] p = new double[bins];

	int i = 0;
	List<SummaryStatistics> stats = binStats(bins);

	for(SummaryStatistics s : stats) {
	    long n = s.getN();
	    if(n == 0) continue;
	    p[i] = n / N;
	    x[i] = s.getMean();
	    i++;
	}

	double[] x1 = new double[i];
	double[] p1 = new double[i];
	for(int j = 0; j < i; j++) {x1[j] = x[j]; p1[j] = p[j];}

	return NumericRandomVariable.Discrete(x1, p1);
    }

    public NumericRandomVariable getContinuousRandomVariable(int bins) {
	double N = (double) size();

	int i = 0;
	List<SummaryStatistics> stats = binStats(bins);

	//   0    1    2    3    4   
	// L q RL q RL q RL q RL q R  i=5
	// 0   01   12   23   34   4

	double[] xL = new double[bins];
	double[] xR = new double[bins];
	double[] q  = new double[bins];

	for(SummaryStatistics s : stats) {
	    long n = s.getN();
	    if(n == 0) continue;
	    xL[i] = s.getMin();
	    xR[i] = s.getMax();
	    q[i]  = n / N;
	    i++;
	}

	//   0   1   2   3   4
	// | p | p | p | p | p |   i=5
	// 0   1   2   3   4   5
	double[] x = new double[i+1];
	double[] p = new double[i+1];

	x[0] = xL[0];   p[0] = q[0];
	x[i] = xR[i-1]; p[i] = 0;
	for(int j = 1; j < i; j++) {
	    p[j] = q[j];
	    x[j] = (xR[j-1] + xL[j])/2;
	}

	return NumericRandomVariable.Continuous(x, p);
    }

    public NumericRandomVariable getContinuousRandomVariable() {
	return getContinuousRandomVariable(DEFAULT_BINS);
    }

    ////////////////////////////////////////////////////////////////////////////////
    // Operations

    public RandomSample sort() {
	RandomSample B = new RandomSample(sample());
	double[] b = B.sample();
	Arrays.sort(b);
	return B;
    }

    public RandomSample negative() {
	double[] a = sample();
	double[] c = new double[a.length];
	for(int i = 0; i < a.length; i++) c[i] = -a[i];
	return new RandomSample(c);
    }

    public RandomSample reciprocal() {
	double[] a = sample();
	double[] c = new double[a.length];
	for(int i = 0; i < a.length; i++) c[i] = 1/a[i];
	return new RandomSample(c);
    }

    public RandomSample add(double b) {
	double[] a = sample();
	double[] c = new double[a.length];
	for(int i = 0; i < a.length; i++) c[i] = a[i] + b;
	return new RandomSample(c);
    }

    public RandomSample add(RandomSample B) {
	double[] a = sample();
	double[] b = B.sample();
	double[] c = new double[a.length];
	for(int i = 0; i < a.length; i++) c[i] = a[i] + b[i]; // assume parallel
	return new RandomSample(c);
    }

    public RandomSample multiply(double b) {
	double[] a = sample();
	double[] c = new double[a.length];
	for(int i = 0; i < a.length; i++) c[i] = a[i] * b;
	return new RandomSample(c);
    }

    public RandomSample multiply(RandomSample B) {
	double[] a = sample();
	double[] b = B.sample();
	double[] c = new double[a.length];
	for(int i = 0; i < a.length; i++) c[i] = a[i] * b[i]; // assume parallel
	return new RandomSample(c);
    }

    // Comparisons

    public Boolean[] LT(double x) {
	double[]  a = sample();
	Boolean[] c = new Boolean[a.length];
	for(int i = 0; i < a.length; i++) c[i] = a[i] < x;
	return c;
    }

    public Boolean[] LT(RandomSample B) {
	double[]  a = sample();
	double[]  b = B.sample();
	Boolean[] c = new Boolean[a.length];
	for(int i = 0; i < a.length; i++) c[i] = a[i] < b[i];
	return c;
    }

    public Boolean[] LTE(double x) {
	double[]  a = sample();
	Boolean[] c = new Boolean[a.length];
	for(int i = 0; i < a.length; i++) c[i] = a[i] <= x;
	return c;
    }

    public Boolean[] LTE(RandomSample B) {
	double[]  a = sample();
	double[]  b = B.sample();
	Boolean[] c = new Boolean[a.length];
	for(int i = 0; i < a.length; i++) c[i] = a[i] <= b[i];
	return c;
    }

    public Boolean[] EQ(double x) {
	double[]  a = sample();
	Boolean[] c = new Boolean[a.length];
	for(int i = 0; i < a.length; i++) c[i] = a[i] == x;
	return c;
    }

    public Boolean[] EQ(RandomSample B) {
	double[]  a = sample();
	double[]  b = B.sample();
	Boolean[] c = new Boolean[a.length];
	for(int i = 0; i < a.length; i++) c[i] = a[i] == b[i];
	return c;
    }

    public Boolean[] GT(double x) {
	double[]  a = sample();
	Boolean[] c = new Boolean[a.length];
	for(int i = 0; i < a.length; i++) c[i] = a[i] > x;
	return c;
    }

    public Boolean[] GT(RandomSample B) {
	double[]  a = sample();
	double[]  b = B.sample();
	Boolean[] c = new Boolean[a.length];
	for(int i = 0; i < a.length; i++) c[i] = a[i] > b[i];
	return c;
    }

    public Boolean[] GTE(double x) {
	double[]  a = sample();
	Boolean[] c = new Boolean[a.length];
	for(int i = 0; i < a.length; i++) c[i] = a[i] >= x;
	return c;
    }

    public Boolean[] GTE(RandomSample B) {
	double[]  a = sample();
	double[]  b = B.sample();
	Boolean[] c = new Boolean[a.length];
	for(int i = 0; i < a.length; i++) c[i] = a[i] >= b[i];
	return c;
    }


    //////////////////////////////////////////////////////////////////////////////////
    // Display
    
    public String toString() {
	return toString(sample());
    }

    // helpers //

    static String str(double x) {return (new Double(x)).toString();}

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
    
    static int getIndex(int i, int N) {
	// deal with modulo issues
	if(i > 0) return i % N;
	return ((1-i/N)*N + i) % N;
    }

    ///////////////////////////////////////////////////////////////////////////////////
} // end RandomSample