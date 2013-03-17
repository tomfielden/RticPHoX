package com.tomfielden.random;

import java.util.*;
import java.io.*;

import org.apache.commons.math.special.Gamma;

/* Usage: phox >>>
from com.tomfielden.random import *
alpha = 1.5
beta  = 0.51
L = LevyStableDistribution(alpha, beta)
L.cdf(1)
*/

public class LevyStableDistribution {

    static final double MODE_TOL = 1E-5;
    static final double EPSILON  = 1E-12;

    public TanhSinhParameters m_params_pos;
    public TanhSinhParameters m_params_neg;

    // Assume symmetry so only store non-negative indices.
    static public int halflevelsize(int level)  {if(level == 0) return 2; return pow2(level-1);}
    ArrayList<ArrayList<Double>> m_phi  = new ArrayList<ArrayList<Double>>();  // anti-symmetric
    ArrayList<ArrayList<Double>> m_dphi = new ArrayList<ArrayList<Double>>();  // symmetric

    public LevyStableDistribution(double alpha, double beta) {
	m_params_pos = new TanhSinhParameters(this, alpha,  beta);
	m_params_neg = new TanhSinhParameters(this, alpha, -beta);
    }

    // This is the distribution-based construction for Numeric Random Variables.
    public static double[][] Density(double alpha,double beta,double location,double scale,int N) {

	double[]   z  = partition(alpha, beta, N);
	double[]   x  = new double[N];
	double[]   p  = new double[N];
	double[][] xp = new double[2][];
	xp[0] = x;
	xp[1] = p;
	
	LevyStableDistribution LSD = new LevyStableDistribution(alpha, beta);
	
	double c         = LSD.cdf(z[0]);          // first cumulative probability
	double d         = LSD.cdf(z[N-1]);        // last cumulative probability  (!!)
	double MP        = (c + 1-d);              // "missing" probability
	
	System.out.println("LevyStable missing probability: "+MP);  // DEBUG
	
	double m = LSD.mode();

	x[0]   = (z[0] + m) * scale + location;
	p[N-1] = 0;
	
	// The probability is applied to last node, position to new node.
	for(int i = 1; i < N; i++) {
	    //System.out.println("LevyStable F(z = "+z+")"); // DEBUG
	    double c1 = LSD.cdf(z[i] + m);
	    p[i-1]    = (c1 - c)/(1-MP);    // new probability scaled to include endpoint fragments
	    x[i]      = (z[i] + m) * scale + location; // shift by mode
	    c         = c1;
	}

	//System.out.println("Finished computing LevyStable Density");
	
	return xp;
    }

    // Uses density not distribution function.
    public static double[][] pdf(double alpha,double beta,double location,double scale,int N) {
	double[]   z  = partition(alpha, beta, N);
	double[]   x  = new double[N];
	double[]   p  = new double[N];
	double[][] xp = new double[2][];
	xp[0] = x;
	xp[1] = p;
	
	LevyStableDistribution LSD = new LevyStableDistribution(alpha, beta);

	for(int i = 0; i < N; i++) {
	    x[i] = z[i] * scale + location;
	    p[i] = LSD.pdf(z[i]);
	}

	return xp;
    }

    public double cdf(double x) {
	double alpha = m_params_pos.alpha();
	double beta  = m_params_pos.beta();
	double zeta  = m_params_pos.zeta();
	if(alpha == 1) {
	    if(beta  > 0) return m_params_pos.cdf( x);
	    if(beta == 0) return 0.5 + Math.atan(x)/pi;
	    if(beta  < 0) return 1 - m_params_neg.cdf(x);
	}
	if(x < zeta) return 1 - m_params_neg.cdf(-x);
	else         return     m_params_pos.cdf( x);
    }

    public double pdf(double x) {
	// NB: not yet updated to handle alpha == 1, beta != 0
	double zeta = m_params_pos.zeta();
	if(x < zeta) return m_params_neg.pdf(-x);
	else         return m_params_pos.pdf( x);
    }

    // Partition Calculation

    public static double[] partition(double alpha, double beta, int N) {
	double[] z = new double[N];

	double   knee  = 0.85;         // percentage points within quasi-linear sampling range
	
	// See notes: 2/28/12 p.5. Solved for parameters to acheive empirical exponents
	// alpha = .25 => L = 1E50, alpha = 1   => L = 1E18, alpha = 2 => L = 1E5
	final double a = -14.7191011236;
	final double b =  42.0224719101;
	final double c =   0.516853932584;

	double W = (a*alpha + b) / (alpha + c);  // use W as the width and exponent.
	double L = Math.pow(10, W);

	for(int i = 0; i < N; i++) z[i]  = sample(i, N, knee, W, L);

	return z;
    }

    private static double sample(int i, int N, double knee, double W, double L) {
	double z = 2.0 * i / (N-1) - 1;  // z(i=0) = -1, z(i=N-1) = 1
	
	if(z < 0) return -sample_value(-z, knee, W, L);
	else      return  sample_value( z, knee, W, L);
    }

    public static double sample_value(double z, double knee, double W, double L) {
	// See notes: 1/23/2012 p.3 && 2/28/2012 p.3 (copy)
	// knee is percent breakpoint for quasi-linear region
	// w is the target value for end of quasi-linear region. Will actually be 2w. (why?)
	
	double ik   = 1/knee;
	double z_k  = ik*z;
	double z_k2 = z_k * z_k;
	double LWik = L/W - ik*ik;

	if(LWik < 1) return 2 * W * z_k2;  // restrict to quasi-linear region only
	
	double y = Math.log(LWik) / Math.log(ik);

	return W * (z_k2 + Math.pow(z_k,y));
    }

    // Mode-finding density maximization /////

    public double mode() {if(m_params_pos.beta() == 0) return 0; return bisection();}

    public double bisection() {
	double[][] triple = new double[2][3];
	triple[0][0] = -2; triple[1][0] = pdf(-2);
	triple[0][1] =  0; triple[1][1] = pdf( 0);
	triple[0][2] =  2; triple[1][2] = pdf( 2);
	final int COUNT_MAX = 500;
	int count = 0;
	while(triple[0][2] - triple[0][0] > MODE_TOL && count++ < COUNT_MAX) bisection_step(triple);
	return triple[0][1];
    }

    // See Notes 2/29/12 p.2
    public double[][] bisection_step(double[][] triple) {
	double x1 = triple[0][0]; double f1 = triple[1][0];
	double x2 = triple[0][1]; double f2 = triple[1][1];
	double x3 = triple[0][2]; double f3 = triple[1][2];

	double x12 = (x1+x2)/2; double f12 = pdf(x12);
	double x23 = (x2+x3)/2;

	if(f12 < f1) {
	    System.out.println("bisection_step WARNING: "+f12+" = f12 < f1,f2"); 
	} else if(f12 > f2) {
	    triple[0][1] = x12; triple[1][1] = f12;
	    triple[0][2] = x2;  triple[1][2] = f2;
	} else {
	    triple[0][0] = x12;  triple[1][0] = f12;
	    double f23 = pdf(x23);
	    if(f23 < f3) {
		System.out.println("bisection_step WARNING: "+f23+" = f23 < f2,f3"); 
	    } else if(f23 > f2) {
		triple[0][0] = x2;  triple[1][0] = f2;
		triple[0][1] = x23; triple[1][1] = f23;
	    } else {
		triple[0][2] = x23;  triple[1][2] = f23;

		// We've moved in both endpoints, but not the middle. Try again.
		if(triple[0][2] - triple[0][0] > MODE_TOL) 
		    bisection_step(triple);
	    }
	}
	//System.out.println("Bisection Step: x* = "+triple[0][1]); // DEBUG
	return triple;
    }

    // The Integration Quadrature //////////////////////////////////////////////////////////////

    static public class TanhSinhParameters {
	static final int    LEVEL_MIN   = 4;
	static final int    LEVEL_LIMIT = 12;  // should be 12
	static final double LEVEL_TOL   = 1E-15;

	static public int fulllevelsize(int level)  {if(level == 0) return 3; return pow2(level);}

	LevyStableDistribution  m_distribution;

	double m_alpha;
	double m_beta;
	double m_zeta;
	double m_theta0;
	double m_c0;
	double m_c1;
	double m_c3;
	double m_slope;
	double m_intercept;
	
	public TanhSinhParameters(LevyStableDistribution distribution, double alpha, double beta) {
	    m_distribution = distribution;
	    m_alpha        = alpha;
	    m_beta         = beta;
	    m_zeta         = zeta     (alpha, beta);
	    m_theta0       = theta0   (alpha, beta);
	    m_c0           = c0       (alpha, beta);
	    m_c1           = c1       (alpha, beta);
	    m_c3           = c3       (alpha);
	    m_slope        = slope    (alpha, beta);
	    m_intercept    = intercept(alpha, beta);
	}

	// The Distribution Function

	public double cdf(double x) {
	    // Assume x >= zeta
	    if(alpha() == 1 && beta() == 0) return 0.5 + Math.atan(x) / pi;
	    if(x == zeta()) return c0();

	    int    level     = 0;   // must always start at zero or miss lower-level values
	    double I         = 0.0;
	    double ICUM      = 0.0;
	    double ICUM_last = 0.0;

	    while(level <= LEVEL_LIMIT) {
		double value = 0.0;
		int    n     = fulllevelsize(level);
		for(int i = 0; i < n; i++) {
		    double theta = theta(level, i);
		    double dphi  = dphi (level, i);
		    double expg  = expg(theta, x);
		    //System.out.println("level= "+level+" i= "+i+" theta= "+theta+" x= "+x+" expg= "+expg+" dphi= "+dphi); // DEBUG
		    value += expg * dphi;
		}
		I += value * slope();

		ICUM_last = ICUM;
		double h  = tsN / pow2(level);
		ICUM      = h * I;

		if(level >= LEVEL_MIN && Math.abs(ICUM - ICUM_last) < LEVEL_TOL) break;

		//System.out.println("level = "+level+" I = "+ICUM);

		level ++;
	    }
	    
	    return c1() + c3() * ICUM;
	}

	// The Density Function

	public double pdf(double x) {
	    // Assume x >= zeta
	    if(alpha() == 1 && beta() == 0) return 1 / (pi * (1 + x*x));
	    if(x == zeta()) return c4();

	    int    level     = 0;   // must always start at zero or miss lower-level values
	    double I         = 0.0;
	    double ICUM      = 0.0;
	    double ICUM_last = 0.0;

	    while(level <= LEVEL_LIMIT) {
		double value = 0.0;
		int    n     = fulllevelsize(level);
		for(int i = 0; i < n; i++) {
		    double theta = theta(level, i);
		    double dphi  = dphi (level, i);
		    double expg  = expg(theta, x);
		    double g     = g(theta, x);
		    double dval  = dphi * expg * g;
		    if(g < 1E44) value += dval;  // i.e. really big
		    //System.out.println("level= "+level+" i = "+i+" theta= "+theta+" dval= "+dval);
		}
		I += value * slope();

		ICUM_last = ICUM;
		double h  = tsN / pow2(level);
		ICUM      = h * I;

		if(level >= LEVEL_MIN && Math.abs(ICUM - ICUM_last) < LEVEL_TOL) break;

		//System.out.println("level = "+level+" I = "+ICUM);

		level ++;
	    }
	    
	    return c2(x) * ICUM;
	}

	// Access functions

	public double alpha()     {return m_alpha;}
	public double beta()      {return m_beta;}
	public double zeta()      {return m_zeta;}
	public double theta0()    {return m_theta0;}
	public double c0()        {return m_c0;}
	public double c1()        {return m_c1;}
	public double c3()        {return m_c3;}
	public double slope()     {return m_slope;}
	public double intercept() {return m_intercept;}

	// Setting Functions
	
	public static double zeta(double alpha, double beta) {
	    if(alpha == 1) return 0;
	    return - beta * Math.tan(pi_2 * alpha);
	}
	
	public static double theta0(double alpha, double beta) {
	    if(alpha == 1) return pi_2;
	    return 1.0/alpha * Math.atan(beta * Math.tan(pi_2 * alpha));
	}
	
	public static double c0(double alpha, double beta) {
	    return 0.5 - theta0(alpha, beta) / pi;
	}
	
	public static double c1(double alpha, double beta) {
	    if(alpha == 1) return 0.0;
	    if(alpha >  1) return 1.0;
	    return c0(alpha, beta);
	}

	public double c2(double x) {
	    double alpha = alpha();
	    double beta  = beta();
	    double zeta  = zeta();
	    if(alpha == 1) return 0.5 / Math.abs(beta);
	    return alpha / (pi * Math.abs(alpha - 1) * (x - zeta));
	}

	public static double c3(double alpha) {
	    if(alpha == 1) return 1/pi;
	    return sign(1 - alpha)/pi;
	}

	public double c4() {
	    double ialpha = 1/alpha();
	    double zeta2   = zeta() * zeta();
	    double theta0 = theta0();
	    return gamma(1 + ialpha) * Math.cos(theta0) / (pi * Math.pow(1 + zeta2, ialpha/2));
	}

	public static double gamma(double x) {
	    return Math.exp(Gamma.logGamma(x));
	}

	// NB: We don't want to integrate all the way to pi/2 since it's really hard numerically.
	public static double slope(double alpha, double beta) {
	    return (pi_2 - EPSILON + theta0(alpha, beta)) / 2;
	}
	public static double intercept(double alpha, double beta) {
	    return (pi_2 - EPSILON - theta0(alpha, beta)) / 2;
	}

	////////////////////////////////////////////////////////////////////////////////////////////
	// Tanh-Sinh computation level caching scheme. See: Notes 2/26/12 p1 && 2/27/12 p1
	// Level 0: t-values = (-E 0 E),    index/2 = (-1 0 1),       index = (0 1 2)
	// Level 1: t-values = (-2 2),      index/2 = (-0 0),         index = (0 1)
	// Level 2: t-values = (-3 -1 1 3), index/2 = (-1, -0, 0, 1), index = (0 1 2 3)
	// NB: The sign of index/2 computed separately 

	// Accepts fulllevelsize-based indicies 0,1,..,2^level-1 and accounts for symmetry
	public double theta(int level, int index) {
	    return slope() * phi(level, index) + intercept();
	}
	public double  phi(int level, int index) {
	    int i = full_to_half_level(level, index);
	    return  phi(level).get(i) * sign_level_index(level, index);
	}
	public double dphi(int level, int index) {
	    int i = full_to_half_level(level, index);
	    return dphi(level).get(i);
	}
	public static int full_to_half_level(int level, int index) {
	    if(level == 0) return Math.abs(index-1);
	    int i =  index - pow2(level-1);
	    return i < 0 ? -(i + 1) : i;
	}
	public static int sign_level_index(int level, int index) {
	    if(level == 0) return index == 0 ? -1 : 1;
	    return sign(index - pow2(level-1));
	}
	
	// Intended to be internal computation of minimally cached results.
	public ArrayList<Double> phi(int level) {
	    while(phi().size() <= level) {
		ArrayList<Double> p = new ArrayList<Double>();
		if(phi().size() == 0) {
		    p.add( phi(0.0) );
		    p.add( phi(tsN) );
		} else {
		    int    n = halflevelsize(level);
		    double h = tsN / Math.pow(2, level);
		    for(int i = 0; i < n; i++) p. add( phi( h*(2*i + 1) ));
		}
		phi().add(p);
	    }
	    return phi().get(level);
	}
	public ArrayList<Double> dphi(int level) {
	    while(dphi().size() <= level) {
		ArrayList<Double> dp = new ArrayList<Double>();
		if(dphi().size() == 0) {
		    dp.add( dphi(0.0) );
		    dp.add( dphi(tsN) );
		} else {
		    int    n = halflevelsize(level);
		    double h = tsN / pow2(level);
		    for(int i = 0; i < n; i++) dp.add(dphi( h*(2*i + 1) ));
		}
		dphi().add(dp);
	    }
	    return dphi().get(level);
	}

	public ArrayList<ArrayList<Double>>  phi() {return m_distribution. m_phi;}
	public ArrayList<ArrayList<Double>> dphi() {return m_distribution.m_dphi;}

	static public double  phi(double t) {
	    return Math.tanh(pi_2 * Math.sinh(t));
	}
	static public double dphi(double t) {
	    return pi_2 * Math.cosh(t) / Math.pow( Math.cosh( pi_2 * Math.sinh(t) ), 2 );
	}
	
	// Integration Functions

	public double V(double theta) {

	    double alpha    = alpha();
	    double beta     = beta();
	    double theta0   = theta0();
	    double costheta = Math.cos(theta);

	    if(alpha == 1) {
		double D = (1 +  beta * theta / pi) / costheta;
		double E = Math.exp((pi_2 / beta + theta) * Math.tan(theta));
      		//System.out.println("V("+theta+")= "+D+"*"+E); // DEBUG
		return D * E;
	    } else {
		if(theta < -theta0) return Double.POSITIVE_INFINITY;
		double A1 = alpha - 1.0;

		// old method
		/*
		double D  = Math.pow(Math.cos(alpha * theta0), 1.0 / A1);
		double E  = Math.pow(costheta /
				     Math.sin(alpha * (theta0 + theta)),
				     alpha / A1);
		double F  = Math.cos(alpha * theta0 + A1 * theta) / costheta;
		*/
		// new method

		double G = Math.pow(Math.cos(alpha * theta0) * costheta, 1/A1);
		double H = Math.cos(alpha*(theta0+theta)-theta);
		double J = Math.pow(Math.sin(alpha*(theta0+theta)), alpha/A1);

		//if(Math.abs(H) < ZERO_TOL) H = Math.abs(H); Don't need
		
		//System.out.println("alpha*(theta0+theta)-theta= "+(alpha*(theta0+theta)-theta));
      		//System.out.println("V("+theta+")= "+G+"*"+H+"/"+J); // DEBUG
		//System.out.println("V = "+(D*E*F)+" or "+(G*H/J));    // DEBUG

		//return D * E * F;
		return G*H/J;
	    }
	}

	public double g(double theta, double x) {
	    double alpha = alpha();
	    double beta  = beta();
	    double zeta  = zeta();
	    double V     = V(theta);

	    if(alpha == 1) {
		double E = Math.exp(- pi_2 * x / beta);
		//System.out.println("g("+theta+","+x+")= "+E+"*"+V+" = "+(E*V));
		return E * V;
	    } else {
		double E = Math.pow(x - zeta, alpha / (alpha - 1.0));
		return E * V;
	    }
	}

	public double expg(double theta, double x) {
	    return Math.exp( -g(theta, x) );
	}

	public static final double pi   = Math.PI;
	public static final double pi_2 = Math.PI/2;

	// 0 => +1
	public static int sign(double x) {if(x < 0) return -1; return +1;}  

    } // end TanhSinhParameters ////////////////////////////////////////////////////////////////////

    public static final double pi   = Math.PI;
    public static final double pi_2 = Math.PI/2;
    public static final double tsN  = 4.0;        // radius of the tanh-sinh quadrature

    public static int pow2(int i)    {
	int e = 1; 
	for(int j = 0; j < i; j++) e *= 2; 
	return e;
    }

} // end class