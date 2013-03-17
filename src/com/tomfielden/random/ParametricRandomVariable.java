package com.tomfielden.random;

import java.util.*;
import java.io.*;

/*
ParametricRandomVariable is intended to be a leaf in the LazyRandomVariable parse tree.
Once sampled, the results are cached.
*/

public class ParametricRandomVariable {

    static private int ID_CODE = 1;                 // DOUBLE will always be ID=0
    static private int newID() {return ID_CODE++;}

    public enum VTYPE {DOUBLE,NUMERIC, UNIFORM, BETA, CAUCHY, CHISQ, EXP, 
	    FDIST, LOGNORMAL, NORMAL, POISSON, STUDENT, TRI};

    static private int DEFAULT_SAMPLE_SIZE = 1000;
    private static double ZERO_TOL         = 1E-15;
    private static boolean SHOW_ID_CODE    = true;

    public static void ShowIdCodes() { SHOW_ID_CODE = true;  }
    public static void HideIdCodes() { SHOW_ID_CODE = false; }

    ////////////////////////////////////////////////////////////////////////////////////////////////
    // Members

    private VTYPE                 m_type;
    private int                   m_id;
    private double                m_param1;      // DOUBLE stored here.
    private double                m_param2;
    private double                m_param3;
    private NumericRandomVariable m_nrv = null;  // also caches samples. Implies no resampling.

    ////////////////////////////////////////////////////////////////////////////////////////////////
    // Constructors

    public ParametricRandomVariable(double param1) {
	m_id     = 0;
	m_type   = VTYPE.DOUBLE;
	m_param1 = param1;
    }

    public ParametricRandomVariable(ParametricRandomVariable P) {
	m_id     = newID();
	m_type   = P.m_type;
	m_param1 = P.m_param1;
	m_param2 = P.m_param2;
	m_param3 = P.m_param3;
	m_nrv    = P.m_nrv;
    }

    private ParametricRandomVariable(NumericRandomVariable X) {
	m_id     = newID();
	m_type   = VTYPE.NUMERIC;
	m_nrv    = X;
    } 

    private ParametricRandomVariable(VTYPE type, double param1) {
	m_id     = newID();
	m_type   = type;
	m_param1 = param1;
    } 

    private ParametricRandomVariable(VTYPE type, double param1, double param2) {
	m_id     = newID();
	m_type   = type;
	m_param1 = param1;
	m_param2 = param2;
    } 

    private ParametricRandomVariable(VTYPE type, double param1, double param2, double param3) {
	m_id     = newID();
	m_type   = type;
	m_param1 = param1;
	m_param2 = param2;
	m_param3 = param3;
    } 

    ////////////////////////////////////////////////////////////////////////////////////////////////
    // Static Constructors

    static public ParametricRandomVariable Double(double d) {
	return new ParametricRandomVariable(d);
    }

    static public ParametricRandomVariable Numeric(NumericRandomVariable X) {
	return new ParametricRandomVariable(X);
    }

    static public ParametricRandomVariable Discrete(double[] x, double[] p) {
	return new ParametricRandomVariable(NumericRandomVariable.Discrete(x, p));
    }

    public static ParametricRandomVariable Continuous(double[] x, double[] p) {
	return new ParametricRandomVariable(NumericRandomVariable.Continuous(x, p));
    }

    public static ParametricRandomVariable Uniform(double a, double b) {
	return new ParametricRandomVariable(VTYPE.UNIFORM, a, b);
    }

    public static ParametricRandomVariable Beta(double alpha, double beta) {
	return new ParametricRandomVariable(VTYPE.BETA, alpha, beta);
    }

    public static ParametricRandomVariable Cauchy(double median, double scale) {
	return new ParametricRandomVariable(VTYPE.CAUCHY, median, scale);
    }

    public static ParametricRandomVariable ChiSquared(double df) {
	return new ParametricRandomVariable(VTYPE.CHISQ, df);
    }

    public static ParametricRandomVariable Exponential(double lambda) {
	return new ParametricRandomVariable(VTYPE.EXP, lambda);
    }

    public static ParametricRandomVariable F(double df1, double df2) {
        return new ParametricRandomVariable(VTYPE.FDIST, df1, df2);
    }

    public static ParametricRandomVariable Normal(double mu, double sigma) {
	return new ParametricRandomVariable(VTYPE.NORMAL, mu, sigma);
    }

    public static ParametricRandomVariable LogNormal(double mu, double sigma) {
	return new ParametricRandomVariable(VTYPE.LOGNORMAL, mu, sigma);
    }

    public static ParametricRandomVariable Poisson(double lambda) {
	return new ParametricRandomVariable(VTYPE.POISSON, lambda);
    }
    
    public static ParametricRandomVariable Student(double median, double scale, double nu) {
	return new ParametricRandomVariable(VTYPE.STUDENT, median, scale, nu);
    }

    public static ParametricRandomVariable Triangular(double L, double M, double R) {
	return new ParametricRandomVariable(VTYPE.TRI, L, M, R);
    }

    ////////////////////////////////////////////////////////////////////////////////////////////////
    // Info

    public int    id()     {return m_id;}
    public VTYPE  type()   {return m_type;}
    public double value()  {return m_param1;} 

    public boolean isDouble() { return type() == VTYPE.DOUBLE; }

    public double P() {
	if(type() == VTYPE.NUMERIC) return m_nrv.P();
	else                        return 1.0;
    }

    public double E() {
	switch(type()) {
	case DOUBLE:    return m_param1;
	case UNIFORM:   return (m_param1 + m_param2)/2;
	case NUMERIC:   return get().E();
	case NORMAL:    return m_param1;
	case LOGNORMAL: return Math.exp(m_param1 + m_param2 * m_param2 / 2);
	case STUDENT:   return m_param1;
	}
	System.out.println("ParametricRandomVariable.E() type "+type()+" not supported.");
	return Double.NEGATIVE_INFINITY;
    }

    public double median() {
	switch(type()) {
	case DOUBLE:    return m_param1;
	case UNIFORM:   return E();
	case NUMERIC:   return get().median();
	case NORMAL:    return m_param1;
	case LOGNORMAL: return Math.exp(m_param1);
	}
	System.out.println("ParametricRandomVariable.median() type "+type()+" not supported.");
	return Double.NEGATIVE_INFINITY;
    }

    public double stdev() {
	switch(type()) {
	case DOUBLE:    return 0;
	case UNIFORM:   return (m_param2 - m_param1) / 3.4641016151377544;  // sqrt(12)
	case NUMERIC:   return get().stdev();
	case NORMAL:    return m_param2;
	case LOGNORMAL: 
	    double s2 = m_param2 * m_param2;
	    return Math.sqrt((Math.exp(s2) - 1) * Math.exp(2*m_param1 + s2));
	}
	System.out.println("ParametricRandomVariable.median() type "+type()+" not supported.");
	return Double.NEGATIVE_INFINITY;
    }

    //public double[][] getData() {return get().getData();}  // ensure numeric rv present.

    ////////////////////////////////////////////////////////////////////////////////////////////////
    // Methods

    public NumericRandomVariable get() {return sample(DEFAULT_SAMPLE_SIZE);}
   
    public NumericRandomVariable sample() {return sample(DEFAULT_SAMPLE_SIZE);}  // DEPRECATED

    public NumericRandomVariable sample(int n) {
	if(m_nrv != null) return m_nrv; // NUMERIC will always hit this. CACHE
	switch(type()) {
	case UNIFORM:   m_nrv = NumericRandomVariable.Uniform    (m_param1, m_param2);              break;
	case BETA:      m_nrv = NumericRandomVariable.Beta       (m_param1, m_param2, n);           break;
	case CAUCHY:    m_nrv = NumericRandomVariable.Cauchy     (m_param1, m_param2, n);           break;
	case CHISQ:     m_nrv = NumericRandomVariable.ChiSquared (m_param1,           n);           break;
	case DOUBLE:    m_nrv = NumericRandomVariable.Discrete   (m_param1);                        break;
	case EXP:       m_nrv = NumericRandomVariable.Exponential(m_param1,           n);           break;
	case FDIST:     m_nrv = NumericRandomVariable.F          (m_param1, m_param2, n);           break;
	case LOGNORMAL: m_nrv = NumericRandomVariable.LogNormal  (m_param1, m_param2, n);           break;
	case NORMAL:    m_nrv = NumericRandomVariable.Normal     (m_param1, m_param2, n);           break;
	case POISSON:   m_nrv = NumericRandomVariable.Poisson    (m_param1);                        break;
	case STUDENT:   m_nrv = NumericRandomVariable.Student    (m_param1, m_param2, m_param3, n); break;
	case TRI:       m_nrv = NumericRandomVariable.Triangular (m_param1, m_param2, m_param3, n); break;
	}
	return m_nrv;   // cached sample
    }

    //////////////////////////////////////////////////////////////////////////////////
    // Aritmetic Operations

    private void NotImplemented(String s) {
	System.out.println("ParametricRandomVariable Function Not Implemented: "+s);
    }

    // TODO: Compute symbolic values and fall back to numeric if not available.
    public double LT (double x) {return get().LT(x);}
    public double LTE(double x) {return get().LTE(x);}
    public double EQ (double x) {return get().EQ(x);}
    public double GT (double x) {return get().GT(x);}
    public double GTE(double x) {return get().GTE(x);}

    public double LT (ParametricRandomVariable B) {NotImplemented("LT");  return Double.NEGATIVE_INFINITY;}
    public double LTE(ParametricRandomVariable B) {NotImplemented("LTE"); return Double.NEGATIVE_INFINITY;}
    public double EQ (ParametricRandomVariable B) {NotImplemented("EQ");  return Double.NEGATIVE_INFINITY;}
    public double GT (ParametricRandomVariable B) {NotImplemented("GT");  return Double.NEGATIVE_INFINITY;}
    public double GTE(ParametricRandomVariable B) {NotImplemented("GTE"); return Double.NEGATIVE_INFINITY;}

    public ParametricRandomVariable neg() {
	switch(type()) {
	case DOUBLE:
	    return Double(-m_param1);
	case NORMAL:
	    return Normal(-m_param1, m_param2);
	}
	return Numeric(get().neg());
    }

    public ParametricRandomVariable sqrt() {
	switch(type()) {
	case DOUBLE:
	    return Double(Math.sqrt(m_param1));
	case LOGNORMAL:
	    return LogNormal(m_param1 / 2, 0.7071067811865476 * m_param2); // sqrt(2)
	}
	return Numeric(get().sqrt());
    }

    public ParametricRandomVariable log() {
	switch(type()) {
	case DOUBLE:
	    return Double(Math.log(m_param1));
	case LOGNORMAL:
	    return Normal(m_param1,m_param2);
	}
	return Numeric(get().log());
    }

    public ParametricRandomVariable exp() {
	switch(type()) {
	case DOUBLE:
	    return Double(Math.exp(m_param1));
	case NORMAL:
	    return LogNormal(m_param1,m_param2);
	}
	return Numeric(get().exp());
    }

    public ParametricRandomVariable reciprocal() {
	switch(type()) {
	case DOUBLE:
	    return Double(1/m_param1);
	case LOGNORMAL:
	    return LogNormal(-m_param1,m_param2);
	}
	return Numeric(get().reciprocal());
    }

    public ParametricRandomVariable add(ParametricRandomVariable B) {
	switch(type()) {
	case DOUBLE:
	    if(B.type() == VTYPE.DOUBLE) return Double(m_param1 + B.m_param1);
	case NUMERIC:
	    switch(B.type()) {
	    case DOUBLE:
		return Numeric(get().add(B.m_param1));
	    }
	    break;
	case LOGNORMAL:
	    switch(B.type()) {
	    case DOUBLE:
		return Numeric(get().add(B.m_param1));
	    case LOGNORMAL:
		if(id() == B.id()) return multiply(Double(2.0));
	    }
	    break;
	case NORMAL:
	    switch(B.type()) {
	    case DOUBLE:
		return Normal(m_param1 + B.m_param1, m_param2);
	    case NORMAL:
		if(id() == B.id()) return multiply(Double(2.0));
	    }
	    break;
	}
	return Numeric(get().add(B.get()));
    }

    public ParametricRandomVariable sub(ParametricRandomVariable B) {
	switch(type()) {
	case DOUBLE:
	    if(B.type() == VTYPE.DOUBLE) return Double(m_param1 - B.m_param1);
	case NUMERIC:
	    switch(B.type()) {
	    case DOUBLE:
	    case NUMERIC:
		if(id() == B.id()) return Double(0.0);
	    }
	    break;
	case LOGNORMAL:
	    switch(B.type()) {
	    case DOUBLE:
		return Numeric(get().sub(B.m_param1));
	    case NUMERIC:
		return Numeric(B.get().sub(get()));
	    case LOGNORMAL:
		if(id() == B.id()) return Double(0.0);
	    }
	    break;
	case NORMAL:
	    switch(B.type()) {
	    case DOUBLE:
		return Normal(m_param1 - B.m_param1, m_param2);
	    }
	    break;
	}
	return Numeric(get().sub(B.get()));
    }

    public ParametricRandomVariable multiply(ParametricRandomVariable B) {
	switch(type()) {
	case DOUBLE:
	    if(B.type() == VTYPE.DOUBLE) return Double(m_param1 * B.m_param1);
	    return B.multiply(this);
	case EXP:
	    if(B.type() == VTYPE.DOUBLE && B.m_param1 > 0) return Exponential(m_param1 * B.m_param1);
	    break;
	case NUMERIC:
	    switch(B.type()) {
	    case NUMERIC:
		if(id() == B.id()) return Numeric(get().pow(2));
	    }
	    break;
	case LOGNORMAL:
	    switch(B.type()) {
	    case DOUBLE:
		return LogNormal(m_param1 + Math.log(B.m_param1), m_param2);
	    case LOGNORMAL:
		if(id() == B.id()) return pow(Double(2));
		return LogNormal(m_param1 + B.m_param1, Math.sqrt(m_param2 * m_param2 + B.m_param2 * B.m_param2));
	    }
	    break;
	case NORMAL:
	    switch(B.type()) {
	    case DOUBLE:
		return Normal(m_param1 * B.m_param1, m_param2 * B.m_param1);
	    case NORMAL:
		if(id() == B.id()) return pow(Double(2));
	    }
	    break;
	}
	return Numeric(get().multiply(B.get()));
    }

    public ParametricRandomVariable pow(ParametricRandomVariable B) {
	switch(type()) {
	case DOUBLE:
	    if(B.type() == VTYPE.DOUBLE) return Double(Math.pow(m_param1,B.m_param1));
	    return B.pow(this);
	case LOGNORMAL:
	    switch(B.type()) {
	    case DOUBLE:
		return LogNormal(B.m_param1 * m_param1, Math.sqrt(B.m_param1)*m_param2);
	    }
	    break;
	case NORMAL:
	    switch(B.type()) {
	    case DOUBLE:
		return Numeric(get().pow(B.m_param1));
	    }
	    break;
	}
	return Numeric(get().pow(B.get()));
    }

    //////////////////////////////////////////////////////////////////////////////////
    // Output

    public String toString() {
	if(SHOW_ID_CODE && id() > 0) return String.format("%s_%d", toStringUndecorated(), id());
	else                         return toStringUndecorated();
    }

    private String sparam(int i) { 
	if(i == 1) return toStringDouble(m_param1); 
	if(i == 2) return toStringDouble(m_param2); 
	if(i == 3) return toStringDouble(m_param3); 
	return "Unknown param: "+i;
    }

    private String toStringUndecorated() {
	switch(type()) {
	case NUMERIC:   return String.format("Numeric");
	case UNIFORM:   return String.format("U[%s,%s]",          sparam(1), sparam(2));
	case BETA:      return String.format("Beta[%s,%s]",       sparam(1), sparam(2));
	case CAUCHY:    return String.format("Cauchy[%s,%s]",     sparam(1), sparam(2));
	case CHISQ:     return String.format("ChiSq[%s]",         sparam(1));
	case EXP:       return String.format("Exp[%s]",           sparam(1));
	case DOUBLE:    return toStringDouble(m_param1);
	case FDIST:     return String.format("F[%s,%s]",          sparam(1), sparam(2));
	case NORMAL:    return String.format("N[%s,%s]",          sparam(1), sparam(2));
	case LOGNORMAL: return String.format("LogN[%s,%s]",       sparam(1), sparam(2));
	case POISSON:   return String.format("Poisson[%s]",       sparam(1));
	case STUDENT:   return String.format("Student[%s,%s,%s]", sparam(1), sparam(2), sparam(3));
	case TRI:       return String.format("Tri[%s,%s,%s]",     sparam(1), sparam(2), sparam(3));
	}
	return "Unknown ParametricRandomVariable";
    }

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

} // end interface
