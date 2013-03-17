package com.tomfielden.util;

import org.python.core.*;
import org.python.util.PythonInterpreter;
import org.apache.commons.math.analysis.interpolation.*; 
import org.apache.commons.math.distribution.*;

import java.util.*;
import java.io.*;

import com.tomfielden.random.LazyRandomVariable;
import com.tomfielden.random.NumericRandomVariable;
import com.tomfielden.random.RandomSample;

/*
 * This class is inteded to facilitate access from Jython to interesting functions. Numpy-like.
 * NB: We assume for performance reasons that Variant does NOT have a copy constructor.
 *
 */

public class Functions {

    /////////////////////////////////////////////////////////////////////////////////////
    // Constants

    public static final double pi = 3.1415926535897932384626433832795028841971693993751058209;

    /////////////////////////////////////////////////////////////////////////////////////
    // Creation

    public static Variant createIndex(PyObject x) {
	return Variant.createIndex(new Variant(x));
    }

    public static Variant createMask(PyObject x) {
	return Variant.createMask(new Variant(x));
    }

    public static Variant createPartition(PyObject x) {
	return Variant.createPartition(new Variant(x));
    }

    public static Variant createVector(PyObject x) {
	return Variant.createVector(new Variant(x));
    }

    public static Variant createSparseDiagonal(PyObject x) {
	return Variant.createSparseDiagonal(new Variant(x));
    }

    public static Variant createSparseDiagonals(PyObject indicies, PyObject diags) {
	return Variant.createSparseDiagonal(new Variant(indicies), new Variant(diags));
    }

    public static Variant createSparseMatrix(PyObject r, PyObject c, PyObject b) {
	return Variant.createSparseMatrix(new Variant(r), new Variant(c), new Variant(b));
    }

    public static Variant createSparseMatrixHorizontal(PyObject x) {
	return Variant.createSparseMatrixHorizontal(new Variant(x));
    }

    public static Variant createSparseMatrixVertical(PyObject  x) {
	return Variant.createSparseMatrixVertical(new Variant(x));
    }

    public static Variant createSparseMatrixDiagonal(PyObject x) {
	return Variant.createSparseMatrixDiagonal(new Variant(x));
    }

    /////////////////////////////////////////////////////////////////////////////////////
    // Math

    public static Variant exp(PyObject a) {
	return Variant.exp(new Variant(a));
    }

    public static Variant log(PyObject a) {
	return Variant.log(new Variant(a));
    }

    public static Variant sqrt(PyObject a) {
	return Variant.sqrt(new Variant(a));
    }

    public static Variant min(PyObject a, PyObject b) {
	return Variant.min(new Variant(a), new Variant(b));
    }

    public static Variant max(PyObject a, PyObject b) {
	return Variant.max(new Variant(a), new Variant(b));
    }

    public static Variant diff (PyObject A) { return Variant.diff (new Variant(A)); }
    public static Variant diff0(PyObject A) { return Variant.diff0(new Variant(A)); }
    public static Variant diff1(PyObject A) { return Variant.diff1(new Variant(A)); }

    public static Variant sum(PyObject A) {
	return (new Variant(A)).sum();
    }
    
    public static Variant cumsum(PyObject A) {
	return (new Variant(A)).cumsum();
    }

    public static Variant DCF(PyObject R, PyObject CF) {
	return Variant.DCF(new Variant(R), new Variant(CF));
    }

    // Discounts start in first year of cash flows unless CF = number of years.
    public static Variant NPV(PyObject R, PyObject CF) {
	return Variant.NPV(new Variant(R), new Variant(CF));
    }

    // Cash Flow Value under Depreciation
    public static Variant CFV(PyObject CF, PyObject Depreciation) { 
	return Variant.CFV(new Variant(CF), new Variant(Depreciation));
    }

    // Cash Flow corrected for Depreciation
    public static Variant CFD(PyObject CF, PyObject D) {
	return Variant.CFD(new Variant(CF), new Variant(D));
    }

    // DEPRECATED. Use CFD.
    public static Variant asset_service(PyObject A, PyObject D, int N) {
	return Variant.asset_service(new Variant(A), new Variant(D), N);
    }

    // Only works under certain circumstances.
    public static Variant asset_service_steady_state(PyObject A, PyObject D) {
	return Variant.asset_service_steady_state(new Variant(A), new Variant(D));
    }

    public static Variant PMT(PyObject R, PyObject lifetime) {
	return Variant.PMT(new Variant(R), new Variant(lifetime));
    }

    public static Variant PLDCF(PyObject R, PyObject lifetime, 
				PyObject periods, PyObject accums) {
	return Variant.PLDCF(new Variant(R), new Variant(lifetime), 
			     new Variant(periods), new Variant(accums));
    }

    // Cumulative Cash Flows under growth rate R.
    public static Variant CCF(PyObject R, PyObject CF) {
	return Variant.CCF(new Variant(R), new Variant(CF));
    }

    public static Variant convolve(PyObject F, PyObject G) {
	return Variant.convolve(new Variant(F), new Variant(G));
    }

    public static Variant choose(PyObject N, PyObject K) {
	return Variant.choose(new Variant(N), new Variant(K));
    }

    public static Variant where(PyObject mask, PyObject LHS, PyObject RHS) {
	return Variant.where(new Variant(mask), new Variant(LHS), new Variant(RHS));
    }

    public static Variant div0(PyObject N, PyObject D) {
	return Variant.div0(new Variant(N), new Variant(D));
    }

    public static Variant cases(PyObject field, PyObject map, PyObject dflt) {
	return Variant.cases(new Variant(field), new Variant(map), new Variant(dflt));
    }

    public static Variant ones(PyObject size) {
	return Variant.ones(new Variant(size));
    }

    public static Variant zeros(PyObject size) {
	return Variant.zeros(new Variant(size));
    }

    public static Variant arange(PyObject stop) {
	return Variant.arange(new Variant(stop));
    }

    public static Variant arange2(PyObject start, PyObject stop) {
	return Variant.arange(new Variant(start), new Variant(stop));
    }

    public static Variant outer(PyObject A, PyObject B) {
	return Variant.outer(new Variant(A), new Variant(B));
    }

    public static Variant hadamard(PyObject x, PyObject y) {
	return Variant.hadamard(new Variant(x), new Variant(y));
    }

    public static Variant linspace(PyObject start, PyObject stop, PyObject count) {
	return Variant.linspace(new Variant(start), new Variant(stop), new Variant(count));
    }

    public static Variant Loess(Variant x, Variant y)
	throws org.apache.commons.math.MathException {
	return RegressionFactory.Loess(x, y);
    }

    public static double sin(double x) {
	return Math.sin(x);
    }

    public static double cos(double x) {
	return Math.cos(x);
    }

    public static double random() {
	return Math.random();
    }

    /////////////////////////////////////////////////////////////////////////////////////
    // Shaping

    public static Variant cat(PyObject A) {
	return (new Variant(A)).cat();
    }

    public static Variant transpose(PyObject A) {
	return (new Variant(A)).transpose();
    }

    public static Variant rshift(PyObject A) {
	return Variant.rshift(new Variant(A));
    }

    /////////////////////////////////////////////////////////////////////////////////////
    // File/Database handling

    public static Variant open(String filename) {
	Datafile table = new Datafile(filename);
	return table.getColumnMap();
    }

    public static String read(String filename) throws IOException {
	StringBuilder text    = new StringBuilder();
	String        NL      = System.getProperty("line.separator");
	Scanner       scanner = new Scanner(new FileInputStream(filename), "UTF-8");
	try {
	    while (scanner.hasNextLine()) text.append(scanner.nextLine() + NL);
	} finally {
	    scanner.close();
	}
	return text.toString();
    }

    /////////////////////////////////////////////////////////////////////////////////////
    // Jython interpretation

    public static Variant call(String filename) throws IOException {
	return interpret2(read(filename), new Registry());
    }

    public static Variant call2(String filename, PyObject X) throws IOException {
	return interpret2(read(filename), new Registry(new Variant(X)));
    }

    public static Variant interpret(String code) {
	return interpret2(code, new Registry());
    }

    public static Variant interpret2(String code, Registry R) {
	try {
	    PythonInterpreter interp = new PythonInterpreter();

	    // inject given variables
	    Set<String> keys = R.keys().toStringSet();
	    for(String k: keys) interp.set(k, R.get(k));

	    interp.exec("from com.tomfielden.util import *");
	    interp.exec("from com.tomfielden.util.Functions import *");

	    Variant orig_locals = new Variant(interp.getLocals());
	    interp.exec(code);
	    Variant cur_locals  = new Variant(interp.getLocals());

	    return Variant.sub(cur_locals, orig_locals);   // CAPITAL variables always excluded.

	} catch(Exception e) {
	    System.out.println(e);
	    return new Variant();
	}
    }

    /////////////////////////////////////////////////////////////////////////////////////
    // LazyRandomVariables

    public static Variant ChiSquared(PyObject df) {
	return new Variant(LazyRandomVariable.ChiSquared((new Variant(df)).toDouble()));
    }

    public static Variant Cauchy(PyObject median, PyObject scale) {
	return new Variant(LazyRandomVariable.Cauchy((new Variant(median)).toDouble(), 
						     (new Variant(scale)).toDouble()));
    }

    public static Variant Exponential(PyObject lambda) {
	return new Variant(LazyRandomVariable.Exponential((new Variant(lambda)).toDouble()));
    }

    public static Variant Normal(PyObject mu, PyObject sigma) {
	return new Variant(LazyRandomVariable.Normal((new Variant(mu)).toDouble(), 
						     (new Variant(sigma)).toDouble()));
    }

    public static Variant NormalEquivalent(Variant X) {
	return new Variant(LazyRandomVariable.Normal(X.mean().toDouble(), 
						     X.stdev().toDouble()));
    }

    public static Variant LogNormal(PyObject mu, PyObject sigma) {
	return new Variant(LazyRandomVariable.LogNormal((new Variant(mu)).toDouble(), 
							(new Variant(sigma)).toDouble()));
    }

    public static Variant LogNormalEquivalent(Variant X) {
	double median = X.median().toDouble();
	double stdev  = X.stdev().toDouble();
	// Given an observed median and stdev we find the paremeters for the fitted LogNormal
	double mu    = Math.log(median);
	double a     = 2*Math.log(stdev) - 2*mu;
	double sigma = Math.sqrt( Math.log(0.5*(Math.sqrt(1+4*Math.exp(a))+1)) );
	System.out.println("mu: "+mu+", sigma: "+sigma); // DEBUG
	return new Variant(LazyRandomVariable.LogNormal(mu, sigma));
    }

    public static Variant Student(PyObject median, PyObject scale, PyObject nu) {
	return new Variant(LazyRandomVariable.Student((new Variant(median)).toDouble(), 
						      (new Variant(scale)).toDouble(),
						      (new Variant(nu)).toDouble() ));
    }

    public static Variant Uniform(PyObject a, PyObject b) {
	return new Variant(LazyRandomVariable.Uniform((new Variant(a)).toDouble(), 
						      (new Variant(b)).toDouble()));
    }

    /////////////////////////////////////////////////////////////////////////////////////
    // NumericRandomVariables

    public static Variant DiscreteNumeric(PyObject x, PyObject p) {
	return new Variant(NumericRandomVariable.Discrete((new Variant(x)).getVector(),
							  (new Variant(p)).getVector()));
    }

    public static Variant ContinuousNumeric(PyObject x, PyObject p) {
	return new Variant(NumericRandomVariable.Continuous((new Variant(x)).getVector(),
							    (new Variant(p)).getVector()));
    }

    public static Variant CauchyNumeric(PyObject median, PyObject scale, PyObject n) {
	return new Variant(NumericRandomVariable.Cauchy((new Variant(median)).toDouble(), 
						        (new Variant(scale)). toDouble(), 
						        (new Variant(n)).     toInteger()));
    }

    public static Variant ChiSquaredNumeric(PyObject df, PyObject n) {
	return new Variant(NumericRandomVariable.ChiSquared((new Variant(df)).toDouble(), 
							    (new Variant(n)). toInteger()));
    }

    public static Variant ExponentialNumeric(PyObject lambda, PyObject n) {
	return new Variant(NumericRandomVariable.Exponential((new Variant(lambda)).toDouble(), 
							     (new Variant(n)).     toInteger()));
    }

    public static Variant NormalNumeric(PyObject mu, PyObject sigma, PyObject n) {
	Variant Mu    = new Variant(mu);
	Variant Sigma = new Variant(sigma);
	int     N     = new Variant(n).toInteger();
	
	if(Mu.isNumber()) {
	    if(Sigma.isNumber())
		return new Variant(NumericRandomVariable.Normal(Mu.toDouble(),
								Sigma.toDouble(),
								N));
	    if(Sigma.isRV())
		return new Variant(NumericRandomVariable.Normal(Mu.toDouble(),
								Sigma.getNumericRandomVariable(),
								N));
	}
	if(Mu.isRV()) {
	    if(Sigma.isNumber())
		return new Variant(NumericRandomVariable.Normal(Mu.getNumericRandomVariable(),
								Sigma.toDouble(),
								N));
	}
	System.out.println("NormalNumeric("+Mu.type()+", "+Sigma.type()+") unsupported");
	return null;
    }

    public static Variant LogNormalNumeric(PyObject mu, PyObject sigma, PyObject n) {
	return new Variant(NumericRandomVariable.LogNormal((new Variant(mu)).   toDouble(), 
							   (new Variant(sigma)).toDouble(), 
							   (new Variant(n)).    toInteger()));
    }

    public static Variant LevyStableNumeric(PyObject alpha, PyObject beta, 
					    PyObject location, PyObject scale, PyObject n) {
	return new Variant(NumericRandomVariable.LevyStable((new Variant(alpha)).   toDouble(),
							    (new Variant(beta)).    toDouble(),
							    (new Variant(location)).toDouble(),
							    (new Variant(scale)).   toDouble(),
							    (new Variant(n)).       toInteger()));
    }

    public static Variant PoissonNumeric(PyObject lambda) {
	return new Variant(NumericRandomVariable.Poisson((new Variant(lambda)).toDouble()));
    }

    public static Variant StudentNumeric(PyObject median, PyObject scale, PyObject nu, PyObject n) {
	return new Variant(NumericRandomVariable.Student((new Variant(median)).toDouble(),
							 (new Variant(scale)). toDouble(),
							 (new Variant(nu)).    toDouble(), 
							 (new Variant(n)).     toInteger()));
    }

    public static Variant UniformNumeric(PyObject a, PyObject b) {
	return new Variant(NumericRandomVariable.Uniform((new Variant(a)).toDouble(), 
						         (new Variant(b)).toDouble()));
    }

    /////////////////////////////////////////////////////////////////////////////////////
    // RandomSample

    public static Variant CreateSample(PyObject array) {
	return new Variant(new RandomSample((new Variant(array)).getVector()));
    }

    public static Variant ExponentialSample(PyObject lambda, PyObject n) {
	return new Variant(RandomSample.Exponential((new Variant(lambda)).toDouble(), 
						    (new Variant(n)). toInteger()));
    }

    public static Variant ChiSquaredSample(PyObject df, PyObject n) {
	return new Variant(RandomSample.ChiSquared((new Variant(df)).toDouble(), 
						   (new Variant(n)). toInteger()));
    }

    public static Variant NormalSample(PyObject mu, PyObject sigma, PyObject n) {
	return new Variant(RandomSample.Normal((new Variant(mu)).   toDouble(), 
					       (new Variant(sigma)).toDouble(), 
					       (new Variant(n)).    toInteger()));
    }

    /////////////////////////////////////////////////////////////////////////////////////
    // Special Functions

    public static double BlackScholes(double So, double K, double sigma) {
	try {
	    NormalDistributionImpl Normal = new NormalDistributionImpl(); // Apache Math-Commons
	    double d1 = Math.log(So / K)/sigma + sigma/2;
	    double d2 = Math.log(So / K)/sigma - sigma/2;
	    return So * Normal.cumulativeProbability(d1) - K * Normal.cumulativeProbability(d2);
	} catch (org.apache.commons.math.MathException e) {
	    System.out.println("BlackScholes exceptiont: "+e);
	    return Double.NEGATIVE_INFINITY;
	}
    }


    ////////////////////////////////////////////////////////////////////////////////////////////
} // End Class