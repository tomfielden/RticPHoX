package com.tomfielden.regression;

import java.util.*;
import java.io.*;

import org.apache.commons.math.stat.regression.*;
import org.apache.commons.math.distribution.*;

public class OLSRegression  {

    private int                         m_N;
    private int                         m_K;
    private double[]                    m_y        = null;
    private OLSMultipleLinearRegression m_multi    = null;  // only if multi
    private double[][]                  m_multi_x  = null;  // only if multi
    private SimpleRegression            m_simple   = null;  // only if simple
    private double[]                    m_simple_x = null;  // only if simple

    public OLSRegression(double[] y, double[] x) {
	// Simple regression
	m_N          = y.length;
	m_K          = 2;
	m_y          = y;
	m_simple_x   = x;
	m_simple     = new SimpleRegression();
	double[][] d = new double[y.length][2];
	for(int i = 0; i < y.length; i++) {
	    d[i][0] = x[i];
	    d[i][1] = y[i];
	}
	m_simple.addData(d);
    }

    public OLSRegression(double[] y, double[][] x) {
	// Assume x ~ R[N,K], i.e. (row, col)
	// PrintArray(x); // DEBUG
	m_N       = x.length;
	m_K       = x[0].length + 1;  // for the ones column we'll add
	m_y       = y;
	m_multi_x = x;
	m_multi   = new OLSMultipleLinearRegression();
	m_multi.newSampleData(y, x);
    }

    // Methods //

    public int N() {return m_N;}
    public int K() {return m_K;}

    public double[] Beta() {
	if(m_multi != null) {
	    return m_multi.estimateRegressionParameters();
	} else {
	    double[] b = new double[2]; 
	    b[0] = m_simple.getIntercept();
	    b[1] = m_simple.getSlope();
	    return b;
	}
    }
    public double varError() {
	// s^2 = RSS/(N-k)  where RSS = u'u and u = y - x Beta (residuals)
	if(m_multi != null)  return m_multi.estimateErrorVariance();
	else                 return m_simple.getMeanSquareError();
    }
    public double seError() {
	return java.lang.Math.sqrt(varError());
    }
    public double[][] varBeta() {
	// var(Beta) = s^2 * (x'x)^-1
	if(m_multi != null) {
	    double[][] xxi = m_multi.estimateRegressionParametersVariance(); // (X'X)^{-1}
	    double     s2  = m_multi.estimateErrorVariance();
	    return in_mult(xxi, s2);
	} else {
	    int    n   = m_simple_x.length;
	    double s2  = m_simple.getMeanSquareError();
	    double sx  = getSum(m_simple_x);
	    double sx2 = getSumSquares(m_simple_x);
	    double det = n*sx2 - sx*sx;
	    double[][] xxi = new double[2][2];
	    xxi[0][0] = s2 * sx2 / det;
	    xxi[0][1] = s2 * (-sx) / det;
	    xxi[1][0] = xxi[0][1];
	    xxi[1][1] = s2 * n / det;
	    return xxi;
	}
    }
    public double[] seBeta() {
	if(m_multi != null) {
	    return m_multi.estimateRegressionParametersStandardErrors();
	} else {
	    double[] se = new double[2]; 
	    se[0] = m_simple.getInterceptStdErr();
	    se[1] = m_simple.getSlopeStdErr();
	    return se;
	}
    }
    public double seY() {
	return java.lang.Math.sqrt(TSS()/(N() -1));
    }
    public double TSS() { // Total Sum of Squares
	if(m_multi != null) return m_multi.calculateTotalSumOfSquares();
	else                return m_simple.getTotalSumSquares();
    }
    public double ESS() { // Estimated Sum of Squares	
	if(m_multi != null)  return TSS() - RSS();
	else                 return m_simple.getRegressionSumSquares();
    }
    public double RSS() { // Residual Sum of Squares
	if(m_multi != null) return m_multi.calculateResidualSumOfSquares();
	else                return m_simple.getSumSquaredErrors();
    }
 
    public double R2() {
	if(m_multi != null) return m_multi.calculateRSquared();
	else                return m_simple.getRSquare();
    }

    // Information //////////////////////////////////////////////////////////

    // Don't use. See: Variant.correlation()
    public double[][] correlationData() {
	if(m_multi == null) {
	    double[][] cor = new double[1][1];
	    cor[0][0] = 1;
	    return cor;
	} else
	    return correlation(m_multi_x);
    }

    // Significance Tests ///////////////////////////////////////////////////

    public double[] tScores() {	return div(Beta(), seBeta()); }

    public double[] pValues() throws org.apache.commons.math.MathException {
	// two-tailed probability => 2*.
	int           df = N() - K();
	TDistribution T  = new TDistributionImpl(df);
	double[]      t  = tScores();
	double[]      p  = new double[t.length];
	for(int i = 0; i < t.length; i++) 
	    p[i] = 2*(1 - T.cumulativeProbability(java.lang.Math.abs(t[i])));
	return p;
    }

    public double FScore() { return ESS() / (K()-1)/(RSS()/(N()-K())); }

    public double probF() throws org.apache.commons.math.MathException {
	FDistributionImpl F = new FDistributionImpl(K()-1, N()-K());
	return 1 - F.cumulativeProbability(FScore());
    } 
    
    // Model Selection //////////////////////////////////////////////////////////////

    public double R2adjusted() {
	if(m_multi != null) return m_multi.calculateAdjustedRSquared();
	else                return 1 - (RSS()/(N()-K()))/(TSS()/(N()-1));
    }
    public double lnAPC() { 
	return java.lang.Math.log((RSS()/(N()-K()))/(1+K()/N())); }
    public double lnAIC() { 
	return 2*K()/N() + java.lang.Math.log(RSS()/N()); }
    public double lnBIC() { 
	return java.lang.Math.log(RSS()/N()) + K()/N()*java.lang.Math.log(N()); }

    // Model Description ////////////////////////////////////////////////////////////

    public double TheilsMeasure() {
	// We compute the current R^2, then compute R^2 for one missing factor in turn
	// ans = R^2 - sum(R^2 - Rj^2) for j = 1..K-1

	if(K() <= 2) return 0;  // need more cols

	double RS = 0;  // accumulate sum of deficient R^2's 
	for(int j = 1; j < K(); j++)  {  // skip ones column.
	    OLSRegression ols = new OLSRegression(m_y, data_minus_col(j));
	    RS += ols.R2();
	}
	return (2 - K())*R2() + RS;   // Notice: 1-(K-1) = 2-K
    }

    // Private Methods //////////////////////////////////////////////////////////////

    private void PrintArray(double[][] x) {
	for(int i = 0; i < x.length; i++) {
	    for(int j = 0; j < x[i].length; j++)
		System.out.print(x[i][j] + ", ");
	    System.out.println();
	}
    }

    private double getSum(double[] x) {
	double s = 0;
	for(int i = 0; i < x.length; i++) s += x[i];
	return s;
    }

    private double getSumSquares(double[] x) {
	double s = 0;
	for(int i = 0; i < x.length; i++) s += x[i]*x[i];
	return s;
    }

    private double[][] in_mult(double[][] x, double c) {  // in-place multiply
	for(int i = 0; i < x.length; i++)
	    for(int j = 0; j < x[i].length; j++)
		x[i][j] *= c;
	return x;                                        // convenience
    }

    private double[] div(double [] x, double[] y) { // parallel
	double[] D = new double[x.length];
	for(int i = 0; i < x.length; i++) D[i] = x[i] / y[i];
	return D;
    }

    private double[][] data_minus_col(int e) {
	// only call if multi-regression.
	double[][] x = m_multi_x;   
	int n = x.length;
	int k = x[0].length;
	double[][] z = new double[n][k-1];
	for(int i = 0; i < n; i++)
	    for(int j = 0; j < k; j++) {
		if(j < e) z[i][j  ] = x[i][j];
		if(j > e) z[i][j-1] = x[i][j];     
	    }
	return z;
    }

    private double[][] correlation(double[][] x) {
	// Assume x ~ R[N,K]
	// return R[K,K]
	int N = x.length;
	int K = x[0].length;
	double[][] y = standardize(x);
	double[][] cor = new double[K][K];
	for(int j = 0; j < K; j++) {
	    for(int k = 0; k < K; k++) {
		cor[j][k] = 0;
		for(int i = 0; i < N; i++) cor[j][k] += y[i][j] * y[i][k];
		cor[j][k] /= N-1;
	    }
	}
	//PrintArray(cor); // DEBUG
	return cor;
    }

    private double[][] standardize(double[][] x) {
	// Assume x ~ R[N,K]
	int N = x.length;
	int K = x[0].length;
	double[] mean = new double[K];
	for(int j = 0; j < K; j++) {
	    mean[j] = 0;
	    for(int i = 0; i < N; i++) mean[j] += x[i][j];
	    mean[j] /= N;
	}
	double[][] y = new double[N][K];
	for(int j = 0; j < K; j++)
	    for(int i = 0; i < N; i++)
		y[i][j] = x[i][j] - mean[j];
	double[] stdev = new double[K];
	for(int j = 0; j < K; j++) {
	    stdev[j] = 0;
	    for(int i = 0; i < N; i++) stdev[j] += y[i][j]*y[i][j];
	    stdev[j] /= (N-1);
	    stdev[j] = java.lang.Math.sqrt(stdev[j]);
	}
	for(int j = 0; j < K; j++)
	    for(int i = 0; i < N; i++)
		y[i][j] /= stdev[j];
	//PrintArray(y); // DEBUG
	return y;
    }

} // end class
