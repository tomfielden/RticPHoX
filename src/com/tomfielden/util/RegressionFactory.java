package com.tomfielden.util;

import java.util.*;
import java.io.*;
import com.tomfielden.util.*;
import com.tomfielden.util.Variant.VTYPE;
import com.tomfielden.regression.*;

import org.apache.commons.math.analysis.interpolation.*;

public abstract class RegressionFactory {

    public static OLSRegression OLS(double[] Y, double[] X) {
	return OLS(Variant.createDoubleArray(Y), Variant.createDoubleArray(X));
    }
    public static OLSRegression OLS(double[] Y, double[][] X) {
	return OLS(Variant.createDoubleArray(Y), Variant.createDoubleArray(X));
    }
    public static OLSRegression OLS(Variant Y, Variant X) {
	double[] y = Y.getDoubleArray();
	switch(X.type()) {
	case ARRAY:
	case INDEX:
	    Variant shape = X.shape();
	    if(shape.size() == 1) {
		double[] x = X.getDoubleArray();
		return new OLSRegression(y,x);
	    } else {
		double[][] x = new double[X.size()][];
		for(int i = 0; i < x.length; i++) 
		    x[i] = X.get(i).getDoubleArray();
		return new OLSRegression(y,transpose(x));
	    }
	default:
	    System.out.println("RegressionFactory:OLS("+Y.type()+", "+X.type()+") unsupported");
	}
        return null;
    }

    public static Variant Loess(double[] x, double[] y) 
	throws org.apache.commons.math.MathException {
	return Loess(Variant.createDoubleArray(x), Variant.createDoubleArray(y));
    }
    public static Variant Loess(Variant x, Variant y) 
	throws 	org.apache.commons.math.MathException {
	LoessInterpolator loess = new LoessInterpolator();
	Variant i = x.isort();
	double[] sorted_x = x.get(i).getDoubleArray();
	double[] sorted_y = y.get(i).getDoubleArray();
	double[] z = loess.smooth(sorted_x, sorted_y);
        return Variant.createDoubleArray(z);
    }

    /// Private Utilities ///

    private static double[][] transpose(double[][] x) {
	int R = x.length;
	int C = x[0].length;
	double[][] data = new double[C][R];
	for(int i = 0; i < R; i++) 
	    for(int j = 0; j < C; j++)
		data[j][i] = x[i][j];
	return data;
    }

} // end class

