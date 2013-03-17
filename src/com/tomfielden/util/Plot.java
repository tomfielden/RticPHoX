package com.tomfielden.util;

import org.python.core.*;

import java.util.*;
import java.io.*;
import java.text.*;

import com.tomfielden.plotting.*;
import com.tomfielden.util.Variant.VTYPE;

/*
  This class knows about Variants and Jython so that .plotting doesn't have to.
 */

public class Plot {

    private static final String default_style = "-";
    private static final String default_title = "data";

    private static final char CUMULATIVE_STYLE = 'c';

    private SimplePlot m_plot;

    public Plot() { 
	m_plot = PlotFactory.createSimplePlot();
    }

    /* I think this collides with the PyObject version
    public Plot(Variant Y) {
	m_plot = PlotFactory.createSimplePlot();
	show(Y);
    }
    */

    public Plot(PyObject X) {
	m_plot = PlotFactory.createSimplePlot();
	show(new Variant(X));
    }

    public Plot(PyObject Y, String style) {
	m_plot = PlotFactory.createSimplePlot();
	show(new Variant(Y), style);
    }

    public Plot(PyObject X, PyObject Y) {
	m_plot = PlotFactory.createSimplePlot();
	show(new Variant(X), new Variant(Y));
    }

    public Plot(PyObject X, PyObject Y, String style) {
	m_plot = PlotFactory.createSimplePlot();
	show(new Variant(X), new Variant(Y), style);
    }

    // Methods ///////////////////////////////////////////////

    public Plot plot(Variant Y) {
	return plot(Y,default_style);
    }

    public Plot plot(Variant Y, String style) {
	return plot(Y, style, default_title);
    }

    // workhorse //
    public Plot plot(Variant Y, String style, String title) {
	switch(Y.type()) {
	case ARRAY:
	    if(Y.shape().size() == 1 && Y.isNumbers()) {
		int     N = Y.size();
		Variant X = Variant.linspace(1, N, N);
		plot(X,Y,style,title);
	    } else {  // peel-off first layer
		for(int i = 0; i < Y.size(); i++)
		    plot(Y.get(i), style);
	    }
	    break;
	case STRING_MAP:
	    Variant K = Y.keys().toArray();
	    for(int i = 0; i < K.size(); i++) {
		String new_title = K.get(i).toString();
		plot(Y.get(K.get(i)), style, new_title);
	    }
	    break;
	case LRV:
	case NRV:
	case NRS:  // XXX
	    if(hasStyle(style, CUMULATIVE_STYLE)) {
		if(style.length() == 1) style = default_style;
		double[][] M = Y.getCumulativeData();
		double[] cx = M[0];
		double[] cy = M[1];
		plot(cx,cy,style,title+" (Cumulative)");
	    } else {
		double[][] M = Y.getMatrix();
		double[] dx = M[0];
		double[] dy = M[1];
		double[] cx = M[2];
		double[] cy = M[3];
		
		if(dx.length >  0 && cx.length == 0) 
		    plot(dx, dy,"i",title);
		else if(dx.length == 0 && cx.length >  0) 
		    plot(cx,cy,style,title);
		else if(dx.length >  0 && cx.length >  0) {
		    plot(dx,dy,"i",title+" (Discrete)");
		    plot(cx,cy,style,title+" (Continuous)");
		} else 
		    System.out.println("No graphs found to plot.");
	    }
	}
	return this;
    }

    public Plot plot(Variant X, Variant Y) {
	return plot(X,Y,default_style,default_title);
    }

    public Plot plot(Variant X, Variant Y, String style) {
	return plot(X,Y,style,default_title);
    }

    public Plot plot(Variant X, Variant Y, String style, String title) {
	// Assume X & Y are parallel
	double[]   x = X.getDoubleArray();
	double[]   y = Y.getDoubleArray();
	return plot(x, y, style, title);
    }

    public Plot plot(double[] x, double[] y, String style, String title) {
	// Assume x & y are parallel
	double[][] data = transpose(x,y);
	m_plot.plot(data, style, title);
	return this;
    }

    public Plot xrange(double min, double max) {
	m_plot.xrange(min, max);
	return this;
    }

    public Plot yrange(double min, double max) {
	m_plot.yrange(min, max);
	return this;
    }

    public void show()                                   { m_plot.show(); }
    public void show(Variant Y)                          { plot(Y);           show(); }
    public void show(Variant Y, String style)            { plot(Y,   style);  show(); }
    public void show(Variant X, Variant Y)               { plot(X, Y);        show(); }
    public void show(Variant X, Variant Y, String style) { plot(X, Y, style); show(); }

    public void save(String filename)                    { m_plot.save(filename); }

    /// Private Utilities ///////////////////////////////////////////////////////////////////

    private static double[][] transpose(double[] x, double[] y) {
	int N = x.length;
	double[][] data = new double[N][2];
	for(int i = 0; i < N; i++) {
	    data[i][0] = x[i];
	    data[i][1] = y[i];
	}
	return data;
    }

    private static double[] rshift(double[] x) {
	double[] y = new double[x.length];
	if(x.length == 0) return y;
	y[0] = 0;
	for(int i = 1; i < x.length; i++) y[i] = x[i-1];
	return y;
    }

    private static boolean hasStyle(String style, char s) {
	CharacterIterator it = new StringCharacterIterator(style);

	// Iterate over the characters in the forward direction
	for (char ch=it.first(); ch != CharacterIterator.DONE; ch=it.next()) 
	    if(ch == s) return true;
	return false;
    }

} // end class

