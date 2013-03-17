package com.tomfielden.plotting;

import java.util.*;
import java.io.*;
import java.text.*;

import com.panayotis.gnuplot.JavaPlot;
import com.panayotis.gnuplot.terminal.FileTerminal;
import com.panayotis.gnuplot.plot.DataSetPlot;
import com.panayotis.gnuplot.style.PlotStyle;
import com.panayotis.gnuplot.style.Style;

/*
  This class is a wrapper for JavaPlot. Only this class and it's progeny know about it.
 */

public class SimplePlot {

    private JavaPlot m_plot;

    public SimplePlot() {
	m_plot = new JavaPlot();
    }

    public void plot(double[][] data, String style, String title) {
        DataSetPlot ds = new DataSetPlot(data);  // So we can retrieve style?
	PlotStyle stl = ds.getPlotStyle();
        setStyle(stl, style);
	ds.setTitle(title); 
        m_plot.addPlot(ds);
    }

    public void xrange(double min, double max) {
	m_plot.getAxis("x").setBoundaries(min, max);
    }

    public void yrange(double min, double max) {
	m_plot.getAxis("y").setBoundaries(min, max);
    }

    public void show() {
	m_plot.plot();
    }

    public void save(String filename) {
	m_plot.setTerminal(new FileTerminal("png", filename));
	m_plot.plot();
    }

    /// Private ///

    private static void setStyle(PlotStyle stl, String style) {
	CharacterIterator it = new StringCharacterIterator(style);

	// Iterate over the characters in the forward direction
	for (char ch=it.first(); ch != CharacterIterator.DONE; ch=it.next()) 
	    setStyle(stl, ch);
    }

    private static void setStyle(PlotStyle stl, char style) {
	switch(style) {
	case '-':  stl.setStyle(Style.LINES); break;
	case '.':  stl.setStyle(Style.DOTS); break;
	case '+':  stl.setStyle(Style.POINTS); break;
	case 'L':  stl.setStyle(Style.LINESPOINTS); break;
	case 'i':  stl.setStyle(Style.IMPULSES); break;
	case 'b':  stl.setStyle(Style.BOXES); break;
	case 'f':  stl.setStyle(Style.FSTEPS); break;
	case 'h':  stl.setStyle(Style.HISTEPS); break;
	case 's':  stl.setStyle(Style.STEPS); break;
	}
    }
    
} // end class
