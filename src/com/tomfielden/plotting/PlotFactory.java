package com.tomfielden.plotting;

import java.util.*;
import java.io.*;
import com.tomfielden.util.*;

/*
  This class provides a level of indirection that is not clear at this time.
 */

public class PlotFactory {

    public static SimplePlot createSimplePlot() {
	return new SimplePlot();
    }

} // end class
