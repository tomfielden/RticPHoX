package com.tomfielden.util;

import org.python.core.*;
import org.python.util.PythonInterpreter;

import java.util.*;
import java.io.*;

public class Interpreter {


    public static void Interpret(String code) {
	try {
	    PythonInterpreter interp = new PythonInterpreter();
	    interp.exec(code);
	    System.out.println(interp.getLocals());
	} catch(Exception e) {
	    System.out.println(e);
	}
    }

    public static Variant Interpret2(Registry R, String code) {
	try {
	    PythonInterpreter interp = new PythonInterpreter();

	    // inject given variables
	    Set<String> keys = R.keys().toStringSet();
	    for(String k: keys) interp.set(k, R.get(k));

	    interp.exec(code);
	    Variant locals = new Variant(interp.getLocals());

	    return Variant.sub(locals, R);   // remove locals we put in (thru Registry).
	} catch(Exception e) {
	    System.out.println(e);
	    return new Variant();
	}
    }

    ////////////////////////////////////////////////////////////////////////////////////////////
} // class