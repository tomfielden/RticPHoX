/* This is the minimal embedding of Jython into Java
   All we need is for CLASSPATH to include jython.jar
   and the ability of the user to write to /usr/share/java/cachedir/packages
 */

//import java.util.Iterator;

import org.python.core.*;
import org.python.util.*;

public class JythonMinimal {
    public static void main(String[] argv) {
	PySystemState.initialize();
	PythonInterpreter interp = new PythonInterpreter();
	PyObject value = interp.eval("2+2");
	System.out.println(value);

	// beyond minimal 
	PyString s = new PyString("Foo String");
	//String t = s.asString();  // sure enough this doesn't exist
	//String t = (String)s;       // instead use a cast? Nope.
	String t = s.toString();    // instead use the usual toString().
	System.out.println("PyString: " + s);      // This should work normally.
	System.out.println("  String: " + t);

	// Variant
	PyList list = new PyList();
	list.add(4);
	list.add("foo");
	list.add(5.2);
	Object p = list.get(1);
	System.out.println("PyList: "+p);

	// Variant needs to set the dictionary.
	PyDictionary dict = new PyDictionary();
	dict.put("foo","bar");
	dict.put(3,72);
	PyList items = dict.items();
	System.out.println(items);
	for(int i = 0; i < items.size(); i++) {
	    PyTuple pair = (PyTuple)(items.get(i));
	    System.out.println(pair.get(0) + ", " + pair.get(1));
	}
    }
}