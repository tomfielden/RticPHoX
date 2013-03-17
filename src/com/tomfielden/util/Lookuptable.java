package com.tomfielden.util;

import java.util.*;
import java.io.*;

/*  The job of this class is to parse a file into a rectangular array of strings.
 *  Various formattings are supported to convert the array to useful structures.
 */

public class Lookuptable {

    public ArrayList<ArrayList<String>> m_array = new ArrayList<ArrayList<String>>(); // (<col>,<row>)

    /////////////////////////////////////
    // Constructors

    public Lookuptable(String filename) {
	try {
	    BufferedReader in = new BufferedReader(new FileReader(filename));
		
	    String str;
	    while((str = in.readLine()) != null) {
		ArrayList<String> tokens = parseLine(str);
		for(int i = 0; i < tokens.size(); i++) {
		    if(m_array.size() == 0) 
			for(int j = 0; j < tokens.size(); j++)
			    m_array.add(new ArrayList<String>());
		    m_array.get(i).add(tokens.get(i));
		}
	    }
	    
	} catch(java.io.FileNotFoundException ex) {
	    System.out.println("Must specify a valid file: " + ex);
	} catch(java.io.IOException ex) {
	    System.out.println("The file is messed up: " + ex);
	}
    }

    public Variant get2DColumnLookup() {
	Variant lookup = Variant.createStringMap();
	if(m_array.size() <= 1) return lookup;

	ArrayList<String>    col_0       = m_array.get(0);
	Map<String, Variant> lookup_map  = lookup.toStringMap();

	for(int i = 1; i < m_array.size(); i++) {
	    ArrayList<String>    col_i   = m_array.get(i);
	    Variant              rows    = Variant.createStringMap();
	    Map<String, Variant> rows_map = rows.toStringMap();
	    lookup_map.put(col_i.get(0), rows);
	    for(int j = 1; j < col_i.size(); j++) 
		rows_map.put(col_0.get(j), getVariant( col_i.get(j) ));
	}

	return lookup;
    }

    public Variant getColumnMap() {
	Variant string_map = Variant.createStringMap();
	for(int i = 0; i < m_array.size(); i++) setColumn(string_map, m_array.get(i));
	return string_map;
    }

    ///////////////////////////////////////////////////////////////////////////////
    // Private File Parsing

    private boolean regular_mode = true;
    private static final String fWHITESPACE_AND_QUOTES = ",\t\r\n\"";
    private static final String fQUOTES_ONLY ="\"";
    
    private ArrayList<String> parseLine(String sLine) {
	ArrayList<String> result = new ArrayList<String>();
	
	boolean returnTokens = true;
	String currentDelims = fWHITESPACE_AND_QUOTES;
	StringTokenizer parser = new StringTokenizer(sLine, currentDelims, returnTokens);
	
	String token = null;
	while ( parser.hasMoreTokens() ) {
	    if (regular_mode) currentDelims = fWHITESPACE_AND_QUOTES;
	    else              currentDelims = fQUOTES_ONLY;
	    token = parser.nextToken(currentDelims);
	    if ( token.equals(fQUOTES_ONLY) )
		regular_mode = !regular_mode;
	    else 
		addNonTrivialWordToResult( token, result );
	}
	
	return result;
    }
    
    private void addNonTrivialWordToResult( String aToken, ArrayList<String> aResult ){
	if (aToken == null) return;
	aToken = aToken.trim();
	if ( aToken.equals("") || aToken.equals(",") ) return;
	aResult.add( aToken );
    }

    ///////////////////////////////////////////////////////////////////////////////
    // Private Cell Parsing
    
    private Variant getVariant(String value) {
	try {
	    Double d = Double.valueOf(value);
	    return new Variant(d);
	} catch (NumberFormatException nfe) {
	    String lower = value.toLowerCase();
	    if(lower == "true")  return new Variant(Boolean.TRUE);
	    if(lower == "false") return new Variant(Boolean.FALSE);
	}
	return new Variant(value);
    }

    private void setColumn(Variant string_map, ArrayList<String> column) {
	// First element is column key (String)
	Variant key = new Variant(column.get(0));
	double [] dcolumn = new double [column.size() - 1];
	for(int i = 1; i < column.size(); i++) {
	    String svalue = column.get(i);
	    try {
		dcolumn[i-1] = Double.valueOf(svalue);
	    } catch (NumberFormatException nfe) {
		// give up: it's a string column.
		String[] scolumn = new String[column.size() - 1];
		for(int j = 1; j < column.size(); j++) 
		    scolumn[j-1] = column.get(j);
		Variant array = new Variant(scolumn);
		string_map.set(key, array);
		return;
	    }
	}
	Variant array = Variant.createDoubleArray(dcolumn); // Python wiped out the simple version.
	string_map.set(key, array);
    }

    ///////////////////////////////////////////////////////////////////////////////
} // end Lookuptable
