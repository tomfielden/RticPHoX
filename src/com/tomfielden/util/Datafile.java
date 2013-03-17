package com.tomfielden.util;

import java.util.*;
import java.io.*;

// OBSOLETE. USE Lookuptable.

/* Assume first row contains column names separated by whitespace
 * Text may be quoted to include embedded whitespace
 * No handling of comments within datafile
 *
 */

public class Datafile {

    private HashMap<String, Variant> m_columns =  new HashMap<String, Variant>();

    public Datafile(String filename) {
	try {
	    BufferedReader in = new BufferedReader(new FileReader(filename));
		
	    String str;
	    ArrayList<ArrayList<String>> array = new ArrayList<ArrayList<String>>();
		
	    while((str = in.readLine()) != null) {
		ArrayList<String> tokens = parseLine(str);
		for(int i = 0; i < tokens.size(); i++) {
		    if(array.size() == 0) 
			for(int j = 0; j < tokens.size(); j++)
			    array.add(new ArrayList<String>());
		    array.get(i).add(tokens.get(i));
		}
	    }
	    
	    for(ArrayList<String> column: array) {
		if(column.size() == 0) break;
		m_columns.put(column.get(0), getVariant(column));
	    }
	    
	} catch(java.io.FileNotFoundException ex) {
	    System.out.println("Must specify a valid file: " + ex);
	} catch(java.io.IOException ex) {
	    System.out.println("The file is messed up: " + ex);
	}
    }

    public Variant getColumnMap() {return Variant.createStringMap(m_columns);}

    ///////////////////////////////////////////////////////////////////
    // Private
    
    private Variant getVariant(ArrayList<String> column) {
	// ignore first item (header): Guaranteed to be present.
	if(column.size() == 1) return new Variant();
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
		return new Variant(scolumn);
	    }
	}
	return Variant.createDoubleArray(dcolumn);   // Python wiped out the simple version.
    }
    
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
    
    /////////////////////////////////////////////////////////////////////////////////////
} // End Class
