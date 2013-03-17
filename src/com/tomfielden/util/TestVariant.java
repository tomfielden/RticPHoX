package com.tomfielden.util;                                                 

import junit.framework.*;
//import com.tomfielden.util.Variant;

public class TestVariant extends TestCase {

    public void testEquals() {
	Variant A1 = new Variant();
	Variant B1 = new Variant();
	assertTrue("NULL type not set", A1.type() == Variant.VTYPE.NULL);
	assertTrue("NULL types not equal", A1.equals(B1));

	Variant A2 = new Variant(true);
	Variant B2 = new Variant(true);
	assertTrue("Expected BOOL", A2.type() == Variant.VTYPE.BOOL);
	assertTrue("BOOL types not equal", A2.equals(B2));

	Variant A3 = new Variant(7);
	Variant B3 = new Variant(7);
	assertTrue("Expected INT", A3.type() == Variant.VTYPE.INT);
	assertTrue("INTs not equal: "+A3+", "+B3, A3.equals(B3));

	Variant A4 = new Variant(3.14);
	Variant B4 = new Variant(3.14);
	assertTrue("Expected DOUBLE", A4.type() == Variant.VTYPE.DOUBLE);
	assertTrue("DOUBLEs not equal", A4.equals(B4));

	Variant A5 = new Variant("Now is the time");
	Variant B5 = new Variant("Now is the time");
	assertTrue("Expected STRING", A5.type() == Variant.VTYPE.STRING);
	assertTrue("STRINGs not equal", A5.equals(B5));

	Variant A6 = new Variant(new Variant[]{A1,A2,A3,A4,A5});
	Variant B6 = new Variant(new Variant[]{B1,B2,B3,B4,B5});
	assertTrue("Expected ARRAY", A6.type() == Variant.VTYPE.ARRAY);
	assertTrue("ARRAYs not equal", A6.equals(B6));

	Variant A7 = Variant.createMask(new Boolean[]{true, false, true});
	Variant B7 = Variant.createMask(new Boolean[]{true, false, true});
	assertTrue("Expected MASK", A7.type() == Variant.VTYPE.MASK);
	assertTrue("MASKs not equal: "+A7+", "+B7, A7.equals(B7));

	Variant A8 = Variant.createSlice("1:2:-3");
	Variant B8 = Variant.createSlice("1:2:-3");
	assertTrue("Expected SLICE", A8.type() == Variant.VTYPE.SLICE);
	assertTrue("SLICEs not equal", A8.equals(B8));

	Variant A9 = Variant.createIntSet(new Variant(new Integer[]{1,2,3}));
	Variant B9 = Variant.createIntSet(new Variant(new Integer[]{1,2,3}));
	assertTrue("Expected INT_SET", A9.type() == Variant.VTYPE.INT_SET);
	assertTrue("INT_SETs not equal", A9.equals(B9));

	Variant A10 = Variant.createStringSet(new Variant(new String[]{"foo","bar"}));
	Variant B10 = Variant.createStringSet(new Variant(new String[]{"foo","bar"}));
	assertTrue("Expected STRING_SET", A10.type() == Variant.VTYPE.STRING_SET);
	assertTrue("STRING_SETs not equal", A10.equals(B10));

	Variant A11 = Variant.createIntMap(new Variant(new Integer[]{1,2,3}), 
					   new Variant(new Variant[]{A2,A3,A4}));
	Variant B11 = Variant.createIntMap(new Variant(new Integer[]{1,2,3}), 
					   new Variant(new Variant[]{A2,A3,A4}));
	assertTrue("Expected INT_MAP", A11.type() == Variant.VTYPE.INT_MAP);
	assertTrue("INT_MAPSs not equal", A11.equals(B11));

	Variant A12 = Variant.createStringMap(new Variant(new String[]{"A","B","C"}),
					      new Variant(new Variant[]{A8,A9,A10}));
	Variant B12 = Variant.createStringMap(new Variant(new String[]{"A","B","C"}),
					      new Variant(new Variant[]{A8,A9,A10}));
	assertTrue("Expected STRING_MAP", A12.type() == Variant.VTYPE.STRING_MAP);
	assertTrue("STRING_MAPSs not equal", A12.equals(B12));
    }

    public void testCreateMask() {
	Variant A1 = Variant.createMask(true);
	assertTrue("Expected MASK", A1.type() == Variant.VTYPE.MASK);
	assertTrue("Expected unit sized MASK: "+A1, A1.size() == 1);
	assertTrue("Expected truth: " + A1, A1.get(0).isTrue());

	Variant A2 = Variant.createMask(5);
	assertTrue("Expected MASK", A2.type() == Variant.VTYPE.MASK);
	assertTrue("Expected sizable MASK: "+A2, A2.size() == 5);
	assertTrue("Expected truth: "+A2, A2.sum().equals(5));

	Variant A3 = Variant.createMask(new Variant(new Integer[]{1,0,1,1}));
	assertTrue("Expected MASK", A3.type() == Variant.VTYPE.MASK);
	assertTrue("Expected unit sized MASK: "+A3, A3.size() == 4);
	assertTrue("Expected some truth: "+A3, A3.sum().equals(3));

	/*
	>>> D = Variant({'Region':('CA','CA','CA','WCI','WCI','WCI'), 
	                 'Sector':('Cement','CMM','Forestry','CMM','Forestry', 'Agriculture')})
	*/

	Variant R = new Variant(new String[]{"CA","CA", "CA", "WCI", "WCI", "WCI"});
	Variant S = new Variant(new String[]{
		"Cement","CMM", "Forestry", "CMM", "Forestry", "Agriculture"});
	Variant D = Variant.createStringMap(new Variant(new String[]{"Region", "Sector"}),
					    new Variant(new Variant[]{R, S}));

	// >>> D[{'Sector': ('Cement', 'Forestry')}]
	Variant X1 = new Variant(new String[]{"Cement", "Forestry"});
	Variant S1 = Variant.createStringMap(new Variant(new String[]{"Sector"}),
					     new Variant(new Variant[]{X1}));
	Variant M1 = D.get(S1);
	Variant N1 = Variant.createMask(new Boolean[]{true, false, true, false, true, false});
	assertTrue("Select 1: "+D+"\n"+S1+"\n"+M1, M1.equals(N1));

	// >>> D[{'Sector': ('Cement', 'Forestry'), 'Region':'WCI'}]
	Variant X2 = new Variant(new String[]{"Cement", "Forestry"});
	Variant S2 = Variant.createStringMap(new Variant(new String[]{"Sector", "Region"}),
					     new Variant(new Variant[]{X1, new Variant("WCI")}));
	Variant M2 = D.get(S2);
	Variant N2 = Variant.createMask(new Boolean[]{true, false, true, true, true, true});
	assertTrue("Select 2: "+D+"\n"+S2+"\n"+M2, M2.equals(N2));
    }

    public void testShape() {
	Variant a123  = new Variant(new Integer[]{1,2,3});
	Variant a2    = new Variant(new Integer[]{2});
	Variant a3    = new Variant(new Integer[]{3});
	Variant a4    = new Variant(new Integer[]{4});
	Variant a23   = new Variant(new Integer[]{2,3});
	Variant a223  = new Variant(new Integer[]{2,2,3});
	Variant a45   = new Variant(new Integer[]{4,5});
	Variant a456  = new Variant(new Integer[]{4,5,6});
	Variant a4567 = new Variant(new Integer[]{4,5,6,7});
	Variant b32   = new Variant(new Variant[]{a123,a45});
	Variant b33   = new Variant(new Variant[]{a123,a456});
	Variant c     = new Variant(new Variant[]{b32, a4567});
	Variant c32   = new Variant(new Variant[]{a3,a2});
	Variant d     = new Variant(new Variant[]{b33, a4567});
	Variant e     = new Variant(new Variant[]{b33, b33});

	Variant Q1 = a123.shape();                            // [1,2,3]
	Variant A1 = a3.toArrayDeeply();                      // (3)
	assertTrue("Shape 1: "+Q1+" != "+A1, Q1.equals(A1));  
	Variant Q2 = b32.shape();                             // ((1,2,3),(4,5))
	Variant A2 = c32.toArrayDeeply();                     // ((3),(2))
	assertTrue("Shape 2: "+Q2+" != "+A2, Q2.equals(A2));  
	Variant Q3 = b33.shape();                             // ((1,2,3),(4,5,6))
	Variant A3 = a23.toArrayDeeply();                     // (2, 3)
	assertTrue("Shape 3: "+Q3+" != "+A3, Q3.equals(A3));  
	Variant Q4 = c.shape();                               // (((1,2,3),(4,5)), (4,5,6,7))
	Variant A4 = (new Variant(new Variant[]{c32, a4})).toArrayDeeply(); // (((3), (2)), (4))
	assertTrue("Shape 4: "+Q4+" != "+A4, Q4.equals(A4));  
	Variant Q5 = d.shape();                               // (((1,2,3),(4,5,6)), (4,5,6,7))
	Variant A5 = (new Variant(new Variant[]{a23, a4})).toArrayDeeply(); // ((2,3), (4))
	assertTrue("Shape 5: "+Q5+" != "+A5, Q5.equals(A5));  
	Variant Q6 = e.shape();                          // (((1,2,3),(4,5,6)),((1,2,3),(4,5,6)))
	Variant A6 = a223.toArrayDeeply();                    // (2,2,3)
	assertTrue("Shape 6: "+Q6+" != "+A6, Q6.equals(A6));  
    }

    public void testAdd() {
        Variant num1  = new Variant(3);
        Variant num2  = new Variant(2);
        Variant total = new Variant(5);
        Variant sum   = Variant.add(num1, num2);
        assertTrue(Variant.cmp(sum, total) == 0);
    }
    
    public void testSubtract() {
        assertTrue("This is supposed to pass", true);
        //assertTrue("This is supposed to fail", false);
    }
    
}