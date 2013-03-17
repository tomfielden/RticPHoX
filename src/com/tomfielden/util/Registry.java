package com.tomfielden.util;

import java.util.*;
import java.io.*;

public class Registry extends Variant {

    //////////////////////////////////////////////////////////////////////////////////////
    
    public Registry()           { super(Variant.createStringMap()); }
    public Registry(Variant R)  { super(R);}
    public Registry(Registry r) { super(r.copy()); }

    ///////////////////////////////////////////////////////////////////////////////////

} // class
