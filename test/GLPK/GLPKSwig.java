import org.gnu.glpk.GLPK;
import org.gnu.glpk.GLPKConstants;
import org.gnu.glpk.glp_iocp;
import org.gnu.glpk.glp_prob;
import org.gnu.glpk.glp_tran;

public class GLPKSwig
{
	static
	{ 
		try
		{
		      // try to load Linux library
			System.loadLibrary("glpk_java"); 
		}
		catch (UnsatisfiedLinkError e)
		{
			// try to load Windows library
			System.loadLibrary("glpk_4_40_java"); 
		}
	}

	public static void main(String[] arg)
	{
		glp_prob lp;
		glp_tran tran;
		glp_iocp iocp;
		
		String fname;
		int skip = 0;
		int ret;
		
		if ( 1 != arg.length ) {
			System.out.println( "Usage: java GLPKSwig model.mod" );
			return;
		}
		
		fname = new String(arg[0]);
		
		lp = GLPK.glp_create_prob();
		System.out.println("Problem created");
		
		tran = GLPK.glp_mpl_alloc_wksp();
		ret = GLPK.glp_mpl_read_model(tran, fname, skip);
		if (ret != 0)
		{
			GLPK.glp_mpl_free_wksp(tran);
			GLPK.glp_delete_prob(lp);
			throw new RuntimeException( "Model file not found: " + fname );
		}
		
		// generate model
		GLPK.glp_mpl_generate(tran, null);
		// build model
		GLPK.glp_mpl_build_prob(tran, lp);
		// set solver parameters
		iocp = new glp_iocp();
		GLPK.glp_init_iocp(iocp);
		iocp.setPresolve(GLPKConstants.GLP_ON);
		// solve model
		ret = GLPK.glp_intopt(lp, iocp);
	    // postsolve model	
		if (ret == 0)
		{
			GLPK.glp_mpl_postsolve(tran, lp, GLPKConstants.GLP_MIP);
		}
		// free memory
		GLPK.glp_mpl_free_wksp(tran);
		GLPK.glp_delete_prob(lp);
	}
}
