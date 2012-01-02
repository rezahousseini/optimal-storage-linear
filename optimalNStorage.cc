#include <octave/oct.h>
#include <stdio.h>
#include <stdlib.h>
#include <glpk.h>

DEFUN_DLD(optimalNStorage, args, nargout, "filename, numN, numS")
{
	octave_value_list retval;
	int nargin = args.length();
	if (nargin != 3)
		print_usage();
	else
	{
		charMatrix name = args(0).char_matrix_value();
		int N = args(1).int_value();
		int A = args(2).int_value();

		charMatrix suffixDat = ".dat";
		charMatrix suffixMod = ".mod";
		charMatrix suffixMps = ".mps";
		charMatrix suffixSol = ".sol";
		octave_idx_type len = name.nelem();

		glp_prob *lp;
		glp_tran *tran;
		glp_smcp parm;
		int ret;

		lp = glp_create_prob();
		tran = glp_mpl_alloc_wksp();
		glp_init_smcp(&parm);

		name.insert(suffixMod, 0, len);
		ret = glp_mpl_read_model(tran, name.row_as_string(0).c_str(), 1);
		if (ret != 0)
		{
			fprintf(stderr, "Error on translating model\n");
			goto skip;
		}

		name.insert(suffixDat, 0, len);
		ret = glp_mpl_read_data(tran, name.row_as_string(0).c_str());
		if (ret != 0)
		{
			fprintf(stderr, "Error on translating data\n");
			goto skip;
		}
		
		ret = glp_mpl_generate(tran, NULL);
		if (ret != 0)
		{
			fprintf(stderr, "Error on generating model\n");
			goto skip;
		}

		glp_mpl_build_prob(tran, lp);

		name.insert(suffixMps, 0, len);
		ret = glp_write_mps(lp, GLP_MPS_FILE, NULL, name.row_as_string(0).c_str());
		if (ret != 0)
		{
			fprintf(stderr, "Error on writing MPS file\n");
		}

		parm.msg_lev = GLP_MSG_ALL;

		glp_simplex(lp, &parm);

//		name.insert(suffixSol, 0, len);
//		ret = glp_print_sol(lp, name.row_as_string(0).c_str());
//		if (ret != 0)
//		{
//			fprintf(stderr, "Error on writing MPS file\n");
//		}

		skip: 
		Matrix uc(A,N);
		Matrix ud(A,N);

		for (int l=0; l<N; l++)
		{
			for (int k=0; k<A; k++)
			{
				uc.elem(k,l) = glp_get_col_prim(lp, A*l+k+1);
				ud.elem(k,l) = glp_get_col_prim(lp, A*l+k+1+A*N);
			}
		}
		
		retval(0) = octave_value(uc);
		retval(1) = octave_value(ud);

		glp_mpl_free_wksp(tran);
		glp_delete_prob(lp);
	}
return retval;
}
