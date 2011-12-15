#include <octave/oct.h>
#include <stdio.h>
#include <stdlib.h>
#include <glpk.h>

DEFUN_DLD(optimal1Storage, args, nargout, "filename, N")
{
	octave_value_list retval;
	int nargin = args.length();
	if (nargin != 2)
		print_usage();
	else
	{
		charMatrix name = args(0).char_matrix_value();
		int N = args(1).int_value();

		charMatrix suffixDat = ".dat";
		charMatrix suffixMod = ".mod";
		charMatrix suffixMps = ".mps";
		charMatrix suffixSol = ".sol";
		octave_idx_type len = name.nelem();

		glp_prob *lp;
		glp_tran *tran;
		glp_smcp parm;
		int ret, k, nvar, m;

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

//		name.insert(suffixMps, 0, len);
//		ret = glp_write_mps(lp, GLP_MPS_FILE, NULL, name.row_as_string(0).c_str());
//		if (ret != 0)
//		{
//			fprintf(stderr, "Error on writing MPS file\n");
//		}

		parm.msg_lev = GLP_MSG_ALL;

		glp_simplex(lp, &parm);

//		name.insert(suffixSol, 0, len);
//		ret = glp_print_sol(lp, name.row_as_string(0).c_str());
//		if (ret != 0)
//		{
//			fprintf(stderr, "Error on writing MPS file\n");
//		}

		skip: 
		nvar = glp_get_num_cols(lp)/N;
		RowVector retVec(N);
		for(m=0; m<nvar; m++)
		{
			for(k=1; k<=N; k++)
			{
				retVec.elem(k-1) = glp_get_col_prim(lp, k+m*N);
			}
			retval(m) = octave_value(retVec);
		}

		glp_mpl_free_wksp(tran);
		glp_delete_prob(lp);
	}
return retval;
}
