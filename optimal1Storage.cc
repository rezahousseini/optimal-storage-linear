#include <octave/oct.h>
#include <stdio.h>
#include <stdlib.h>
#include <glpk.h>

DEFUN_DLD(optimal1StorageV2_2, args, nargout, "optimal1StorageV2_2")
{
	octave_value_list retval;
	int nargin = args.length();
	if (nargin != 1)
		print_usage();
	else
	{
		charMatrix name = args(0).char_matrix_value();

		charMatrix suffixDat = ".dat";
		charMatrix suffixMod = ".mod";
		charMatrix suffixMps = ".mps";
		charMatrix suffixSol = ".sol";
		octave_idx_type len = name.nelem();

		glp_prob *lp;
		glp_tran *tran;
		int ret, N, k;

		lp = glp_create_prob();
		tran = glp_mpl_alloc_wksp();

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

		glp_simplex(lp, NULL);

//		name.insert(suffixSol, 0, len);
//		ret = glp_print_sol(lp, name.row_as_string(0).c_str());
//		if (ret != 0)
//		{
//			fprintf(stderr, "Error on writing MPS file\n");
//		}

		skip: 
		N = glp_get_num_cols(lp)/5;
		RowVector uc(N);
		RowVector ud(N);
		RowVector gs(N);
		RowVector rs(N);
		RowVector Q(N);
		for(k=1; k<=N; k++)
		{
			uc.elem(k-1) = glp_get_col_prim(lp, k);
			ud.elem(k-1) = glp_get_col_prim(lp, k+N);
			gs.elem(k-1) = glp_get_col_prim(lp, k+2*N);
			rs.elem(k-1) = glp_get_col_prim(lp, k+3*N);
			Q.elem(k-1) = glp_get_col_prim(lp, k+4*N);
		}

		retval(0) = octave_value(uc);
		retval(1) = octave_value(ud);
		retval(2) = octave_value(gs);
		retval(3) = octave_value(rs);
		retval(4) = octave_value(Q);

		glp_mpl_free_wksp(tran);
		glp_delete_prob(lp);
	}
return retval;
}
