#include <octave/oct.h>
#include <stdio.h>
#include <stdlib.h>
#include <glpk.h>
#include <octave/ov-struct.h>

DEFUN_DLD(optimalNStorage, args, nargout, "filename, Set, Param")
{
	octave_value_list retval;
	int nargin = args.length();
	if (nargin != 3)
		print_usage();
	else
	{
		charMatrix name = args(0).char_matrix_value();
		octave_scalar_map Set = args(1).scalar_map_value();
		octave_scalar_map Param = args(2).scalar_map_value();
		
		int32NDArray set_fin = Set.contents("F").array_value();
		int numN = Set.contents("I").array_value().nelem();
		int numSfin = Set.contents("F").array_value().nelem();
		int numSinf = Set.contents("E").array_value().nelem();
		int numS = numSfin+numSinf;
		
		FloatNDArray g = Param.contents("g").array_value();
		FloatNDArray r = Param.contents("r").array_value();
		float T = Param.contents("T").float_value();
		FloatNDArray Qmax = Param.contents("Qmax").array_value();
		FloatNDArray Qmin = Param.contents("Qmin").array_value();
		FloatNDArray q0 = Param.contents("q0").array_value();
		FloatNDArray C = Param.contents("C").array_value();
		FloatNDArray D = Param.contents("D").array_value();
		FloatNDArray nul = Param.contents("nul").row_vector_value();
		FloatNDArray nuc = Param.contents("nuc").row_vector_value();
		FloatNDArray nud = Param.contents("nud").row_vector_value();
		FloatNDArray pg = Param.contents("pg").array_value();
		FloatNDArray pr = Param.contents("pr").array_value();
		FloatNDArray pc = Param.contents("pc").array_value();
		FloatNDArray pd = Param.contents("pd").array_value();
		
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
		
		FloatNDArray q(dim_vector(numSfin, numN));
		Matrix uc(numS, numN);
		Matrix ud(numS, numN);
		FloatNDArray cost(dim_vector(numN, 1), 0);
		
		for (int l=0; l<numN; l++)
		{
			for (int k=0; k<numS; k++)
			{
				uc(k, l) = glp_get_col_prim(lp, numS*l+k+1);
				ud(k, l) = glp_get_col_prim(lp, numS*l+k+1+numS*numN);
				
				cost(l) = cost(l)+uc(k, l)*pc(k, l)+ud(k, l)*pd(k, l);
			}
			cost(l) = cost(l)+g(l)*pg(l)+r(l)*pr(l);
		}
		
		for (int m=0; m<numSfin; m++)
		{
			q(m, 0) = q0(m);
			
			for(int k=1; k<numN; k++)
			{
				q(m, k) = nul(m)*q(m, k-1)+T*(nuc(m)*uc(m, k)-(1/nud(m))*ud(m, k));
			}
		}
		
		retval(0) = octave_value(q);
		retval(1) = octave_value(uc);
		retval(2) = octave_value(ud);
		retval(3) = octave_value(cost);
		
		glp_mpl_free_wksp(tran);
		glp_delete_prob(lp);
	}
return retval;
}
