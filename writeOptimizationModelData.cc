#include <octave/oct.h>
#include <iostream>
#include <fstream>
#include <sstream>

using namespace std;

DEFUN_DLD(writeOptimizationModelData, args, nargout, "Writes the optimization model data")
{
	octave_value_list retval;
	int k, m, nargin = args.length();
	charMatrix filename = args(0).char_matrix_value();
	int N = args(1).int_value();

	octave_idx_type len = filename.nelem();
	charMatrix suffix = ".dat";
	filename.insert(suffix, 0, len);
	ofstream myfile;
	myfile.open (filename.row_as_string(0).c_str());
	ostringstream s;

	RowVector n(N);
	for (k=0; k<N; k++)
	{
		n.elem(k) = k+1;
	}

	s << "data;\n\n";
	s << "set I :=" << n << ";\n";

	for(k=2; k<nargin-1; k=k+2)
	{
		s << "param ";
		charMatrix name = args(k).char_matrix_value();
		s << name.row_as_string(0).c_str() << " :=";

		RowVector fVec = args(k+1).row_vector_value();

		if(fVec.nelem() == 1)
		{
			s << " " << fVec.elem(0);
		}
		else
		{
			for (m=0; m<N; m++)
			{
				s << " " << m+1 << " " << fVec.elem(m);
			}
		}
		s << ";\n";
	}

	s << "\n" << "end;\n\n";
	myfile << s.str();
	myfile.close();
	return retval;
}
