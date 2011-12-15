#include <octave/oct.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <octave/ov-struct.h>

using namespace std;

DEFUN_DLD(writeOptimizationModelData, args, nargout, "filename, Set, Param")
{
	octave_value_list retval;
	charMatrix filename = args(0).char_matrix_value();
	octave_scalar_map S = args(1).scalar_map_value();
	octave_scalar_map P = args(2).scalar_map_value();

	octave_idx_type len = filename.nelem();
	charMatrix suffix = ".dat";
	filename.insert(suffix, 0, len);
	ofstream myfile;
	myfile.open (filename.row_as_string(0).c_str());
	ostringstream s;
	s << "data;\n\n";

	int nS = S.nfields();
	string_vector keysS = S.fieldnames();
	int nP = P.nfields();
	string_vector keysP = P.fieldnames();

	for(int k=0; k<nS; k++)
	{
		s << "set " << keysS[k] << " :=" << S.contents(k).row_vector_value() << ";\n";
	}

	s << "\n";

	for(int k=0; k<nP; k++)
	{
		s << "param " << keysP[k];
		dim_vector dv = P.contents(k).dims();
		if(dv(0) == 1 and dv(1) == 1)
		{
			s << ":= " << P.contents(k).float_value();
		}
		else if(dv(0) > 1 and dv(1) == 1)
		{
			s << ":= ";
			for(int m=0; m<dv(0); m++)
			{
				if(P.contents(k).array_value().elem(m) == 1.0/0.0)
				{
					s << " " << m+1 << " " << -1;
				}
				else
				{
					s << " " << m+1 << " " << P.contents(k).row_vector_value().elem(m);
				}
			}
		}
		else if(dv(0) == 1 and dv(1) > 1)
		{
			s << ":= ";
			for(int m=0; m<dv(1); m++)
			{
				if(P.contents(k).row_vector_value().elem(m) == 1.0/0.0)
				{
					s << " " << m+1 << " " << -1;
				}
				else
				{
					s << " " << m+1 << " " << P.contents(k).row_vector_value().elem(m);
				}
			}
		}
		else
		{
			s << " : ";
			for(int m=0; m<dv(1); m++)
			{
				s << " " << m+1;
			}
			s << " :=\n";
			for(int m=0; m<dv(0); m++)
			{
				s << m+1;
				for(int n=0; n<dv(1); n++)
				{
					s << " " << P.contents(k).matrix_value().elem(m,n);
				}
				if(m < dv(0)-1)
				{
					s << "\n";
				}
			}
		}
		s << ";\n";
	}

	s << "\n" << "end;\n\n";
	myfile << s.str();
	myfile.close();
	return retval;
}
