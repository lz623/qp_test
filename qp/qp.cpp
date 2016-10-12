// example: construct a quadratic program from data
// the QP below is the first quadratic program example in the user manual
#include <iostream>
#include <cassert>
#include <CGAL/basic.h>
#include <CGAL/QP_models.h>
#include <CGAL/QP_functions.h>
#include <algorithm>
#include "constrain_opt_AL.hpp"
#include "sqp.hpp"
#include <dlib/optimization.h>


// choose exact integral type
#ifdef CGAL_USE_GMP
#include <CGAL/Gmpz.h>
typedef CGAL::Gmpzf ET;
#else
#include <CGAL/MP_Float.h>
typedef CGAL::MP_Float ET;
#endif

// program and solution types
typedef CGAL::Quadratic_program<double> Program;
typedef CGAL::Quadratic_program_solution<ET> Solution;

void load_matrix(std::string filename, dlib::matrix<double> &m)
{
	std::ifstream input(filename);
	unsigned nc, nr;
	input >> nr;
	input >> nc;
	m.set_size(nr, nc);
	m = 0;
	for (unsigned j = 0; j < nr; j++){
		for (unsigned i = 0; i < nc; i++){
			input >> m(j, i);
		}
	}
};


class rosen_model
{

private:

	int num_Cs, num_Var;
	dlib::matrix<double> CS, CS_val, CS_range;
public:
	 unsigned num_fcall, num_gcall;
	rosen_model(unsigned nv, unsigned ncs) :num_Var(nv), num_Cs(ncs){
		num_fcall = 0;
		num_gcall = 0;
		load_cs_gradient();
		load_matrix("input/constrain/constrain.dat", CS);
		load_matrix("input/constrain/constrain_val.dat", CS_val);
		load_matrix("input/constrain/constrain_range.dat", CS_range);
		val.set_size(num_Cs, 1);
	};

	double run(dlib::matrix<double> m){
		const double x = m(0, 0);
		const double y = m(1, 0);
		// compute Rosenbrock's function and return the result
		num_fcall++;
		return 100.0*pow(y - x*x, 2) + pow(1 - x, 2);

	};

	const dlib::matrix<double> grad(const dlib::matrix<double>& m){
		const double x = m(0, 0);
		const double y = m(1, 0);
		dlib::matrix<double> res(2, 1);
		res(0, 0) = -400 * x*(y - x*x) - 2 * (1 - x); // derivative of rosen() with respect to x
		res(1, 0) = 200 * (y - x*x);              // derivative of rosen() with respect to y
		num_gcall++;
		return res;
	}

	void load_cs_gradient()
	{
		std::ifstream input("input/constrain/cs_gradient.dat");
		CS_gradient.resize(num_Cs);
		for (unsigned i = 0; i <CS_gradient.size(); i++){
			CS_gradient[i].set_size(num_Var, 1);
			for (unsigned j = 0; j < num_Var; j++){
				input >> CS_gradient[i](j, 0);
			}
		}
	}
	void compute_constrain(dlib::matrix<double> m)
	{
		//should include simulation run in reservoir simulator
		for (unsigned i = 0; i < num_Cs; i++){
			val(i) = (rowm(CS, i)*m + CS_val(i)) / abs(CS_range(i));
		}
	}
	dlib::matrix<double> val;
	std::vector<dlib::matrix<double>> CS_gradient;
};




int main() {
	dlib::matrix<double> starting_point(2, 1);
	starting_point(0, 0) = -1; // Start with a valid point inside the constraint box.
	starting_point(1, 0) = 2.5;
	AL<rosen_model> test(starting_point, 4, 1);
	dlib::matrix<double> result = test.larangian_funct();
	std::cout << result << std::endl;
	// by default, we have a nonnegative QP with Ax <= b

	SQP<rosen_model> sqp_solver(starting_point,4);
	std::cout << sqp_solver.solve_SQP();

	// solve the program, using ET as the exact type

	// output solution

	return 0;
}
