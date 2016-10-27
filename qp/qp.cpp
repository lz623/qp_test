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
		//load_matrix("input/constrain/constrain.dat", CS);
		//load_matrix("input/constrain/constrain_val.dat", CS_val);
		//load_matrix("input/constrain/constrain_range.dat", CS_range);
		//val.set_size(num_Cs, 1);
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
		f_value=run(m);
		dlib::matrix<double> res(2, 1);
		res(0, 0) = -400 * x*(y - x*x) - 2 * (1 - x); // derivative of rosen() with respect to x
		res(1, 0) = 200 * (y - x*x);              // derivative of rosen() with respect to y
		compute_constrain_grad(m);
		num_gcall++;
		return res;
	}

	void compute_constrain_grad(const dlib::matrix<double>& m)
	{
		CS_gradient[0](0) = -5-2*m(0);
		CS_gradient[0](1) = 6-2*m(1);
		CS_gradient[1](0) = -5 -2 * m(0);
		CS_gradient[1](1) = -6-2*m(1);
	}

	void load_cs_gradient()
	{
		//std::ifstream input("input/constrain/cs_gradient.dat");
		num_Cs = 2;
		num_Var = 2;
		val.set_size(num_Cs,1);
		CS_gradient.resize(num_Cs);
		for (unsigned i = 0; i <CS_gradient.size(); i++){
			CS_gradient[i].set_size(num_Var, 1);
		}
//only for this case
	}
	// linear constrain
	//void compute_constrain(dlib::matrix<double> m)
	//{
	//	//should include simulation run in reservoir simulator
	//	for (int i = 0; i < num_Cs; i++){
	//		val(i) = (rowm(CS, i)*m + CS_val(i)) / abs(CS_range(i));
	//	}
	//}
	// nolinear constrain
	void compute_constrain(dlib::matrix<double> m)
	{
		val(0) = 16 - pow((m(0)+2.5),2) - pow((m(1) - 3),2);
		val(1) = 16 - pow((m(0)+2.5),2) - pow((m(1) + 3), 2);
	}
	

	dlib::matrix<double> val;
	std::vector<dlib::matrix<double>> CS_gradient;
	double f_value;
};



int main() {
	dlib::matrix<double> starting_point(2, 1);
	starting_point(0, 0) = -3; // Start with a valid point inside the constraint box.
	starting_point(1, 0) = 0.8;
	AL<rosen_model> test(starting_point, 2, 1);
	dlib::matrix<double> result = test.larangian_funct();
	std::cout << result << std::endl;
	//// by default, we have a nonnegative QP with Ax <= b
	//dlib::matrix<int> i(1,1);
	//i = 0;
	//std::cout << (result,0);
	starting_point(0, 0) = -3; // Start with a valid point inside the constraint box.
	starting_point(1, 0) = 0.8;
	SQP<rosen_model> sqp_solver(starting_point,2);
	std::cout << sqp_solver.solve_SQP_new();
	// solve the program, using ET as the exact type

	// output solution

	return 0;
}
