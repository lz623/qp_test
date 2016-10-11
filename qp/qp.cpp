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
	unsigned num_Cs, num_Var;
	dlib::matrix<double> CS, CS_val, CS_range;
public:

	rosen_model(unsigned nv, unsigned ncs) :num_Var(nv), num_Cs(ncs){
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
		return 100.0*pow(y - x*x, 2) + pow(1 - x, 2);
	};

	const dlib::matrix<double> grad(const dlib::matrix<double>& m){
		const double x = m(0, 0);
		const double y = m(1, 0);
		dlib::matrix<double> res(2, 1);
		res(0, 0) = -400 * x*(y - x*x) - 2 * (1 - x); // derivative of rosen() with respect to x
		res(1, 0) = 200 * (y - x*x);              // derivative of rosen() with respect to y
		return res;
	}

	void load_cs_gradient()
	{
		std::ifstream input("input/constrain/cs_gradient.dat");
		CS_gradient.resize(num_Cs);
		for (unsigned i = 0; i <CS_gradient.size(); i++)
		{
			CS_gradient[i].set_size(num_Var, 1);
			for (unsigned j = 0; j < num_Var; j++)
			{
				input >> CS_gradient[i](j, 0);
			}
		}
	}
	void compute_constrain(dlib::matrix<double> m)
	{
		//should include simulation run in reservoir simulator
		for (unsigned i = 0; i < num_Cs; i++)
		{
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
	Program qp(CGAL::LARGER, true, 0, false, 0);
	// now set the non-default entries: 
	const double X = 0.1;
	const double Y = 1.5;
	qp.set_a(0, 0, X); qp.set_a(1,0, Y); qp.set_b(0, 7.5);  //  x + y  <= 7
	qp.set_a(0, 1, -1); qp.set_a(1, 1, 2); qp.set_b(1, 4);  // -x + 2y <= 4
	//qp.set_u(Y, true, 4);                                   //       y <= 4
	qp.set_d(0, 0, 2); qp.set_d(1, 1, 8); // !!specify 2D!!    x^2 + 4 y^2
	qp.set_c(1, -32);                                       // -32y
	//qp.set_c0(64);
	// +64
	Solution s = CGAL::solve_quadratic_program(qp, ET());
	assert(s.solves_quadratic_program(qp));
	std::vector<double> x;
	assign_vector(s, x);
	std::cout << s << std::endl;

	QP_problem qp_solver(2,2);
	dlib::matrix<double> A(2, 2), b(2, 1), c(2, 1), D(2, 2);
	A = 1, -1,
		1, 2;
	b = 7,
		4;
	c = 0,
		-32;
	D = 2, 0,
		0, 8;
	std::cout << qp_solver.solve_qp(A, b, D,c) << std::endl;
	// solve the program, using ET as the exact type

	// output solution

	return 0;
}
