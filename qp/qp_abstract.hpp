#pragma once
#include <cassert>
#include <CGAL/basic.h>
#include <CGAL/QP_models.h>
#include <CGAL/QP_functions.h>
#include <dlib/matrix.h>



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
typedef dlib::matrix<double> matrix;
class QP_problem
{
public:
	QP_problem(unsigned nv,unsigned ncs) :num_var(nv),
		num_cs(ncs),qp(CGAL::LARGER, true, 0, false, 0)
	{}
	matrix solve();
	void assign_A(const matrix &A);
	void assign_A(const std::vector<matrix> &A);
	void assign_b(const matrix &b);
	void assign_D(const matrix &D);
	void assign_c(const matrix &c);
	matrix solve_qp(const std::vector<matrix> &A,
		const matrix &b, const matrix &D, const matrix& c);
private:
	void get_solution(Solution &s, std::vector<double> &solution);
	void get_solution(Solution &s, matrix &solution);
	void vector_to_matrix(std::vector<double> &v, matrix &m);
	Program qp;
	unsigned num_var, num_cs;
};

