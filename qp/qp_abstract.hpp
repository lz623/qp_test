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
		num_cs(ncs),qp(CGAL::SMALLER, true, 0, false, 0)
	{}
	matrix solve();
	void assign_A(matrix &A);
	void assign_A(std::vector<matrix> &A);
	void assign_b(matrix &b);
	void assign_D(matrix &D);
	void assign_c(matrix &c);
private:
	void get_solution(Solution &s, std::vector<double> &solution);
	void get_solution(Solution &s, matrix &solution);
	void vector_to_matrix(std::vector<double> &v, matrix &m);
	Program qp;
	unsigned num_var, num_cs;
};

