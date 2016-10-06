#include "qp_abstract.hpp"


matrix QP_problem::solve()
{
	matrix solution;
	Solution s = CGAL::solve_quadratic_program(qp, ET());
	get_solution(s,solution);
	return solution;
};

void QP_problem::get_solution(Solution &s, matrix &solution){
	std::vector<double> v;
	assign_vector(s, v);
	vector_to_matrix(v, solution);
}



void QP_problem::get_solution(Solution &s, std::vector<double> &solution){
	assign_vector(s,solution);
}

void QP_problem::vector_to_matrix(std::vector<double> &v, matrix &m){
	m.set_size(v.size(),1);
	for (int i = 0; i < v.size(); i++)
		m(i, 0) = v[i];
}


void QP_problem::assign_A(matrix &A){
	for (int i = 0; i < A.nr(); i++){
		for (int j = 0; j < A.nc(); j++){
			qp.set_a(i, j, A(i, j));
		}
	}
}
void QP_problem::assign_A(std::vector<matrix> &A){
	for (int i = 0; i < A.size(); i++){
		for (int j = 0; j < A[i].nr(); j++){
			qp.set_a(i, j, A[i](j,0));
		}
	}
}
void QP_problem::assign_b(matrix &b){
	for (int i = 0; i < b.nr(); i++){
		qp.set_b(i, b(i, 0));
	}
}
void QP_problem::assign_D(matrix &D){
	for (int i = 0; i < D.nr(); i++){
		for (int j = 0; j < D.nc(); j++){
			qp.set_d(i, j, D(i, j));
		}
	}
}
void QP_problem::assign_c(matrix &c){
	for (int i = 0; i < c.nr(); i++){
		qp.set_c(i, c(i, 0));
	}
}