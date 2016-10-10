#include "qp_abstract.hpp"
#include "dlib/optimization/optimization_search_strategies.h"

template<typename simulator>
class SQP
{
	typedef dlib::matrix<double> matrix;
private:
	unsigned num_Var, num_CS;
	matrix lamda,x;
	simulator model;
	QP_problem qp_solver;
	dlib::bfgs_search_strategy  BFGS;
	matrix delta;
	double alpha;
	const double coverge_c;
	matrix grad(num_Var, 1);

	void initial_parameters();
	void update_lamda();
	double constrain_funct(const matrix &x);
	matrix constrain_grad(const matrix &x);
	double compute_constrain();
	void compute_delta();
	void comput_alpha();

	matrix extract_active_set();
public:
	SQP(matrix &startpoint, unsigned ncs,double conc=1e-5) :num_Var(startpoint.nr()), num_Cs(ncs)
		, x(startpoint), qp_solver(startpoint.nr(), ncs), coverge_c(conc){
		initial_parameters();
	}
	matrix solve_SQP();
};



template<typename simulator>
matrix SQP<simulator>::constrain_grad(const matrix &x)
{
	std::vector<dlib::matrix<double>> CS_gradient;
	grad = 0;
	compute_constrain(x);
	grad = model.grad(x);
	for (unsigned i = 0; i < num_CS; i++){
		if (model.val(i)<0){
			grad += lamda(i)*model.CS_gradient[i];
		}
	}
	return grad;
}


template<typename simulator>
double SQP<simulator>::constrain_funct(const matrix &x)
{
	double result;
	result = model.run(x) + dot(lamda,model.val);
	return result;
}

template<typename simulator>
void SQP<simulator>::compute_delta()
{
	double f_value = constrain_funct(x);
	matrix Y = BFGS.get_next_direction(x, model.f_value, constrain_grad(x));
	delta=qp_solver.solve_qp(model.CS_gradient,-model.val,Y,grad);
	update_lamda();
}

template<typename simulator>
void SQP<simulator>::update_lamda()
{
	matrix A_a=extract_active_set();
	//check if the matrix computation valid £¡£¡£¡£¡£¡£¡£¡£¡£¡£¡£¡£¡£¡£¡£¡£¡£¡£¡
	lamda = inv(A_a*trans(A_a))*A_a*(Y*delta+grad);
}


template<typename simulator>
matrix SQP<simulator>:: extract_active_set(){

}

template<typename simulator>
void SQP<simulator>::comput_alpha(){
	//Use the method proposed in Practical Optimization Andreas p510

}