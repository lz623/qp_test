#pragma once
#include "qp_abstract.hpp"
#include "dlib/optimization/optimization_search_strategies.h"
#include "dlib/optimization/optimization_stop_strategies.h"
#include "optimize.hpp"
template<typename simulator>
class SQP
{
	typedef dlib::matrix<double> matrix;
private:
	unsigned num_Var, num_CS;
	matrix lamda, x, x_up1;
	simulator model;
	QP_problem qp_solver;
	dlib::bfgs_search_strategy  BFGS;
	dlib::x_delta_stop_strategy stop_criterion;
	std::vector<unsigned> active_index;
	matrix delta;
	double alpha;
	double f_value;
	double last_alpha;
	matrix grad, g,Y; //grad is gradient for constrain function, g is gradient for model
	optimization opt;

	void initial_parameters();
	void update_lamda();
	double compute_constrain();
	void compute_delta();
	void comput_alpha();
	bool active_cs(double x1, double x2, double toler = 1e-5);
	matrix 	extract_active_set(matrix &delta);
public:
	SQP(matrix &startpoint, unsigned ncs,double conc=1e-5) :num_Var(startpoint.nr()), num_CS(ncs)
		, x(startpoint), qp_solver(startpoint.nr(), ncs), model(startpoint.nr(), ncs), stop_criterion(conc){
		initial_parameters();
	}

	matrix solve_SQP();
	double constrain_funct(const matrix &x);
	matrix constrain_grad(const matrix &x);
};

template<typename simulator>
void SQP<simulator>::initial_parameters(){
	alpha = 0.1;
	last_alpha = 0.1;
	lamda.set_size(num_CS, 1);
	lamda = 0;
	x_up1.set_size(num_Var, 1);
	delta.set_size(num_Var, 1);
	grad.set_size(num_Var, 1), g.set_size(num_Var, 1);
	Y.set_size(num_Var, num_Var);
}

template<typename simulator>
matrix SQP<simulator>::constrain_grad(const matrix &x)
{
	std::vector<dlib::matrix<double>> CS_gradient;
	g = 0;
	g = model.grad(x);
	grad = g;
	for (unsigned i = 0; i < num_CS; i++){
		if (model.val(i)<0){
			grad +=-lamda(i)*model.CS_gradient[i];
		}
	}
	return grad;
}


template<typename simulator>
double SQP<simulator>::constrain_funct(const matrix &x)
{
	double result;
	result = model.run(x);
	model.compute_constrain(x);
	result+=dot(lamda, model.val);
	return result;
}

template<typename simulator>
void SQP<simulator>::compute_delta()
{
	f_value = constrain_funct(x);
	grad=constrain_grad(x);
	Y = BFGS.get_Hessian(x, f_value, grad);
	delta=qp_solver.solve_qp(model.CS_gradient,-model.val,Y,g);
}

template<typename simulator>
void SQP<simulator>::update_lamda()
{
	//use x_up1 to compute active set this time
	matrix A_a = extract_active_set(delta);
	lamda = 0;
	matrix lamda_act = inv(A_a*trans(A_a))*A_a*(Y*delta+g);
	for (unsigned i = 0; i < active_index.size();i++){
		lamda(active_index[i], 0) = lamda_act(i,0);
	}
	active_index.clear();
}


template<typename simulator>
matrix SQP<simulator>::extract_active_set(matrix &delta){
	std::vector<matrix> active_A;  //Don't know why no response;
	double c;
	for (unsigned i = 0; i < num_CS; i++){
		matrix cv = model.CS_gradient[i] * trans(delta);
		c = cv(0,0);
		if (active_cs(c,model.val(i))){
			active_A.push_back(model.CS_gradient[i]);
			active_index.push_back(i);
		}
	}
	//N_acs*Nvar matrix
	matrix A_a(active_A.size(),num_Var);
	for (int i = 0; i < active_A.size();i++){
		set_rowm(A_a, i) = trans(active_A[i]);
	}
	return A_a;
}

template<typename simulator>
bool SQP<simulator>::active_cs(double x1, double x2, double toler = 1e-5){
	return (abs(x1+x2)<toler);
}

template<typename simulator>
void SQP<simulator>::comput_alpha(){
	//Use the method proposed in Practical Optimization Andreas p510
	double alpha1;
	alpha1 = opt.backtracking_line_search_new(
		*this, x, delta,
		f_value,
		dot(grad, delta), // compute gradient for the line search
		last_alpha,
		BFGS.get_wolfe_rho(),
		BFGS.get_max_line_search_iterations());
	alpha = 0.95*alpha1;
	//alpha_ = 0.95*min(alpha1,alpha2);
}

template<typename simulator>
matrix SQP<simulator>::solve_SQP(){
	while (stop_criterion.should_continue_search(x))
	{
		compute_delta();
		update_lamda();
		comput_alpha();
		delta = alpha*delta;
		update_lamda();
		x = x + delta;
		last_alpha = alpha;
		std::cout << x << std::endl;
		std::cout << lamda << std::endl;
	}
	return x;
}