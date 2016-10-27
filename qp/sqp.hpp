#pragma once
#include "qp_abstract.hpp"
#include "dlib/optimization/optimization_search_strategies.h"
#include "dlib/optimization/optimization_stop_strategies.h"
#include "optimize.hpp"
#include <algorithm>
#include <stdexcept>
template<typename simulator>
class SQP
{
	typedef dlib::matrix<double> matrix;
private:
	unsigned num_Var, num_CS;
	matrix lamda,qp_lamda,miu ,x;
	simulator model;
	QP_problem qp_solver;
	dlib::bfgs_search_strategy  BFGS;
	dlib::x_delta_stop_strategy stop_criterion;
	std::vector<int> active_index, inactive_index, mactive_index, minactive_index;
	matrix delta,gamma;
	double alpha,min_alpha;
	double f_value,f_original,z_value,z_grad,z_step;
	double nit;
	double last_alpha;
	bool been_used,output_res,powell_BFGS;
	matrix grad, g,Y,g_old; //grad is gradient for constrain function, g is gradient for model
	optimization opt;
	std::vector<dlib::matrix<double>> CSo_gradient;
	std::vector<int> fi;
	void initial_parameters();
	void update_lamda(matrix &lamda);
	void update_lamda(matrix &lamda, std::vector<int> &active_index);
	double compute_constrain();
	void compute_delta();
	void comput_alpha();
	bool active_cs(double x1, double x2, double toler = 1e-5);
	bool active_cs_1(double lamda, double miu, double r);
	matrix extract_active_set(matrix &delta);
	matrix extract_active_set(std::vector<int> &active_index);
	void extract_all_set(const matrix &delta, matrix &A_a
		,matrix &A_i,std::vector<int> &active_index, std::vector<int> &inactive_index);
	void merit();
	double merit_value();
	double merit_grad();
	void compute_penelty();
	void get_Hessian(const matrix& delta, matrix& gamma);
	void initial_Hessian();
public:
	SQP(matrix &startpoint, unsigned ncs, double conc = 1e-5, bool output_flag = false, bool powell_flag = true)
		:num_Var(startpoint.nr()), num_CS(ncs), x(startpoint), qp_solver(startpoint.nr(), ncs),
		model(startpoint.nr(), ncs), stop_criterion(conc), output_res(output_flag), powell_BFGS(powell_flag)
	{
		initial_parameters();
	}

	matrix solve_SQP();
	matrix solve_SQP_new();
	double constrain_funct(const matrix &x);
	double constrain_funct(double fvalue);
	matrix constrain_grad(const matrix &x,double &fvalue);
};

template<typename simulator>
void SQP<simulator>::initial_parameters(){
	alpha = 1;
	min_alpha = 1e-6;
	last_alpha = 1;
	lamda.set_size(num_CS, 1);
	qp_lamda.set_size(num_CS, 1);
	miu.set_size(num_CS, 1);
	miu = 1e-2;
	lamda = 0;
	qp_lamda = 0;
	delta.set_size(num_Var, 1);
	gamma.set_size(num_Var, 1);
	grad.set_size(num_Var, 1), g.set_size(num_Var, 1);
	Y.set_size(num_Var, num_Var);
	CSo_gradient.resize(num_CS);
	for (unsigned i = 0; i < num_CS; i++)
		CSo_gradient[i].set_size(num_Var, 1);
	g_old.set_size(num_Var, 1);
	fi.resize(num_CS);
	for (unsigned i = 0; i < num_CS; i++) fi[i] = i;
}

template<typename simulator>
matrix SQP<simulator>::constrain_grad(const matrix &x,double &f_value)
{
	g = model.grad(x);
	f_original = model.f_value;
	f_value = f_original + constrain_funct(x);
	gamma = g - g_old;
	grad = g;
	std::cout << "qp_lamda: " << std::endl << qp_lamda << std::endl;
	for (unsigned i = 0; i < num_CS; i++){
		if (model.val(i) < 0){
			grad += -qp_lamda(i)*model.CS_gradient[i];
			gamma += -qp_lamda(i)*(model.CS_gradient[i] - CSo_gradient[i]);
		}
		CSo_gradient[i] = model.CS_gradient[i];
	}
	g_old = g;
	return grad;
}


template<typename simulator>
double SQP<simulator>::constrain_funct(const matrix &x)
{
	double result;
	f_original = model.run(x);
	model.compute_constrain(x);
	result = f_original-dot(lamda, model.val);
	return result;
}
template<typename simulator>
double SQP<simulator>::constrain_funct(double fvalue)
{
	double result(fvalue);
	model.compute_constrain(x);
	result -= dot(lamda, model.val);
	return result;
}


template<typename simulator>
void SQP<simulator>::update_lamda(matrix &lamda)
{
	//use x_up1 to compute active set this time
	matrix A_a=	extract_active_set(delta);
	lamda = 0;
	if (active_index.size()!=0)
	{	matrix lamda_act = inv(A_a*trans(A_a))*A_a*(Y*delta + g);
		for (unsigned i = 0; i < active_index.size(); i++){
			lamda(active_index[i], 0) = abs(lamda_act(i, 0));
		}
	}

}

template<typename simulator>
void SQP<simulator>::update_lamda(matrix &lamda,std::vector<int> &active_index)
{
	//use x_up1 to compute active set this time
	matrix A_a = extract_active_set(active_index);
	lamda = 0;
	if (active_index.size() != 0)
	{
		matrix lamda_act = inv(A_a*trans(A_a))*A_a*(Y*delta + g);
		for (unsigned i = 0; i < active_index.size(); i++){
			lamda(active_index[i], 0) = abs(lamda_act(i, 0));
		}
	}

}

template<typename simulator>
matrix SQP<simulator>::extract_active_set(matrix &delta){
	std::vector<matrix> active_A;  //Don't know why no response;
	double c;
	active_index.clear();
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
matrix SQP<simulator>::extract_active_set(std::vector<int> &active_index){
	std::vector<matrix> active_A;  //Don't know why no response;
	int index;
	matrix A_a(active_index.size(), num_Var);
	for (unsigned i = 0; i < active_index.size();i++)
	{
		index = active_index[i];
		set_rowm(A_a, i) = trans(CSo_gradient[i]);
	}
	return A_a;
}

template<typename simulator>
bool SQP<simulator>::active_cs(double x1, double x2, double toler = 1e-5){
	return ((x1+x2)<toler);
}
template<typename simulator>
bool SQP<simulator>::active_cs_1(double lamda,double miu,double c){
	return (lamda/miu>=c);
}

template<typename simulator>
void SQP<simulator>::comput_alpha(){
	//Use the method proposed in Practical Optimization Andreas p510
	double rho(1e-4),a_2;
	alpha = 1;
	matrix x_n = x + alpha*delta;
	matrix lamda_old = lamda;
	lamda = qp_lamda;
	constrain_funct(x_n);
	z_step=merit_value();
	int niit = 0;
	while (z_step>z_value+rho*alpha*z_grad
		&& alpha>min_alpha)
	{
		a_2 = ((z_step - z_value) / alpha - z_grad) / alpha;
		alpha = std::max(0.1*alpha, -z_grad/(2*a_2));
		x_n = x + alpha*delta;
		constrain_funct(x_n);
		z_step = merit_value();
		lamda = lamda_old + alpha*(qp_lamda - lamda_old);
	} 

}

//template<typename simulator>
//matrix SQP<simulator>::solve_SQP(){
//	unsigned i=0;
//	nit=0;
//	while (stop_criterion.should_continue_search(x))
//	{
//		nit++;
//		compute_delta();
//		comput_alpha();
//		delta = alpha*delta;
//		x = x + delta;
//		std::cout << i++ << std::endl << z_value << std::endl;
//	}
//	return x;
//}


template<typename simulator>
matrix SQP<simulator>::solve_SQP(){
	unsigned i = 0;
	nit = 1;
	initial_Hessian();
	grad = constrain_grad(x, f_value);
	while (stop_criterion.should_continue_search(x))
	{
		delta = qp_solver.solve_qp(model.CS_gradient, -model.val, Y, g);//solve to get delta;
		update_lamda(qp_lamda,fi);
		merit();
		comput_alpha();
		delta = alpha*delta;
		x = x + delta;
		grad = constrain_grad(x, f_value);
		//std::cout << "cs:" << std::endl << model.val << std::endl;
		//std::cout << "delta:" << std::endl << delta << std::endl;
		//std::cout << "gamma:" << std::endl << gamma << std::endl;
		get_Hessian(delta, gamma);
		//std::cout << i++ << std::endl << z_value << std::endl;
	}
	return x;
}

template<typename simulator>
void SQP<simulator>::compute_delta()
{
	grad = constrain_grad(x, f_value);
	get_Hessian(delta, gamma);
	delta = qp_solver.solve_qp(model.CS_gradient, -model.val, Y, g);
	update_lamda(qp_lamda);
	merit();
}

template<typename simulator>
void SQP<simulator>::extract_all_set(const matrix &delta, matrix &A_a
	, matrix &A_i, std::vector<int> &active_index, std::vector<int> &inactive_index){
	double c;
	active_index.clear();
	inactive_index.clear();
	for (unsigned i = 0; i < num_CS; i++){
		matrix cv = model.CS_gradient[i] * trans(delta);
		c = cv(0, 0);
		if (active_cs(lamda(i),miu(i),model.val(i))){
			active_index.push_back(i);
		}
		else{
			inactive_index.push_back(i);
		}
	}
	//N_acs*Nvar matrix
	A_a.set_size(active_index.size(), num_Var);
	A_i.set_size(inactive_index.size(), num_Var);
	for (int i = 0; i < active_index.size(); i++){
		set_rowm(A_a, i) = trans(model.CS_gradient[active_index[i]]);
	}
	for (int i = 0; i < inactive_index.size(); i++){
		set_rowm(A_i, i) = trans(model.CS_gradient[inactive_index[i]]);
	}
}


template<typename simulator>
void SQP<simulator>::merit()
{
	matrix Ai, Aa;
	compute_penelty();
	z_value=merit_value();
	z_grad=merit_grad();
	if (z_grad > 0)
		throw std::invalid_argument("Bad search direction from QP");
}



template<typename simulator>
double SQP<simulator>::merit_grad()
{
	double tmp(0),tmp1(0);
	matrix grad2 = g;
	int a, ia;
	for (unsigned i = 0; i < mactive_index.size(); i++){
		a = mactive_index[i];
		grad2 += (lamda(a)-miu(a)*model.val(a))*(-model.CS_gradient[a]);
		tmp = model.val(a)*(qp_lamda(a)-lamda(a));
	}
	for (unsigned i = 0; i <minactive_index.size(); i++){
		ia = minactive_index[i];
		tmp1 = lamda(ia) / miu(ia)*(qp_lamda(ia) - lamda(ia));
	}
	return trans(grad2)*delta - tmp1 - tmp;
};

template<typename simulator>
double SQP<simulator>::merit_value()
{
	int a, ia;
	double res = f_original;
	mactive_index.clear();
	minactive_index.clear();
	for (unsigned i = 0; i <num_CS; i++)
	{
		if (active_cs_1(lamda(i), miu(i), model.val(i))){
			mactive_index.push_back(i);
		}
		else{
			minactive_index.push_back(i);
		}
	}

	for (unsigned i = 0; i <mactive_index.size() ;i++){
		a = mactive_index[i];
		double act_res1 = lamda(a)*model.val(a);
		double act_res=- 0.5*miu(a)*model.val(a)*model.val(a);
		res -= (act_res + act_res1);
	}
	for (unsigned i = 0; i <minactive_index.size(); i++){
		ia = minactive_index[i];
		res -= 0.5*lamda(ia)*lamda(ia) / miu(ia);
	}
	return res;
};

template<typename simulator>
void SQP<simulator>::compute_penelty()
{
	matrix sigma(num_CS,1);
	matrix sHs = trans(delta)*Y*delta + 1;
	double dHd = sHs(0, 0);
	for (unsigned i = 0; i < num_CS;i++)
	{
		sigma(i) = std::min(miu(i),nit*pow(miu(i),0.5));
		miu(i) = std::min(1e9, 2 * num_CS*
			(qp_lamda(i) - lamda(i))*(qp_lamda(i) - lamda(i)) / dHd);
		miu(i) = std::max(sigma(i), miu(i));
	}
}

template<typename simulator>
void SQP<simulator>::get_Hessian(
	const matrix& delta,
	 matrix&gamma
	)
{
	double theta = 0;
	if (been_used == false)
	{
		been_used = true;
		Y = dlib::identity_matrix<double>(delta.size());
	}
	else
	{
		double dg = dot(delta, gamma);
		matrix	dHd_m = trans(delta)*Y*delta;
		matrix  Hd = Y*delta;
		double dHd = dHd_m(0, 0);
		if (powell_BFGS)
		{
			if (dg > 0.2*dHd)
				theta = 1;
			else
				theta = (0.8*dHd) / (dHd - dg);
			gamma = theta*gamma + (1 - theta)*Hd;
		}
		std::cout << "delta " << std::endl << delta << std::endl;
		std::cout << "gamma: " << std::endl << gamma << std::endl;
		//std::cout << "Hd: " << std::endl << Hd << std::endl;
		//std::cout << "dg " << std::endl << dg << std::endl;
		if (dot(delta, gamma)>0)
			Y = Y + gamma*trans(gamma) / dg - (Hd*trans(Hd)) / dHd;

	}
}

template<typename simulator>
void SQP<simulator>::initial_Hessian()
{
	double theta = 0;
	been_used = true;
	Y = dlib::identity_matrix<double>(num_Var);
}