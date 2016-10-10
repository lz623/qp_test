#pragma once
#include <iostream>
#include <fstream>
#include "optimize.hpp"
#include <algorithm> 
#include <cmath>
#include <limits>
#include "dlib/optimization/optimization_abstract.h"
#include "dlib/optimization/optimization_search_strategies.h"
#include "dlib/optimization/optimization_stop_strategies.h"
#include "dlib/optimization/optimization_line_search.h"

class constrain_fun;
template <typename simulator>
class AL
{
	typedef dlib::matrix<double> matrix;
	typedef double(constrain_fun::*funprt)(const dlib::matrix<double> &x);
	typedef const dlib::matrix<double> (constrain_fun::*gradprt)(const dlib::matrix<double> &x);
private:
	matrix  lamda, sign_val; //constrain information
	simulator model;
	optimization opt;
/////////////////////////////////////////////constant parameters//////////////////////////////////
	const double omega_f_final;
	const double omega_x_final;
	const double eta_final;
	const double alpha;
	const double beta;
	const double w;
	const unsigned num_CS,num_Var;
	bool output_flag;
/////////////////////////////////////////////modified parameters//////////////////////////////////
	double miu, eta0, alpha_bar = 0.1;
	double eta, omega_f, omega_x, delta_cv;
	matrix x;
	//fun cs_fun = &constrain_fun::constrain_funct;

	double compute_constrain(const matrix &x);

	void initial_parameter();
	double compute_validation();
	void update_larange_lamda();
	void update_larange_miu();
	void load_matrix(std::string filename, matrix &m);



	double sign_function(const double &lamda, const double &miu, const double &val){
		double res(0);
		if (val<0) //constraint is voliated
		{
			res = (1 / (2 * miu))*val*val - miu*val;
		}
		else//constraint is satisfied
		{
			res = -0.5*lamda*lamda*miu;
		}
		return res;
	}

public:
	AL(matrix &startpoint, unsigned Ncs, bool output_i = false,
		double delta_f=1e-7,double delta_x=1e-5,double eta_final=1e-3) 
		:x(startpoint),num_CS(Ncs), num_Var(startpoint.nr()), output_flag(output_i),
		omega_f_final(delta_f), omega_x_final(delta_x), eta_final(eta_final), alpha(0.1), model(startpoint.nr(), Ncs)
		, beta(0.9), w(0.5){	
		load_matrix("input/constrain/lamda.dat", lamda);
		initial_parameter();
	}
	dlib::matrix<double> larangian_funct();
	double constrain_funct(const matrix &x);
	dlib::matrix<double> constrain_grad(const matrix &x);
};


template<typename simulator>
 double AL<simulator>::constrain_funct(const matrix &x)
{
	double result;
	result = model.run(x) + compute_constrain(x);
	return result;
}






template<typename simulator>
void AL<simulator>::initial_parameter()
{
	////////////////////////////////////////////should modify /////////////////////////////////////////!!!!!!!!!!!!!!!!!!!
	lamda.set_size(num_CS, 1); lamda = 0;
	//initial miu;
	sign_val.set_size(num_CS, 1);
	compute_constrain(x);
	double cx = dlib::sum(model.val);
	miu = 10 * cx / model.run(x);
	//compute initial miu;
	omega_f = std::min(miu, 0.1);
	omega_x = std::min(miu, 0.1);
	eta0 = pow(std::min(miu, 0.1), 0.1);
	eta = eta0;
}

template<typename simulator>
void AL<simulator>::update_larange_lamda()
{
	for (unsigned i = 0; i < num_CS; i++){
		lamda(i, 0) = std::max((lamda(i, 0) - model.val(i, 0) / miu), 0.0);
	}
	omega_f = std::max(omega_f*std::min(pow(miu, beta), w), omega_f_final);
	omega_x = std::max(omega_x*std::min(pow(miu, beta), w), omega_x_final);
	eta = std::max(eta*std::min(pow(miu, beta), w), eta_final);; //correct update
}
template<typename simulator>
void AL<simulator>::update_larange_miu()
{
	miu = 0.1*miu;
	//???????????????????????????????????????????????????? need to be fixed
	omega_f = std::max(omega_f*std::min(pow(miu, alpha), w), omega_f_final);
	omega_x = std::max(omega_x*std::min(pow(miu, alpha), w), omega_x_final);
	//????????????????????????????????????????????????????
	eta = std::max(eta*std::min(pow(miu, alpha), w), eta_final);
}

template<typename simulator>
double AL<simulator>::compute_validation()
{
	double nv(0), result(0); //nv here is number of volaties
	for (unsigned i = 0; i<num_CS; i++){
		if (model.val(i)<0){
			nv++;
			result += model.val(i)*model.val(i);
		}
	}
	if (nv != 0){
		result /= nv;
		return pow(result, 0.5);
	}
	else{
		return 0;
	}
}

template<typename simulator>
dlib::matrix<double>  AL<simulator>::larangian_funct()
{
	while (1){
		opt.find_min_new(dlib::bfgs_search_strategy(),
			dlib::objective_delta_stop_strategy(omega_f),
			dlib::x_delta_stop_strategy(omega_x),
			*this, x, output_flag);
		delta_cv = compute_validation();
		if (delta_cv < eta){
			if (eta = eta_final
				&& omega_f_final == omega_f
				&& omega_x_final == omega_x){
				return x; //AL is converge
			}
			else{
				update_larange_lamda();
				//update lamda with small violation
			}
		}
		else{
			update_larange_miu();
			//update miu with large violation
		}
	}
}


///////////////////////////////////////////////////////////should delete real version//////////////////////////////////////////////////////////////


//should be modify the cs gradient should computed from model
template<typename simulator>
dlib::matrix<double> AL<simulator>::constrain_grad(const matrix &x)
{
	matrix grad(x.nr(), 1);
	std::vector<dlib::matrix<double>> CS_gradient;
	grad = 0;
	compute_constrain(x);
	grad = model.grad(x);
	for (unsigned i = 0; i < num_CS; i++){
		if (model.val(i)<0){
			grad += (model.val(i) / miu - lamda(i))*model.CS_gradient[i];
		}
	}
	return grad;
}

template<typename simulator>
double AL<simulator>::compute_constrain(const matrix &x)
{
	double sum = 0;
	model.compute_constrain(x);
	for (unsigned i = 0; i < num_CS; i++){
		sign_val(i) = sign_function(lamda(i), miu, model.val(i));
	}

	return dlib::sum(sign_val);
}


 template<typename simulator>
 void AL<simulator>::load_matrix(std::string filename, matrix &m)
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

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////