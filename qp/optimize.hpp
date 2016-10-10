#pragma once
#include <cmath>
#include <limits>
#include <fstream>
#include <iostream>
#include "dlib/optimization/optimization_abstract.h"
#include "dlib/optimization/optimization_search_strategies.h"
#include "dlib/optimization/optimization_stop_strategies.h"
#include "dlib/optimization/optimization_line_search.h"


class optimization
{
public:
	optimization()
	{};

	template <
		typename search_strategy_type,
		typename stop_strategy_type,
		typename stop_strategy_type_1,
		typename model,
		typename T
	>
	double find_min_new(
	search_strategy_type search_strategy,
	stop_strategy_type stop_strategy,
	stop_strategy_type_1 stop_strategy_1,
	model& model,
	T& x,
	bool & is_output
	)
	{
		//output1.open("output/obj.dat", ios::app);
		//output2.open("output/iterpoint.dat", ios::app);
		T g, s;
		double f_value = model.constrain_funct(x);
		g = model.constrain_grad(x);
		double last_alpha = 0.1;
		stop_strategy_1.set_used(x);
		while (stop_strategy.should_continue_search(x, f_value, g) 
			|| stop_strategy_1.should_continue_search(x, f_value, g))
		{
			s = search_strategy.get_next_direction(x, f_value, g);
			double alpha = backtracking_line_search_new(
				model,x, s,
				f_value,
				dot(g, s), // compute gradient for the line search
				last_alpha,
				search_strategy.get_wolfe_rho(),
				search_strategy.get_max_line_search_iterations());
			if (alpha == last_alpha)
				last_alpha = std::min(last_alpha * 10, 1.0);
			else
				last_alpha = alpha;
			// Take the search step indicated by the above line search
			x += alpha*s;
			g = model.constrain_grad(x);
			//if (is_output){
			//	std::cout << (stop_strategy.should_continue_search(x, f_value, g)
			//		|| stop_strategy_1.should_continue_search(x, f_value, g)) << std::endl;
			//	std::cout << "x: " << trans(x) << std::endl;
			//	std::cout << "g: " << trans(g) << std::endl;
			//	std::cout << "f_value: " << f_value << std::endl;
			//	std::cout << "-----------------------------------------------------------" << std::endl;
			//}
		}
		return f_value;
	}

	template <
		typename search_strategy_type,
		typename stop_strategy_type,
		typename funct,
		typename funct_der,
		typename T
	>
	double find_min(
	search_strategy_type search_strategy,
	stop_strategy_type stop_strategy,
	funct& f,
	funct_der& der,
	T& x
	)
	{
		T g, s;
		double f_value = f(x);
		g = der(x);
		double last_alpha = 1;
		while (stop_strategy.should_continue_search(x, f_value, g))
		{
			s = search_strategy.get_next_direction(x, f_value, g);

			double alpha = backtracking_line_search(
				make_line_search_function(f, x, s, f_value),
				f_value,
				dot(g, s), // compute gradient for the line search
				last_alpha,
				search_strategy.get_wolfe_rho(),
				search_strategy.get_max_line_search_iterations());

			if (alpha == last_alpha)
				last_alpha = std::min(last_alpha * 10, 1.0);
			else
				last_alpha = alpha;
			// Take the search step indicated by the above line search
			x += alpha*s;
			g = der(x);
		}
		return f_value;
	}





	template <typename model>
	double backtracking_line_search_new(
		model f,
		dlib::matrix<double> x,
		dlib::matrix<double> grad,
		double &f0,
		double d0,
		double alpha,
		double rho,
		unsigned long max_iter
		)
	{
		if ((d0 > 0 && alpha > 0) ||
			(d0 < 0 && alpha < 0))
		{
			alpha *= -1;
		}
		bool have_prev_alpha = false;
		double prev_alpha = 0;
		double prev_val = 0;
		unsigned long iter = 0;
		while (true)
		{
			++iter;
			const double val = f.constrain_funct(alpha*grad + x);
			if (val <= f0 + alpha*rho*d0 || iter >= max_iter)
			{
				f0 = val;
				return alpha;
			}
			else
			{
				// Interpolate a new alpha.  We also make sure the step by which we
				// reduce alpha is not super small.
				double step;
				if (!have_prev_alpha)
				{
					if (d0 < 0)
						step = alpha*put_in_range(0.1, 0.9, poly_min_extrap(f0, d0, val));
					else
						step = alpha*put_in_range(0.1, 0.9, poly_min_extrap(f0, -d0, val));
					have_prev_alpha = true;
				}
				else
				{
					if (d0 < 0)
						step = put_in_range(0.1*alpha, 0.9*alpha, poly_min_extrap(f0, d0, alpha, val, prev_alpha, prev_val));
					else
						step = put_in_range(0.1*alpha, 0.9*alpha, -poly_min_extrap(f0, -d0, -alpha, val, -prev_alpha, prev_val));
				}

				prev_alpha = alpha;
				prev_val = val;

				alpha = step;
			}
		}
	}

	// ----------------------------------------------------------------------------------------
	private:
	double put_in_range(const double& a,const double& b,const double& val)
	{
		if (a < b)
		{
			if (val < a)
				return a;
			else if (val > b)
				return b;
		}
		else
		{
			if (val < b)
				return b;
			else if (val > a)
				return a;
		}

		return val;
	}

	inline double poly_min_extrap(
		double f0,
		double d0,
		double x1,
		double f_x1,
		double x2,
		double f_x2
		)
	{
		DLIB_ASSERT(0 < x1 && x1 < x2, "Invalid inputs were given to this function");
		// The contents of this function follow the equations described on page 58 of the
		// book Numerical Optimization by Nocedal and Wright, second edition.
		dlib::matrix<double, 2, 2> m;
		dlib::matrix<double, 2, 1> v;

		const double aa2 = x2*x2;
		const double aa1 = x1*x1;
		m = aa2, -aa1,
			-aa2*x2, aa1*x1;
		v = f_x1 - f0 - d0*x1,
			f_x2 - f0 - d0*x2;


		double temp = aa2*aa1*(x1 - x2);

		// just take a guess if this happens
		if (temp == 0)
		{
			return x1 / 2.0;
		}

		dlib::matrix<double, 2, 1> temp2;
		temp2 = m*v / temp;
		const double a = temp2(0);
		const double b = temp2(1);

		temp = b*b - 3 * a*d0;
		if (temp < 0 || a == 0)
		{
			// This is probably a line so just pick the lowest point
			if (f0 < f_x2)
				return 0;
			else
				return x2;
		}
		temp = (-b + std::sqrt(temp)) / (3 * a);
		return put_in_range(0, x2, temp);
	}

	inline double poly_min_extrap(
		double f0,
		double d0,
		double f1
		)
	{
		const double temp = 2 * (f1 - f0 - d0);
		if (std::abs(temp) <= d0*std::numeric_limits<double>::epsilon())
			return 0.5;

		const double alpha = -d0 / temp;

		// now make sure the minimum is within the allowed range of (0,1) 
		return put_in_range(0, 1, alpha);
	}
//
//private:
//	std::ofstream output1;
//	std::ofstream output2;

};