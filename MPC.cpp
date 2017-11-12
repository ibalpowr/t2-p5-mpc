#include "MPC.h"
#include <cppad/cppad.hpp>
#include <cppad/ipopt/solve.hpp>
#include "Eigen-3.3/Eigen/Core"
#include <math.h>

using CppAD::AD;  //to find jacobians and hessians

//set the timestep length and duration
//finite and undiscounted horizon = N * dt
size_t N = 10;
double dt = 0.1;
//set reference speed
double ref_v= 30;

//weights for each cost components
double w0_cte = 1;
double w1_epsi = 1;
double w2_v = 1;
double w3_steer = 10;//1;
double w4_a = 10;//1;
double w5_diff_steer = 10;//1;
double w6_diff_a = 10;//100;//10;//1;

//the length from front (bumper) to Center of Gravity of the vehicle
const double Lf = 2.67;

//assume reference values for cte and epsi are both zero
double ref_cte = 0;
double ref_epsi = 0;

//set starting indices for state variables and actuators in "vars" vector
size_t x_start = 0;
size_t y_start = x_start + N;
size_t psi_start = y_start + N;
size_t v_start = psi_start + N;
size_t cte_start = v_start + N;
size_t epsi_start = cte_start + N;
size_t steer_start = epsi_start + N;
size_t a_start = steer_start + N -1;


//FG stands for f(x), g(x) linear functions in Ipopt package
//fg is a vector of the cost and constraints
class FG_eval {
  public:

  //coefficients for 3rd order polynomial
  Eigen::VectorXd coeffs;

  //class constructor
  FG_eval(Eigen::VectorXd coeffs) { this->coeffs = coeffs; }

  //create a type alias for CppAD vector
  typedef CPPAD_TESTVECTOR(AD<double>) ADvector;

  void operator()(ADvector& fg, const ADvector& vars) {

    //fg[0] is chosen to store the sum of all costs
    //initialize it to 0
    fg[0] = 0;

    for (unsigned int t = 0; t < N; ++t) {
      //1st cost: location ... represented by cte
      fg[0] += w0_cte * CppAD::pow(vars[cte_start + t] - ref_cte, 2);
      //2nd cost: vehicle heading
      fg[0] += w1_epsi * CppAD::pow(vars[epsi_start + t] - ref_epsi, 2);
    }

    for (unsigned int t = 0; t < N; ++t) {
		  //3rd cost: for stop-and-go situation
      fg[0] += w2_v * CppAD::pow(vars[v_start + t] - ref_v, 2);
    }

    for (unsigned int t = 0; t < N - 1; ++t) {
      //4th cost: for avoiding sudden jerk in steering
      fg[0] += w3_steer * CppAD::pow(vars[steer_start + t], 2);
      //5th cost: for avoiding sudden jerk in acceleration
      fg[0] += w4_a * CppAD::pow(vars[a_start + t], 2);
    }


		//6th and 7th costs: for differences of actuators values 
    //between adjacent look-ahead steps
    for (unsigned int t = 0; t < N - 2; ++t) {
      fg[0] += w5_diff_steer * CppAD::pow(vars[steer_start + t + 1] - 
               vars[steer_start + t], 2);
      fg[0] += w6_diff_a * CppAD::pow(vars[a_start + t + 1] - 
               vars[a_start + t], 2);
    }

    //adjust starting indices in fg vector for six state variables
    fg[1 + x_start] = vars[x_start];
    fg[1 + y_start] = vars[y_start];
    fg[1 + psi_start] = vars[psi_start];
    fg[1 + v_start] = vars[v_start];
    fg[1 + cte_start] = vars[cte_start];
    fg[1 + epsi_start] = vars[epsi_start];


    //set up model equations for each state variable update
    for (unsigned int t = 1; t < N; t++) {

      //the state at time t+1 .
      AD<double> x1 = vars[x_start + t ];
      AD<double> y1 = vars[y_start + t];
      AD<double> psi1 = vars[psi_start + t];
      AD<double> v1 = vars[v_start + t];
      AD<double> cte1 = vars[cte_start + t];

      //the state at time t.
      AD<double> x0 = vars[x_start + t-1];
      AD<double> y0 = vars[y_start + t -1];
      AD<double> psi0 = vars[psi_start + t-1];
      AD<double> v0 = vars[v_start + t-1];
      AD<double> cte0 = vars[cte_start + t-1];

      //recall the equations for the model:

      //for position x
      // x(t+1) = x(t) + v(t) * cos(psi(t)) * dt
      fg[1 + x_start + t] = x1 - (x0 + v0 * CppAD::cos(psi0) * dt);

      //for position y
      // y(t+1) = y(t) + v(t) * sin(psi(t)) * dt
      fg[1 + y_start + t] = y1 - (y0 + v0 * CppAD::sin(psi0) * dt);

      //for vehicle heading psi
      // psi(t+1) = psi(t) + v(t) / Lf * steer(t) * dt
      AD<double> delta0 = vars[steer_start + t -1];
      fg[1 + psi_start + t] = psi1 - (psi0 + v0 * delta0 / Lf * dt);

      //for speed v
      // v(t+1) = v(t) + a(t) * dt
      AD<double> a0 = vars[a_start + t-1];
      fg[1 + v_start + t] = v1 - (v0 + a0 * dt);

      //for cross track error cte
      // cte(t+1) = poly(x(t)) - y(t) + v(t) * sin(epsi(t)) * dt
      AD<double> epsi0 = vars[epsi_start + t-1];
      AD<double> f0 = coeffs[0] + coeffs[1] * x0 + 
                      coeffs[2] * x0 * x0 + coeffs[3] * x0 * x0 * x0;
      fg[1 + cte_start + t] = cte1 - ((f0 - y0) + 
                                     (v0 * CppAD::sin(epsi0) * dt));

      //for vehicle heading error epsi
      // epsi(t+1) = epsi(t) - psi_des(t) + v(t) / Lf * steer(t) * dt
      AD<double> epsi1 = vars[epsi_start + t];
      AD<double> psides0 = CppAD::atan(3*coeffs[3] * x0 * x0 + 
                           2 * coeffs[2] * x0 + coeffs[1]);
      fg[1 + epsi_start + t] = epsi1 - ((psi0 - psides0) + 
                               v0 * delta0 / Lf * dt);
    }
  }
};

//MPC class definition implementation
MPC::MPC() {}
MPC::~MPC() {}
						
//an important method
vector<double> MPC::Solve(Eigen::VectorXd state, Eigen::VectorXd coeffs) {
  bool ok = true;
  typedef CPPAD_TESTVECTOR(double) Dvector;

  // options for IPOPT solver
  std::string options;
  // Uncomment this if you'd like more print information
  options += "Integer print_level  0\n";
  // NOTE: Setting sparse to true allows the solver to take advantage
  // of sparse routines, this makes the computation MUCH FASTER. If you
  // can uncomment 1 of these and see if it makes a difference or not but
  // if you uncomment both the computation time should go up in orders of
  // magnitude.
  options += "Sparse  true        forward\n";
  options += "Sparse  true        reverse\n";
  // NOTE: Currently the solver has a maximum time limit of 0.5 seconds.
  // Change this as you see fit, especially in large N values
  options += "Numeric max_cpu_time          0.5\n";

	//set the length for vector vars
  //6 state variables (x, y, psi, v, cte, epsi) and 2 actuators
	//N steps ahead ... so in between (N - 1) actuations
  size_t n_vars = 6 * N + 2 * (N - 1);

	//set the number of constraints
  size_t n_constraints = 6 * N;

  //initialize vars
  //should be 0 except the initial state
  Dvector vars(n_vars);	
  for (unsigned int i = 0; i < n_vars; i++) {
    vars[i] = 0;
  }

  //state variables in vars have no lower or upper bounds
  Dvector vars_lowerbound(n_vars);
  Dvector vars_upperbound(n_vars);
  for (unsigned int i = 0; i < steer_start; i++) {
    vars_lowerbound[i] = -1.0e19;  //could be any very large number
    vars_upperbound[i] = 1.0e19;
  }

  //set lower and upper bounds for steer actuator in vars
  //[-25 deg, 25 deg]
  for (unsigned int i = steer_start; i < a_start; i++) {
    vars_lowerbound[i] = -25.0 * M_PI / 180.0; // in unit of [radian]
    vars_upperbound[i] = 25.0 * M_PI / 180.0;
  }

  //set lower and upper bounds for acceleration actuator in vars
  //[-1.0, 1.0]
  for (unsigned int i = a_start; i < n_vars; i++) {
    vars_lowerbound[i] = -1.0;
    vars_upperbound[i] = 1.0;
  }
	
  //set constraints for the Ipopt solver
  Dvector constraints_lowerbound(n_constraints);
  Dvector constraints_upperbound(n_constraints);
  for (unsigned int i = 0; i < n_constraints; i++) {
    constraints_lowerbound[i] = 0;
  constraints_upperbound[i] = 0;
  }

  //assign initial value to each state variable
  double x = state[0];
  double y = state[1];
  double psi = state[2];
  double v = state[3];
  double cte = state[4];
  double epsi = state[5];

  //for the first state
  constraints_lowerbound[x_start] = x;
  constraints_lowerbound[y_start] = y;
  constraints_lowerbound[psi_start] = psi;
  constraints_lowerbound[v_start] = v;
  constraints_lowerbound[cte_start] = cte;
  constraints_lowerbound[epsi_start] = epsi;

  constraints_upperbound[x_start] = x;
  constraints_upperbound[y_start] = y;
  constraints_upperbound[psi_start] = psi;
  constraints_upperbound[v_start] = v;
  constraints_upperbound[cte_start] = cte;
  constraints_upperbound[epsi_start] = epsi;

  //instantiate an object that defines costs and constraints
  FG_eval fg_eval(coeffs);


  //place to return solution
  CppAD::ipopt::solve_result<Dvector> solution;

/*  this is the KEY
*/	
  //solve the problem
  CppAD::ipopt::solve<Dvector, FG_eval>(
          options, vars, vars_lowerbound, vars_upperbound, 
          constraints_lowerbound, constraints_upperbound, fg_eval, solution);

  //check some of the solution values
  ok &= solution.status == CppAD::ipopt::solve_result<Dvector>::success;

  //cost
  auto cost = solution.obj_value;
  std::cout << "Cost " << cost << std::endl;

  //return MPC predicted steering angle and throttle value
  vector<double> result;
  result.push_back(solution.x[steer_start]);
  result.push_back(solution.x[a_start]);

  //for display MPC predicted trajectory
  for (unsigned int i = 0; i < N - 1; i++) {
    result.push_back(solution.x[x_start + i + 1]);
    result.push_back(solution.x[y_start + i + 1]);
  }

  return result;
}
