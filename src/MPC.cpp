#include "MPC.h"
#include <cppad/cppad.hpp>
#include <cppad/ipopt/solve.hpp>
#include "Eigen-3.3/Eigen/Core"

using CppAD::AD;

// This value assumes the model presented in the classroom is used.
//
// It was obtained by measuring the radius formed by running the vehicle in the
// simulator around in a circle with a constant steering angle and velocity on a
// flat terrain.
//
// Lf was tuned until the the radius formed by the simulating the model
// presented in the classroom matched the previous radius.
//
// This is the length from front to CoG that has a similar radius.
const double Lf = 2.67;

double reference_cte = 0;
double reference_epsi = 0;
double reference_velocity = 85;


// variables indexes in vector
size_t start_x = 0;
size_t start_y = start_x + N;
size_t start_psi = start_y + N;
size_t v_start = start_psi + N;
size_t cte_start = v_start + N;
size_t start_epsi = cte_start + N;
size_t start_delta = start_epsi + N;
size_t start_acc = start_delta + N - 1;


class FG_eval {
 public:
  // Fitted polynomial coefficients
  Eigen::VectorXd coeffs;
  vector<double> previous_actuations;

  FG_eval(Eigen::VectorXd coeffs, vector<double> previous_actuations) {
    this->coeffs = coeffs;
    this->previous_actuations = previous_actuations;
 }

  typedef CPPAD_TESTVECTOR(AD<double>) ADvector;

  void operator()(ADvector& fg, const ADvector& vars) {
      // `fg` a vector of the cost constraints, `vars` is a vector of variable values (state & actuators)
     fg[0] = 0;

     for (int i = 0; i < N; i++) {
       // trajectory
       fg[0] += CppAD::pow(vars[cte_start + i] - reference_cte, 2);
       fg[0] += CppAD::pow(vars[start_epsi + i] - reference_epsi, 2);
       fg[0] += CppAD::pow(vars[v_start + i] - reference_velocity, 2);
     }

     // actuators factors
     for (int i = 0; i < N - 1; i++) {
       fg[0] += CppAD::pow(vars[start_delta + i], 2);
       fg[0] += 10*CppAD::pow(vars[start_acc + i], 2);
     }

     // gap between sequential actuations factors
     for (int i = 0; i < N - 2; i++) {
       fg[0] += 600*CppAD::pow(vars[start_delta + i + 1] - vars[start_delta + i], 2);
       fg[0] += CppAD::pow(vars[start_acc + i + 1] - vars[start_acc + i], 2);
     }


     fg[1 + start_x] = vars[start_x];
     fg[1 + start_y] = vars[start_y];
     fg[1 + start_psi] = vars[start_psi];
     fg[1 + v_start] = vars[v_start];
     fg[1 + cte_start] = vars[cte_start];
     fg[1 + start_epsi] = vars[start_epsi];

     // following constraints
     for (int i = 0; i < N - 1; i++) {
       // State at t+1 .
       AD<double> x1 = vars[start_x + i + 1];
       AD<double> y1 = vars[start_y + i + 1];
       AD<double> psi1 = vars[start_psi + i + 1];
       AD<double> v1 = vars[v_start + i + 1];
       AD<double> cte1 = vars[cte_start + i + 1];
       AD<double> epsi1 = vars[start_epsi + i + 1];

       // State at t.
       AD<double> x0 = vars[start_x + i];
       AD<double> y0 = vars[start_y + i];
       AD<double> psi0 = vars[start_psi + i];
       AD<double> v0 = vars[v_start + i];
       AD<double> cte0 = vars[cte_start + i];
       AD<double> epsi0 = vars[start_epsi + i];

       // Actuation at time t.
       AD<double> delta0 = vars[start_delta + i];
       AD<double> a0 = vars[start_acc + i];

       AD<double> f0 = coeffs[0] + coeffs[1] * x0 + coeffs[2]*x0*x0 + coeffs[3]*x0*x0*x0;
       AD<double> psides0 = CppAD::atan(coeffs[1]+2*coeffs[2]*x0 + 3 * coeffs[3]*x0*x0);


       // Equations of the model:
       fg[2 + start_x + i] = x1 - (x0 + v0 * CppAD::cos(psi0) * dt);
       fg[2 + start_y + i] = y1 - (y0 + v0 * CppAD::sin(psi0) * dt);
       fg[2 + start_psi + i] = psi1 - (psi0 + v0 * delta0 / Lf * dt);
       fg[2 + v_start + i] = v1 - (v0 + a0 * dt);
       fg[2 + cte_start + i] =
           cte1 - ((f0 - y0) + (v0 * CppAD::sin(epsi0) * dt));
       fg[2 + start_epsi + i] =
           epsi1 - ((psi0 - psides0) + v0 * delta0 / Lf * dt);
     }
   }
};

//
// MPC class definition implementation.
//
MPC::MPC() {}
MPC::~MPC() {}

Solution MPC::Solve(Eigen::VectorXd x0, Eigen::VectorXd coeffs) {
//vector<double> MPC::Solve(Eigen::VectorXd state, Eigen::VectorXd coeffs) {
  bool ok = true;
  typedef CPPAD_TESTVECTOR(double) Dvector;


  double x = x0[0];
  double y = x0[1];
  double psi = x0[2];
  double v = x0[3];
  double cte = x0[4];
  double epsi = x0[5];


  // Set the number of model variables (includes both states and inputs).
  // For example: If the state is a 4 element vector, the actuators is a 2
  // element vector and there are 10 timesteps. The number of variables is:

  // N timesteps == N - 1 actuations
  size_t n_vars = N * 6 + (N - 1) * 2;
  // Number of constraints
  size_t n_constraints = N * 6;


  Dvector vars(n_vars);
  for (int i = 0; i < n_vars; i++) {
    vars[i] = 0;
  }

  // Variables initialization
  vars[start_x] = x;
  vars[start_y] = y;
  vars[start_psi] = psi;
  vars[v_start] = v;
  vars[cte_start] = cte;
  vars[start_epsi] = epsi;

  Dvector vars_lowerbound(n_vars);
  Dvector vars_upperbound(n_vars);

  for (int i = 0; i < start_delta; i++) {
    vars_lowerbound[i] = -1.0e19;
    vars_upperbound[i] = 1.0e19;
  }

  for (int i = start_delta; i < start_acc; i++) {
    vars_lowerbound[i] = -0.436332;
    vars_upperbound[i] = 0.436332;
  }

  // constraint delta to be the previous control for the latency time
  for (int i = start_delta; i < start_delta + latency_ind; i++) {
    vars_lowerbound[i] = delta_prev;
    vars_upperbound[i] = delta_prev;
  }

  // Acceleration/decceleration upper and lower limits.
  for (int i = start_acc; i < n_vars; i++) {
    vars_lowerbound[i] = -1.0;
    vars_upperbound[i] =  1.0;
  }

  for (int i = start_acc; i < start_acc+latency_ind; i++) {
    vars_lowerbound[i] = a_prev;
    vars_upperbound[i] = a_prev;
  }


  Dvector constraints_lowerbound(n_constraints);
  Dvector constraints_upperbound(n_constraints);
  for (int i = 0; i < n_constraints; i++) {
    constraints_lowerbound[i] = 0;
    constraints_upperbound[i] = 0;
  }
  constraints_lowerbound[start_x] = x;
  constraints_lowerbound[start_y] = y;
  constraints_lowerbound[start_psi] = psi;
  constraints_lowerbound[v_start] = v;
  constraints_lowerbound[cte_start] = cte;
  constraints_lowerbound[start_epsi] = epsi;

  constraints_upperbound[start_x] = x;
  constraints_upperbound[start_y] = y;
  constraints_upperbound[start_psi] = psi;
  constraints_upperbound[v_start] = v;
  constraints_upperbound[cte_start] = cte;
  constraints_upperbound[start_epsi] = epsi;

  vector<double> previous_actuations = {delta_prev,a_prev};
  FG_eval fg_eval(coeffs,previous_actuations);

  // options
  std::string options;
  options += "Integer print_level  0\n";
  options += "Sparse  true        forward\n";
  options += "Sparse  true        reverse\n";
  options += "Numeric max_cpu_time          0.05\n";

  // place to return solution
  CppAD::ipopt::solve_result<Dvector> solution;

  // solve the problem
  CppAD::ipopt::solve<Dvector, FG_eval>(
      options, vars, vars_lowerbound, vars_upperbound, constraints_lowerbound,
      constraints_upperbound, fg_eval, solution);

  // Check some of the solution values
  ok &= solution.status == CppAD::ipopt::solve_result<Dvector>::success;

  Solution sol;
  for (auto i = 0; i < N-1 ; i++){
  	cout << i << ": " << "solution.x[start_x+i]: " << solution.x[start_x+i] << "solution.x[start_y+i]: " << solution.x[start_y+i] << endl;
  	sol.positionX.push_back(solution.x[start_x+i]);
  	sol.positionY.push_back(solution.x[start_y+i]);
  	sol.steeringAngle.push_back(solution.x[start_delta+i]);
  	sol.acceleration.push_back(solution.x[start_acc+i]);
  }

  // Cost
  auto cost = solution.obj_value;
  std::cout << "Cost " << cost << std::endl;

  return sol;
}
