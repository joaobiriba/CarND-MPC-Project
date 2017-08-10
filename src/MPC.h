#ifndef MPC_H
#define MPC_H

#include <vector>
#include "Eigen-3.3/Eigen/Core"

using namespace std;
//GLA
const size_t N = 12;
const double dt = 0.05;
const int latency_ind = 2; //latency  in units of dt (100ms)

struct Solution {

		vector<double> positionX;
		vector<double> positionY;
		vector<double> psi;
		vector<double> velocity;
		vector<double> cte;
		vector<double> epsi;
		vector<double> steeringAngle;
		vector<double> acceleration;
};

//-GLA

class MPC {
 public:
  MPC();

  virtual ~MPC();

  // Solve the model given an initial state and polynomial coefficients.
  // Return the first actuatotions.
  //vector<double> Solve(Eigen::VectorXd state, Eigen::VectorXd coeffs);
  Solution Solve(Eigen::VectorXd state, Eigen::VectorXd coeffs);


  double delta_prev {0};
  double a_prev {0.1};
};

#endif /* MPC_H */
