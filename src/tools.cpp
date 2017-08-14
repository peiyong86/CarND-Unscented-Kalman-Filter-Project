#include <iostream>
#include "tools.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
  /**
  TODO:
    * Calculate the RMSE here.
  */
  VectorXd rmse(4);
  VectorXd ee, eg;
  double n = double(estimations.size());
  for(std::vector<VectorXd>::size_type i = 0; i != estimations.size(); i++) {
    ee = estimations[i];
    eg = ground_truth[i];
    for(int j = 0; j < 4; j++){
      rmse(j) = rmse(j) + (ee(j) - eg(j))*(ee(j) - eg(j));
    }
  }
  for(int i = 0; i < 4; i++){
    rmse(i) = rmse(i) / n;
  }
  return rmse;
}
