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
    * Calculate the RMSE here.
  */
	VectorXd rmse(4);
	rmse << 0,0,0,0;
    int n=estimations.size();

	// check the validity of the following inputs:
	//  * the estimation vector size should not be zero
	//  * the estimation vector size should equal ground truth vector size
	// ... your code here
    if (n==0 || n!=ground_truth.size()){
        cout<<"Error: sizes of estimation and ground truth don't match \n";
    }
	//accumulate squared residuals
	for(int i=0; i < n; ++i){
        VectorXd error=(estimations[i]-ground_truth[i]);
        error=error.array() *error.array();
        rmse+=error;
	}

	//calculate the mean
	// ... your code here
    rmse=rmse/n;
	//calculate the squared root
	// ... your code here
    rmse=rmse.array().sqrt();
	//return the result
	return rmse;
}
