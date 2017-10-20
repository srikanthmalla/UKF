#include "ukf.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/**
 * Initializes Unscented Kalman filter
 */
UKF::UKF() {
  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = true;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;

  // initial state vector
  x_ = VectorXd(5);

  // initial covariance matrix
  P_ = MatrixXd(5, 5);

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 30;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 30;

  // Laser measurement noise standard deviation position1 in m
  std_laspx_ = 0.15;

  // Laser measurement noise standard deviation position2 in m
  std_laspy_ = 0.15;

  // Radar measurement noise standard deviation radius in m
  std_radr_ = 0.3;

  // Radar measurement noise standard deviation angle in rad
  std_radphi_ = 0.03;

  // Radar measurement noise standard deviation radius change in m/s
  std_radrd_ = 0.3;

  //intializing
  n_x_=x_.size();
  n_aug_=n_x_+2;
  n_sig_=2*n_aug_+1;
  lambda_=3-n_aug_;
  weights_=VectorXd(2*n_aug_+1);
  weights_(0)= lambda_/(lambda_+n_aug_);
  double weight=0.5/(n_aug_+lambda_);
  for(int i=1;i<2*n_aug_+1;i++){
    weights_(i)=weight;
  }
  R_radar_=MatrixXd(3,3);
  R_radar_<< std_radr_*std_radr_,0,0,
            0, std_radphi_*std_radphi_,0,
            0,0,std_radrd_*std_radrd_;
  R_lidar_=MatrixXd(2,2);
  R_lidar_<<std_laspx_*std_laspx_,0,
            0,std_laspy_*std_laspy_;

}

//destructor
UKF::~UKF() {}

//normalizing angle
void UKF::NormAng(double& ang){
  if (ang>M_PI) ang-=2.0*M_PI;
  if (ang<M_PI) ang+=2.0*M_PI;
}
void UKF::ProcessMeasurement(MeasurementPackage measurement_pack) {
  //initializing state with first measurement 
  if (!is_initialized_){
    P_=Eigen::MatrixXd::Identity(5,5);
    float px,py,v;
    if(measurement_pack.sensor_type_=MeasurementPackage::RADAR){
      float rho=measurement_pack.raw_measurements_[0];
      float phi=measurement_pack.raw_measurements_[1];
      float rho_dot=measurement_pack.raw_measurements_[2];
      px=rho*cos(phi);
      py=rho*sin(phi);
      v=rho_dot;
    }
    else if(measurement_pack.sensor_type_=MeasurementPackage::LASER){
      px=measurement_pack.raw_measurements_[0];
      py=measurement_pack.raw_measurements_[1];
      v=0;
      if (fabs(px)<0.001 and fabs(py)<0.001){
        px=0.001;
        py=0.001;
      }
    }
    x_<<px,py,v,0,0;
    time_us_=measurement_pack.timestamp_;
    is_initialized_=true;
    return;
  }
  else{
    double dt=(measurement_pack.timestamp_-time_us_)/1000000.0;
    time_us_=measurement_pack.timestamp_;
    Prediction(dt);//predict sigma points
    if (measurement_pack.sensor_type_==MeasurementPackage::RADAR && use_radar_){
      UpdateRadar(measurement_pack);
    }
    else if(measurement_pack.sensor_type_ ==MeasurementPackage::LASER && use_laser_){
      UpdateLidar(measurement_pack);
    }
  }
}

/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 */
void UKF::Prediction(double dt) {
  //TODO:
  double dt_2=dt*dt;
  VectorXd x_aug=VectorXd(n_aug_);
  MatrixXd P_aug=MatrixXd(n_aug_,n_aug_);
  MatrixXd Xsig_aug= MatrixXd(n_aug_,n_sig_);
  //augmenting state and covariance matrices with yawdd and a
  x_aug.fill(0.0);
  x_aug.head(n_x_)=x_;
  P_aug.fill(0.0);
  P_aug.topLeftCorner(n_x_,n_x_)=P_;
  P_aug(5,5)=std_a_*std_a_; 
  P_aug(6,6)=std_yawdd_*std_yawdd_;
  MatrixXd L = P_aug.llt().matrixL();
  //generating sigma points
  Xsig_aug.col(0)=x_aug;
  double scaling_factor=sqrt(lambda_+n_aug_);
  for(int i=0;i<n_aug_;i++){
    Xsig_aug.col(i+1) = x_aug + scaling_factor*L.col(i);
    Xsig_aug.col(i+1+n_aug_) = x_aug - scaling_factor*L.col(i);
  }
  // Predict sigma points using motion model
  for(int i=0;i<n_sig_;i++){
    double px=Xsig_aug(0,i);
    double py=Xsig_aug(1,i);
    double v=Xsig_aug(2,i);
    double yaw=Xsig_aug(3,i);
    double yawd=Xsig_aug(4,i);
    double nu_a=Xsig_aug(5,i);
    double nu_yawdd=Xsig_aug(6,i);  
    double px_p,py_p;
    double arg= yaw+yaw*dt;
    if (fabs(yawd)>0.001){
      px_p=px + (v/yawd)*(sin(arg)-sin(yaw));
      py_p=py + (v/yawd)*(cos(yaw)-cos(arg));
    }
    else{
      px_p=px + v*dt*cos(yaw);
      py_p=py + v*dt*sin(yaw);
    }
    double v_p=v;
    double yaw_p=arg;
    double yawd_p=yawd;
    //Add noise
    px_p +=0.5*dt_2*nu_a*cos(yaw);
    py_p +=0.5*dt_2*nu_a*sin(yaw);
    v_p += nu_a*dt;
    yaw_p += 0.5*nu_yawdd*dt_2;
    yawd_p += nu_yawdd*dt;
    //write predicted siga points to Xsig_pred_
    Xsig_pred_(0,i)=px_p;
    Xsig_pred_(1,i)=py_p;
    Xsig_pred_(2,i)=v_p;
    Xsig_pred_(3,i)=yaw_p;
    Xsig_pred_(4,i)=yawd_p;
  }
  //predicting state mean from the sigma points and weights
  x_=Xsig_pred_ * weights_;
  //predicting state covariance using sigma points and weights
  P_.fill(0.0);
  for (int i=0;i<n_sig_;i++){
    VectorXd x_diff=Xsig_pred_.col(i)=x_;
    NormAng(x_diff(3));
    P_=P_+weights_(i)*x_diff*x_diff.transpose();
  }
}

/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLidar(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Use lidar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the lidar NIS.
  */
}

/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateRadar(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Use radar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the radar NIS.
  */
}
void UKF::UpdateUKF(MeasurementPackage measurement_pack, MatrixXd Zsig, int n_z){
}
