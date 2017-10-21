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
  is_initialized_=false;
  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = true;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;

  // initial state vector
  x_ = VectorXd(5);

  // initial covariance matrix
  P_ = MatrixXd(5, 5);

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 1.5;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 0.57;

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
  Xsig_pred_ = MatrixXd(n_x_, n_sig_); 
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
// void UKF::NormAng(double& ang){
//   if (ang>M_PI) ang-=2.0*M_PI;
//   if (ang<M_PI) ang+=2.0*M_PI;
// }
void UKF::NormAng(double *ang) {
    while (*ang > M_PI) *ang -= 2. * M_PI;
    while (*ang < -M_PI) *ang += 2. * M_PI;
}


void UKF::ProcessMeasurement(MeasurementPackage measurement_pack) {
  //initializing state with first measurement 
  if (!is_initialized_){
    P_=Eigen::MatrixXd::Identity(5,5);
    float px,py,v;
    if(measurement_pack.sensor_type_==MeasurementPackage::RADAR){
      float rho=measurement_pack.raw_measurements_[0];
      float phi=measurement_pack.raw_measurements_[1];
      float rho_dot=measurement_pack.raw_measurements_[2];
      px=rho*cos(phi);
      py=rho*sin(phi);
      v=rho_dot;
    }
    else if(measurement_pack.sensor_type_==MeasurementPackage::LASER){
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

//this is fine
void UKF::Prediction(double dt) {
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
    double arg= yaw+ (yawd*dt);
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
    VectorXd x_diff=Xsig_pred_.col(i)-x_;
    NormAng(&x_diff(3));
    P_=P_ + weights_(i)*x_diff*x_diff.transpose();
  }
}

void UKF::UpdateLidar(MeasurementPackage measurement_pack) {
  int n_z=2;
  //Transform sigma points to measurement space
  MatrixXd Zsig = Xsig_pred_.block(0,0,n_z,n_sig_);
  UpdateUKF(measurement_pack, Zsig, n_z);
}

void UKF::UpdateRadar(MeasurementPackage measurement_pack) {
  int n_z =3;
  MatrixXd Zsig= MatrixXd (n_z, n_sig_);
  //Transform sigma points to measurement space
  for (int i=0;i< n_sig_; i++){
    double p_x = Xsig_pred_(0,i);
    double p_y = Xsig_pred_(1,i);
    double v= Xsig_pred_(2,i);
    double yaw= Xsig_pred_(3,i);
    double v1=cos(yaw)*v;
    double v2=sin(yaw)*v;
    //measuremnt model
    Zsig(0,i) = sqrt(p_x*p_x + p_y*p_y);  //r
    Zsig(1,i) = atan2(p_y,p_x); //phi
    Zsig(2,i) = (p_x*v1+ p_y*v2)/Zsig(0,i); // r_dot
  }
  UpdateUKF(measurement_pack, Zsig, n_z);
}
void UKF::UpdateUKF(MeasurementPackage measurement_pack, MatrixXd Zsig, int n_z){
  //mean predicted measurement
  VectorXd z_pred = VectorXd(n_z);
  z_pred= Zsig* weights_;
  //measurement covariance matrix
  MatrixXd S = MatrixXd(n_z, n_z);
  S.fill(0.0);
  for (int i=0; i< n_sig_;i++){
    VectorXd z_diff = Zsig.col(i) - z_pred;
    NormAng(&z_diff(1));
    S=S+weights_(i)*z_diff*z_diff.transpose();
  }
  //add measurement noise covariance matrix
  MatrixXd R= MatrixXd(n_z, n_z);
  if (measurement_pack.sensor_type_==MeasurementPackage::RADAR){
    R=R_radar_;
  }
  else if(measurement_pack.sensor_type_==MeasurementPackage::LASER){
    R=R_lidar_;
  }
  S=S+R;
  // Cross Correlation matrix
  MatrixXd Tc=  MatrixXd(n_x_, n_z);
  Tc.fill(0.0);
  for(int i=0;i< n_sig_; i++){
    VectorXd z_diff= Zsig.col(i)-z_pred;
    if(measurement_pack.sensor_type_ ==MeasurementPackage::RADAR){
      NormAng(&z_diff(1));
    }
    VectorXd x_diff= Xsig_pred_.col(i)-x_;
    NormAng(&x_diff(3));
    Tc=Tc+weights_(i)*x_diff*z_diff.transpose();
  }
  //measurement
  VectorXd z= measurement_pack.raw_measurements_;
  //kalman gain K;
  MatrixXd K= Tc*S.inverse();
  VectorXd z_diff = z-z_pred;
  if (measurement_pack.sensor_type_ == MeasurementPackage:: RADAR){
    NormAng(&z_diff(1));
  }
  //update mean and covariance
  x_ = x_ + K * z_diff;
  P_ = P_ - K * S* K.transpose();
  //calculate NIS
  if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR){
    NIS_radar_ = z.transpose()* S.inverse() *z;
  }
  else if(measurement_pack.sensor_type_ == MeasurementPackage::LASER){
    NIS_laser_= z.transpose()*S.inverse()*z;
  }
  cout<<"state:"<<x_<<"\n";
  cout<<"covariance:"<<P_<<"\n";
}
