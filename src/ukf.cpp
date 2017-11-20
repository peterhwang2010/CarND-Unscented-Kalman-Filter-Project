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
  std_a_ = 0.5; //30

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 0.15; //30

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

  /**
  TODO:

  Complete the initialization. See ukf.h for other member properties.

  Hint: one or more values initialized above might be wildly off...
  */
  is_initialized_ = false;

  // State dimension
  n_x_ = 5;

  // Augmented state dimension
  n_aug_ = 7;

  // Sigma point spreading parameter
  lambda_ = 3 - n_aug_;
    
  // Vector for weights
  weights_ = VectorXd(2*n_aug_+1);

  //define spreading parameter
  weights_(0) = lambda_ / (lambda_+n_aug_);
  for(int i=0; i<2*n_aug_; i++) weights_(i+1)=1/(2*(lambda_+n_aug_));  
  
}

UKF::~UKF() {}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Make sure you switch between lidar and radar
  measurements.
  */
  if(is_initialized_){
    double delta_t = double(meas_package.timestamp_ - time_us_) / 1000000;
    
    if (meas_package.sensor_type_ == MeasurementPackage::LASER && use_laser_){
      // Predict
      Prediction(delta_t);    
      
      // Update
      UpdateLidar(meas_package);   

      // Update timestamp
      time_us_ = meas_package.timestamp_;
    }
    if(meas_package.sensor_type_ == MeasurementPackage::RADAR && use_radar_){
        // Predict
        Prediction(delta_t);

        // Update
        UpdateRadar(meas_package);

        // Update timestamp
        time_us_ = meas_package.timestamp_;
    }
    
  }else{
    // Initialize
    time_us_ = meas_package.timestamp_;

    x_ = VectorXd(n_x_);
       
    if (meas_package.sensor_type_ == MeasurementPackage::LASER){

      // First Measurement
      x_ << meas_package.raw_measurements_[0], meas_package.raw_measurements_[1], 0, 0, 0;
    }else{
      double rho = meas_package.raw_measurements_[0]; 
      double theta = meas_package.raw_measurements_[1];
      double rhod = meas_package.raw_measurements_[2];

      // First measurement
      x_ << rho * cos(theta), rho * sin(theta), rhod, 0, 0;
    }

    // Initial covariance matrix
    P_ = MatrixXd(5,5);
    P_.fill(0.0);
    P_(0,0)=1;
    P_(1,1)=1;
    P_(2,2)=10;
    P_(3,3)=10;
    P_(4,4)=10;
    
    is_initialized_ = true;
    
    return;
  } 
}

/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(double delta_t) {
  /**
  TODO:

  Complete this function! Estimate the object's location. Modify the state
  vector, x_. Predict sigma points, the state, and the state covariance matrix.
  */

  // Augmented mean vector
  VectorXd x_aug = VectorXd(n_aug_);

  // Augmented state covariance
  MatrixXd P_aug = MatrixXd(n_aug_, n_aug_);

  // Sigma point matrix
  MatrixXd Xsig_aug = MatrixXd(n_aug_, 2 * n_aug_ + 1);

  // Augmented mean state
  x_aug.head(n_x_) = x_;
  x_aug(5) = 0;
  x_aug(6) = 0;
  
  P_aug.fill(0.0);
  P_aug.topLeftCorner( n_x_, n_x_ ) = P_;
  P_aug(5,5) = std_a_ * std_a_;
  P_aug(6,6) = std_yawdd_ * std_yawdd_;

  // Square root matrix
  MatrixXd A = P_aug.llt().matrixL();
  MatrixXd root =  sqrt(lambda_ + n_aug_) * A;
    
  // Augmented sigma points
  Xsig_aug.col(0) = x_aug;
  for(int i = 0; i < n_aug_; i++){
    Xsig_aug.col(1 + i) = x_aug + root.col(i);
    Xsig_aug.col(1 + n_aug_ + i) = x_aug - root.col(i);
  }

  Xsig_pred_ = MatrixXd(n_x_, 2 * n_aug_ + 1);

  double dt2 = delta_t * delta_t;
  for(int i = 0; i < 2 * n_aug_ + 1; i++){
    double p_x = Xsig_aug(0,i);
    double p_y = Xsig_aug(1,i);
    double v = Xsig_aug(2,i);
    double yaw = Xsig_aug(3,i);
    double yaw_d = Xsig_aug(4,i);
    double nu_a = Xsig_aug(5,i);
    double nu_yawdd = Xsig_aug(6,i);
    
  // Predict sigma points
    if(fabs(yaw_d) < 0.001){
        Xsig_pred_(0,i) = p_x + v * cos(yaw) * delta_t + 0.5 * dt2 * cos(yaw) * nu_a;
        Xsig_pred_(1,i) = p_y + v * sin(yaw) * delta_t + 0.5 * dt2 * sin(yaw) * nu_a;
    }else{
        Xsig_pred_(0,i) = p_x + (v / yaw_d) * (sin(yaw + yaw_d * delta_t) - sin(yaw)) + 0.5 * dt2 * cos(yaw) * nu_a;
        Xsig_pred_(1,i) = p_y + (v / yaw_d) * (-cos(yaw + yaw_d * delta_t) + cos(yaw)) + 0.5 * dt2 * sin(yaw) * nu_a;
    }
    Xsig_pred_(2,i) = v + delta_t * nu_a;
    Xsig_pred_(3,i) = yaw + yaw_d * delta_t + 0.5 * dt2 * nu_yawdd;
    Xsig_pred_(4,i) = yaw_d + delta_t * nu_yawdd;
  }
 
  // Predict state mean
  x_ = VectorXd(5);
  x_.fill(0.0);
  for(int j=0; j < 2 * n_aug_ + 1; j++){
      x_ += weights_(j) * Xsig_pred_.col(j);
  }

  // Predict state covariance matrix
  P_ = MatrixXd(5,5);
  P_.fill(0.0);
  for(int j=0; j < 2 * n_aug_ + 1; j++){
      VectorXd diff = Xsig_pred_.col(j) - x_;
      diff(3) = Normalize(diff(3));
      P_ += weights_(j) * diff * diff.transpose();      
  }
  return;  
  
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

  VectorXd z = meas_package.raw_measurements_;
  
    // Measurement dimension
    int n_z = 2;
  
    // Matrix for sigma points in measurement space
    MatrixXd Zsig = MatrixXd(n_z, 2 * n_aug_ + 1);
  
    //transform sigma points into measurement space
    for (int i = 0; i < 2 * n_aug_ + 1; i++) {
      // Measurement model
      Zsig(0, i) = Xsig_pred_(0, i);
      Zsig(1, i) = Xsig_pred_(1, i);
    }
  
    // Mean predicted measurement
    VectorXd z_pred = VectorXd(n_z);
    z_pred.fill(0.0);
    for (int i = 0; i < 2 * n_aug_ + 1; i++) {
      z_pred = z_pred + weights_(i) * Zsig.col(i);
    }
  
    // Covariance matrix S
    MatrixXd S = MatrixXd(n_z, n_z);
    S.fill(0.0);
    for (int i = 0; i < 2 * n_aug_ + 1; i++) {
      // Residual
      VectorXd z_diff = Zsig.col(i) - z_pred;
      S = S + weights_(i) * z_diff * z_diff.transpose();
    }
  
    // Measurement noise covariance matrix
    MatrixXd R = MatrixXd(n_z, n_z);
    R << std_laspx_*std_laspx_, 0,
         0, std_laspy_*std_laspy_;
    S = S + R;
  
    // Matrix for cross correlation Tc
    MatrixXd Tc = MatrixXd(n_x_, n_z);
  
    // Cross correlation matrix
    Tc.fill(0.0);
    for (int i = 0; i < 2 * n_aug_ + 1; i++) {
  
      // Residual
      VectorXd z_diff = Zsig.col(i) - z_pred;
  
      // State difference
      VectorXd x_diff = Xsig_pred_.col(i) - x_;
  
      Tc = Tc + weights_(i) * x_diff * z_diff.transpose();
    }
  
    // Kalman gain K;
    MatrixXd K = Tc * S.inverse();
  
    // Residual
    VectorXd z_diff = z - z_pred;
  
    // NIS
    NIS_L_ = z_diff.transpose() * S.inverse() * z_diff;

    x_ = x_ + K * z_diff;
    P_ = P_ - K*S*K.transpose();
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
  int n_z = 3;
 
  // Matrix for sigma points in measurement space
  MatrixXd Zsig = MatrixXd(n_z, 2 * n_aug_ + 1);

  // Mean predicted measurement
  VectorXd z_pred = VectorXd(n_z);
  
  // Measurement covariance matrix S
  MatrixXd S = MatrixXd(n_z,n_z);

  z_pred.fill(0.0);
  S.fill(0.0);
  S(0,0)=std_radr_*std_radr_;
  S(1,1)=std_radphi_*std_radphi_;
  S(2,2)=std_radrd_*std_radrd_;
  
  for (int i=0; i<2*n_aug_+1; i++) {  
    double p_x = Xsig_pred_(0,i);
    double p_y = Xsig_pred_(1,i);
    double v = Xsig_pred_(2,i);
    double yaw = Xsig_pred_(3,i);
    double yaw_d = Xsig_pred_(4,i);
    
    double rho = sqrt( p_x * p_x + p_y * p_y);
    double psi = atan2(p_y, p_x);
    double rho_d = (p_x * cos (yaw) * v + p_y * sin(yaw) * v ) / rho;
    
    if(fabs(rho) < 0.0001){
      cout << "Division by zero error" << endl;
      cin >> rho;
    }
    
    Zsig(0,i)=rho;
    Zsig(1,i)=psi;
    Zsig(2,i)=rho_d;
    
    z_pred += weights_(i) * Zsig.col(i);
  }

  for (int i=0; i < 2 * n_aug_ + 1; i++) {  
      VectorXd residual = Zsig.col(i) - z_pred;
      residual(1) = Normalize(residual(1));    
      S += weights_(i) * residual * residual.transpose();
  }
  
  // Matrix for cross correlation Tc
  MatrixXd Tc = MatrixXd(n_x_, n_z);

  // Cross correlation matrix
  Tc.fill(0.0);
  for (int i=0; i < 2 * n_aug_+1; i++) {
      VectorXd x_diff = Xsig_pred_.col(i) - x_;
      x_diff(3) = Normalize(x_diff(3));

      VectorXd z_diff = Zsig.col(i) - z_pred;
      z_diff(1) = Normalize(z_diff(1));

      Tc += weights_(i) * x_diff * z_diff.transpose();
  }
  
  // Kalman gain K;
  MatrixXd K = Tc * S.inverse();
  
  VectorXd z = VectorXd(3);
  z(0) = meas_package.raw_measurements_(0);
  z(1) = meas_package.raw_measurements_(1);
  z(2) = meas_package.raw_measurements_(2);

  VectorXd z_diff = z - z_pred;
  z_diff(1) = Normalize(z_diff(1));
  
  x_ = x_ + K * z_diff;
  P_ = P_ - K * S * K.transpose();
  NIS_R_ = z_diff.transpose() * S.inverse() * z_diff;  
    
  return;
}

double UKF::Normalize(double x) {
  while(x < -M_PI) x += 2 * M_PI;
  while(x >  M_PI) x -= 2 * M_PI;
  
  return x;
}
