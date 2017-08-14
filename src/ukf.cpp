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

  /**
  TODO:

  Complete the initialization. See ukf.h for other member properties.

  Hint: one or more values initialized above might be wildly off...
  */
  //set state dimension
  n_x = 5;

  //set augmented dimension
  n_aug = 7;

  is_initialized_ = false;
  Xsig_pred_ = MatrixXd(n_x, 2*n_aug + 1);
  Xsig_aug = MatrixXd(n_aug, 2 * n_aug + 1);
  lambda = 3 - n_aug;
}

UKF::~UKF() {}


void UKF::AugmentedSigmaPoints() {
  //create augmented mean vector
  VectorXd x_aug = VectorXd(7);

  //create augmented state covariance
  MatrixXd P_aug = MatrixXd(7, 7);

  //create sigma point matrix
  MatrixXd Xsig_aug = MatrixXd(n_aug, 2 * n_aug + 1);

  x_aug.head(5) = x;
  P_aug.topLeftCorner(5, 5) = P;
  P_aug(5, 5) = std_a * std_a;
  P_aug(6, 6) = std_yawdd * std_yawdd;
  MatrixXd A = P_aug.llt().matrixL();
  Xsig_aug.col(0) = x_aug;
  for(int i = 0; i < n_aug; i++){
      Xsig_aug.col(1 + i) = x_aug + sqrt(3) * A.col(i);
      Xsig_aug.col(8 + i) = x_aug - sqrt(3) * A.col(i);
  }
}

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
  if (!is_initialized_) {
    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
      float x = measurement_pack.raw_measurements_[0] * cos(measurement_pack.raw_measurements_[1] * PI);
      float y = measurement_pack.raw_measurements_[0] * sin(measurement_pack.raw_measurements_[1] * PI);
      ekf_.x_ << x, y, 0, 0;
    }
    else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
      ekf_.x_ << measurement_pack.raw_measurements_[0], measurement_pack.raw_measurements_[1], 0, 0;
    }

    // done initializing, no need to predict or update
    previous_timestamp_ = measurement_pack.timestamp_;
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
  //create matrix with predicted sigma points as columns
  double px, py, v, psi, psi_rate, noise_v, noise_psi;
  for(int i = 0; i < 2*n_aug+1; i++){
      px = Xsig_aug(0, i);
      py = Xsig_aug(1, i);
      v = Xsig_aug(2, i);
      psi = Xsig_aug(3, i);
      psi_rate = Xsig_aug(4, i);
      noise_v = Xsig_aug(5, i);
      noise_psi = Xsig_aug(6, i);
      
      if (psi_rate == 0){
          Xsig_pred(0, i) = px + v * cos(psi) * delta_t + 0.5 * delta_t * delta_t * cos(psi) * noise_v;
          Xsig_pred(1, i) = py + v * sin(psi) * delta_t + 0.5 * delta_t * delta_t * sin(psi) * noise_v;
          Xsig_pred(2, i) = v + delta_t * noise_v;
          Xsig_pred(3, i) = psi + psi_rate * delta_t + 0.5 * delta_t * delta_t * noise_psi;
          Xsig_pred(4, i) = psi_rate + delta_t * noise_psi;
      }
      else {
          Xsig_pred(0, i) = px + v / psi_rate * (sin(psi + psi_rate * delta_t) - sin(psi)) + 0.5 * delta_t * delta_t * cos(psi) * noise_v;
          Xsig_pred(1, i) = py + v / psi_rate * (-cos(psi + psi_rate * delta_t) + cos(psi)) + 0.5 * delta_t * delta_t * sin(psi) * noise_v;
          Xsig_pred(2, i) = v + delta_t * noise_v;
          Xsig_pred(3, i) = psi + psi_rate * delta_t + 0.5 * delta_t * delta_t * noise_psi;
          Xsig_pred(4, i) = psi_rate + delta_t * noise_psi;
      }
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
