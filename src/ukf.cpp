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
  std_a_ = 1.55;//30; // NEED TO TUNE THIS VALUE

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 0.76;//30; // NEED TO TUNE THIS VALUE

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
  // State dimension
  n_x_ = 5;

  // Augmented state dimension
  n_aug_ = 7;

  // Sigma point spreading parameter
  lambda_= 3 - n_x_;

  // Weights of sigma points
  weights_ = VectorXd(2*n_aug_+1);
  double weight_0 = lambda_/(lambda_+n_aug_);
  weights_(0) = weight_0;
  for (int i=1; i<2*n_aug_+1; i++) {
    double weight = 0.5 / (n_aug_ + lambda_);
    weights_(i) = weight;
  }

  // predicted sigma points matrix
  Xsig_pred_ = MatrixXd(n_x_, 2 * n_aug_ + 1);

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
  /*****************************************************************************
  *  Initialization
  ****************************************************************************/
  if (!is_initialized_) {
    /**
    * Initialize the state x_ with the first measurement.
    * state vector: [pos1 pos2 vel_abs yaw_angle yaw_rate] in SI units and rad
    * Create the covariance matrix.
    */
    if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
      /**
      Convert radar from polar to cartesian coordinates and initialize state.
      */
      float rho = meas_package.raw_measurements_[0];
      float phi = meas_package.raw_measurements_[1];
      float norm_phi = atan2(sin(phi), cos(phi));

      float px = rho * cos(norm_phi);
      float py = rho * sin(norm_phi);
      x_ << px,
            py,
            3,
            0,
            0.95;
    } else if (meas_package.sensor_type_ == MeasurementPackage::LASER) {
      /**
      Initialize state.
      */
      x_ << meas_package.raw_measurements_[0],
            meas_package.raw_measurements_[1],
            3,
            0,
            0.95;
    }

    time_us_ = meas_package.timestamp_;


    // state covariance matrix P
    P_ << 0.06, 0, 0, 0, 0,
        0, 0.00001, 0, 0, 0,
        0, 0, 1, 0, 0,
        0, 0, 0, 0.08, 0,
        0, 0, 0, 0, 0.15;

    // done initializing, no need to predict or update
    is_initialized_ = true;
    return;
  }

  float delta_t = (meas_package.timestamp_ - time_us_) / 1000000.0;
  time_us_ = meas_package.timestamp_;
  Prediction(delta_t);

  if (meas_package.sensor_type_ == MeasurementPackage::RADAR && use_radar_) {
    UpdateRadar(meas_package);
  } else if (meas_package.sensor_type_ == MeasurementPackage::LASER && use_laser_) {
    UpdateLidar(meas_package);
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

  //create augmented mean vector
  VectorXd x_aug = VectorXd(7);

  //create augmented state covariance
  MatrixXd P_aug = MatrixXd(7, 7);

  //create sigma point matrix
  MatrixXd Xsig_aug = MatrixXd(n_aug_, 2 * n_aug_ + 1);

  //create augmented mean state
  x_aug.head(n_x_) = x_;
  x_aug(n_x_) = 0;
  x_aug(n_x_+1) = 0;

  //create augmented covariance matrix
  P_aug.fill(0.0);
  P_aug.topLeftCorner(n_x_, n_x_) = P_;
  MatrixXd Q = MatrixXd(2, 2);
  Q << std_a_*std_a_, 0,
      0, std_yawdd_*std_yawdd_;
  P_aug.bottomRightCorner(2, 2) = Q;

  //create square root matrix
  MatrixXd A_aug = P_aug.llt().matrixL();

  //create augmented sigma points
  Xsig_aug.col(0) = x_aug;
  for (int i = 0; i < n_aug_; i++) {
    Xsig_aug.col(i+1) = x_aug + (sqrt(lambda_ + n_aug_) * A_aug.col(i));
    Xsig_aug.col(i+1+n_aug_) = x_aug - (sqrt(lambda_ + n_aug_) * A_aug.col(i));
  }

  //predict sigma points
  for (int i = 0; i < Xsig_aug.cols(); i++) {
    double p_x = Xsig_aug(0, i);
    double p_y = Xsig_aug(1, i);
    double v = Xsig_aug(2, i);
    double psi = Xsig_aug(3, i);
    double psi_dot = Xsig_aug(4, i);
    double nu_a = Xsig_aug(5, i);
    double nu_psi_dot_dot = Xsig_aug(6, i);

    double dt_2 = delta_t / 2.0;
    double dt_nu_a = delta_t * nu_a;
    double dt2_2_nu_a = dt_2 * dt_nu_a;
    double dt_nu_psi_dot_dot = delta_t * nu_psi_dot_dot;

    //write predicted sigma points into right column
    //avoid division by zero
    if (fabs(psi_dot) > 0.001) {
      double v_psi_dot = v / psi_dot;
      double psi_psidot_dt = psi + psi_dot * delta_t;
      Xsig_pred_(0, i) = p_x + v_psi_dot * (sin(psi_psidot_dt) - sin(psi)) + dt2_2_nu_a * cos(psi);
      Xsig_pred_(1, i) = p_y + v_psi_dot * (-cos(psi_psidot_dt) + cos(psi)) + dt2_2_nu_a * sin(psi);

    } else {
      double v_dt = v * delta_t;
      Xsig_pred_(0, i) = p_x + v_dt * cos(psi) + dt2_2_nu_a * cos(psi);
      Xsig_pred_(1, i) = p_y + v_dt * sin(psi) + dt2_2_nu_a * sin(psi);
    }
    Xsig_pred_(2, i) = v + dt_nu_a;
    Xsig_pred_(3, i) = psi + psi_dot * delta_t + dt_2 * dt_nu_psi_dot_dot;
    Xsig_pred_(4, i) = psi_dot + dt_nu_psi_dot_dot;

  }

  //predict state mean
  x_.fill(0.0);
  for (int i = 0; i < 2*n_aug_+1; i++) {
    x_ += weights_(i) * Xsig_pred_.col(i);
  }

  //predict state covariance matrix
  P_.fill(0.0);
  for (int i = 0; i < 2*n_aug_+1; i++) {
    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    //normalize angle
    x_diff(3) = atan2(sin(x_diff(3)), cos(x_diff(3)));
    P_ += weights_(i) * (x_diff * x_diff.transpose());
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

  /**
  * update the state by using Kalman Filter equations
  */
  //measurement covariance matrix - laser
  MatrixXd R = MatrixXd(2, 2);
  R << std_laspx_, 0,
      0, std_laspy_;

  MatrixXd H = MatrixXd(2, 5);
  // measurement matrix - laser
  H << 1, 0, 0, 0, 0,
      0, 1, 0, 0, 0;

  VectorXd z_pred = H * x_;
  VectorXd y = meas_package.raw_measurements_ - z_pred;
  MatrixXd Ht = H.transpose();
  MatrixXd S = H * P_ * Ht + R;
  MatrixXd Si = S.inverse();
  MatrixXd PHt = P_ * Ht;
  MatrixXd K = PHt * Si;

  //new estimate
  x_ = x_ + (K * y);
  long x_size = x_.size();
  MatrixXd I = MatrixXd::Identity(x_size, x_size);
  P_ = (I - K * H) * P_;

  double NIS = y.transpose() * Si * y;
  cout << "Lidar NIS: " << NIS << endl;
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
  //set measurement dimension, radar can measure r, phi, and r_dot
  int n_z = 3;

  //create vector for incoming radar measurement
  VectorXd z = VectorXd(n_z);
  z << meas_package.raw_measurements_[0],
      meas_package.raw_measurements_[1],
      meas_package.raw_measurements_[2];

  //create matrix for sigma points in measurement space
  MatrixXd Zsig = MatrixXd(n_z, 2 * n_aug_ + 1);

  //mean predicted measurement
  VectorXd z_pred = VectorXd(n_z);

  //measurement covariance matrix S
  MatrixXd S = MatrixXd(n_z,n_z);

  //transform sigma points into measurement space
  for (int i = 0; i < 2*n_aug_+1; i++) {
    double p_x = Xsig_pred_(0,i);
    double p_y = Xsig_pred_(1,i);
    double v = Xsig_pred_(2,i);
    double psi = Xsig_pred_(3,i);
    double rho = sqrt(p_x*p_x + p_y*p_y);
    Zsig(0,i) = rho;
    Zsig(1,i) = atan2(p_y, p_x);
    Zsig(2,i) = (p_x*cos(psi)*v + p_y*sin(psi)*v) / rho;
  }

  //calculate mean predicted measurement
  z_pred.fill(0.0);
  for (int i = 0; i < 2*n_aug_+1; i++) {
    z_pred += weights_(i) * Zsig.col(i);
  }

  //calculate measurement covariance matrix S
  S.fill(0.0);
  for (int i = 0; i < 2*n_aug_+1; i++) {
    VectorXd z_diff = Zsig.col(i) - z_pred;
    //normalize angle
    z_diff(1) = atan2(sin(z_diff(1)), cos(z_diff(1)));
    S += weights_(i) * z_diff * z_diff.transpose();
  }
  MatrixXd R = MatrixXd(3, 3);
  R << std_radr_*std_radr_, 0, 0,
      0, std_radphi_*std_radphi_, 0,
      0, 0, std_radrd_*std_radrd_;
  S += R;

  //create matrix for cross correlation Tc
  MatrixXd Tc = MatrixXd(n_x_, n_z);
  //calculate cross correlation matrix
  Tc.fill(0.0);
  for (int i = 0; i < 2*n_aug_+1; i++) {
    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    VectorXd z_diff = Zsig.col(i) - z_pred;
    //normalize angles
    x_diff(3) = atan2(sin(x_diff(3)), cos(x_diff(3)));
    z_diff(1) = atan2(sin(z_diff(1)), cos(z_diff(1)));
    Tc += weights_(i) * x_diff * z_diff.transpose();
  }

  //calculate Kalman gain K;
  MatrixXd Si = S.inverse();
  MatrixXd K = Tc * Si;

  //update state mean and covariance matrix
  VectorXd z_diff = z - z_pred;
  //normalize angle
  z_diff(1) = atan2(sin(z_diff(1)), cos(z_diff(1)));

  x_ += K * z_diff;
  P_ -= K * S * K.transpose();

  double NIS = z_diff.transpose() * Si * z_diff;
  cout << "Radar NIS: " << NIS << endl;
}
