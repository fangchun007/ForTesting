#include "ukf.h"
#include "tools.h"
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
    
    ///* initially set to false, set to true in first call of ProcessMeasurement
    is_initialized_ = false;
  
    ///* time when state is true, in us / previous time stamp
    time_us_ = 0;
  
    // if this is false, laser measurements will be ignored (except during init)
    use_laser_ = true;
    
    // if this is false, radar measurements will be ignored (except during init)
    use_radar_ = true;
    
    ///* State dimension
    n_x_ = 5;
    
    // initial state vector
    x_ = VectorXd(n_x_);
    
    // initial covariance matrix
    P_ = MatrixXd(n_x_, n_x_);

    ///* Augmented state dimension
    n_aug_ = n_x_+2;
    
    ///* Sigma point spreading parameter
    lambda_ = 3-n_x_;
    
    ///* Weights of sigma points
    weights_ = VectorXd(2*n_aug_+1);
    // set weights
    double weight_0 = lambda_/(lambda_+n_aug_);
    weights_(0) = weight_0;
    for (int i=1; i<2*n_aug_+1; i++) {
      double weight = 0.5/(n_aug_+lambda_);
      weights_(i) = weight;
    }
    
    ///* the current NIS for radar
    NIS_radar_ = 0.0;
    
    ///* the current NIS for laser
    NIS_laser_ = 0.0;
    
    // Process noise standard deviation longitudinal acceleration in m/s^2
    std_a_ = 1;
    
    // Process noise standard deviation yaw acceleration in rad/s^2
    std_yawdd_ = 0.5;
    
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
    
  
    
    ///* linear measurement matrix for lidar update
    H_laser_ = MatrixXd(2,5);
    H_laser_ << 1, 0, 0, 0, 0,
                0, 1, 0, 0, 0;
  
    ///* dimension of measurement space for radar
    dim_radar_measurement_ = 3;
  
    ///* dimension of measurement space for lidar
    dim_laser_measurement_ = 2;
  
    ///* lidar measurement covariance matrix
    R_laser_ = MatrixXd(2,2);
    R_laser_ << std_laspx_ * std_laspx_, 0,
                0, std_laspy_ * std_laspy_;
  
    ///* Radar measurement covariance matrix
    R_radar_ = MatrixXd(dim_radar_measurement_, dim_radar_measurement_);
    R_radar_ << std_radr_ * std_radr_, 0,                         0,
                0,                     std_radphi_ * std_radphi_, 0,
                0,                     0,                         std_radrd_ * std_radrd_;
  
    //create matrix with augmented sigma points as columns
    Xsig_aug_ = MatrixXd(n_aug_, 2*n_aug_+1);

  
    ///* sigma points in matrix form
    Xsig_pred_ = MatrixXd(n_x_, 2*n_aug_+1);
    
    /**
     TODO:
     Complete the initialization. See ukf.h for other member properties.
     Hint: one or more values initialized above might be wildly off...
     Yes, on the first glance, at leat std_a_ and std_yawdd_ are too big for a bicycle.
     */
    
    // See the above codes which start with "///*"
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
  
    /***********************************************************************
     *                        Initialization                               *
     ***********************************************************************/
  
    if (!is_initialized_){
      
        cout << "Unscented Kalman Filter Intialization: " << endl;
        if (meas_package.sensor_type_ == MeasurementPackage::RADAR){
          x_ = tools.ConvertPolarToCartesian(meas_package.raw_measurements_);
        }
        else if (meas_package.sensor_type_ == MeasurementPackage::LASER) {
            x_ << meas_package.raw_measurements_[0], meas_package.raw_measurements_[1], 0, 0, 0;
        }
        P_ << 10, 0, 0, 0, 0,
              0, 10, 0, 0, 0,
              0, 0, 10, 0, 0,
              0, 0, 0, 10, 0,
              0, 0, 0, 0, 10;
      
        // set time stamp value
        time_us_ = meas_package.timestamp_;
      
        is_initialized_ = true;
      
        cout << "UKF Initialization Completed!" << endl;
      
        return;
    }
  
  
    /***********************************************************************
     *                        Prediction                                   *
     ***********************************************************************/
  
    float dt = (meas_package.timestamp_ - time_us_) / 1000000.0;	//dt - expressed in seconds
    time_us_ = meas_package.timestamp_;
    
    Prediction(dt);
  
  
    /***********************************************************************
     *                        Update                                       *
     ***********************************************************************/
  
    // Here, variables use_radar_ and use_laser_ are useful for calculation of NISes
  
    if (meas_package.sensor_type_ == MeasurementPackage::RADAR)
    {
      if (use_radar_) {
        UpdateRadar(meas_package);
      }
    }
    else if (meas_package.sensor_type_ == MeasurementPackage::LASER)
    {
      if (use_laser_) {
        UpdateLidar(meas_package);

      }
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
  
    // Generate Augmented Sigma Points
    Xsig_aug_.fill(0.0);
    Xsig_aug_ = AugmentedSigmaPoints(x_, P_);
  
    // Generate Predict Sigma Points
    Xsig_pred_.fill(0.0);
    Xsig_pred_ = SigmaPointPrediction(Xsig_aug_, delta_t);
  
    // Predict Mean and Covariance
    PredictMeanAndCovariance(&x_, &P_, Xsig_pred_, n_x_);

}


MatrixXd UKF::AugmentedSigmaPoints(VectorXd x, MatrixXd P) {
  
  //create augmented mean vector
  VectorXd x_aug = VectorXd(n_aug_);
  //create augmented state covariance
  MatrixXd P_aug = MatrixXd(n_aug_, n_aug_);
  //create augmented sigma point matrix
  MatrixXd Xsig_aug = MatrixXd(n_aug_, 2*n_aug_+1);
  
  //set augmented mean state
  x_aug.head(5) = x;
  x_aug(5)      = 0;
  x_aug(6)      = 0;
  //set augmented covariance matrix
  P_aug.fill(0.0);
  P_aug.topLeftCorner(n_x_, n_x_) = P;
  P_aug(5,5)                      = std_a_ * std_a_;
  P_aug(6,6)                      = std_yawdd_ * std_yawdd_;
  
  //create square root matrix
  MatrixXd L = P_aug.llt().matrixL();
  
  //create augmented sigma points
  Xsig_aug.col(0) = x_aug;
  for (int i=0; i<n_aug_; i++)
  {
    Xsig_aug.col(i+1)        = x_aug + sqrt(lambda_ + n_aug_) * L.col(i);
    Xsig_aug.col(i+1+n_aug_) = x_aug - sqrt(lambda_ + n_aug_) * L.col(i);
  }
  
  //return augmented sigma points
  return Xsig_aug;
}

MatrixXd UKF::SigmaPointPrediction(MatrixXd Xsig_in, double delta_t) {
  
  //create matrix with predicted sigma points as columns
  MatrixXd Xsig_pred = MatrixXd(n_x_, 2*n_aug_+1);
  
  //predict sigma points
  for (int i=0; i<2*n_aug_+1; i++)
  {
    // extract values for better readability
    double p_x      = Xsig_in(0,i);
    double p_y      = Xsig_in(1,i);
    double v        = Xsig_in(2,i);
    double yaw      = Xsig_in(3,i);
    double yawd     = Xsig_in(4,i);
    double nu_a     = Xsig_in(5,i);
    double nu_yawdd = Xsig_in(6,i);
    
    //predict state values
    double px_p, py_p;
    
    //avoid division by zero
    if (fabs(yawd) > 0.001) {
      px_p = p_x + v/yawd * (sin(yaw + yawd*delta_t) - sin(yaw));
      py_p = p_y + v/yawd * (cos(yaw) - cos(yaw + yawd*delta_t));
    }
    else {
      px_p = p_x + v*delta_t*cos(yaw);
      py_p = p_y + v*delta_t*sin(yaw);
    }
    
    double v_p    = v;
    double yaw_p  = yaw + yawd*delta_t;
    double yawd_p = yawd;
    
    //add noise
    px_p = px_p + 0.5*nu_a*delta_t*delta_t*cos(yaw);
    py_p = py_p + 0.5*nu_a*delta_t*delta_t*sin(yaw);
    v_p  = v_p + nu_a*delta_t;
    
    yaw_p  = yaw_p + 0.5*nu_yawdd*delta_t*delta_t;
    yawd_p = yawd_p + nu_yawdd*delta_t;
    
    //write predicted sigma point into right column
    Xsig_pred(0,i) = px_p;
    Xsig_pred(1,i) = py_p;
    Xsig_pred(2,i) = v_p;
    Xsig_pred(3,i) = yaw_p;
    Xsig_pred(4,i) = yawd_p;
  }
  
  //write result
  return Xsig_pred;
}

void UKF::PredictMeanAndCovariance(VectorXd* x_out, MatrixXd* P_out, MatrixXd X_in, int dim) {
  
  //create vector for predicted state
  VectorXd x = VectorXd(dim);
  //create covariance matrix for prediction
  MatrixXd P = MatrixXd(dim, dim);
  
  //predicted state mean
  x.fill(0.0);
  for (int i=0; i<2*n_aug_+1; i++) {
    //iterate over sigma points
    x = x + weights_(i) * X_in.col(i);
  }
  
  //predicted state covariance matrix
  P.fill(0.0);
  for (int i=0; i<2*n_aug_+1; i++) {
    //iterate over sigma points
    
    // residual
    VectorXd x_diff = X_in.col(i) - x;
    //angle normalization
    if (dim == n_x_) {
      while (x_diff(3)> M_PI) x_diff(3) -= 2.*M_PI;
      while (x_diff(3)<-M_PI) x_diff(3) += 2.*M_PI;
    }
    else if (dim == dim_radar_measurement_) {
      while (x_diff(1)> M_PI) x_diff(1) -= 2.*M_PI;
      while (x_diff(1)<-M_PI) x_diff(1) += 2.*M_PI;
    }
    
    
    P = P + weights_(i) * x_diff * x_diff.transpose();
  }
  
  //write result
  *x_out = x;
  *P_out = P;
  
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
    //linear update
    VectorXd z = meas_package.raw_measurements_;
    VectorXd z_pred = H_laser_ * x_;
    VectorXd y = z - z_pred;
    MatrixXd Ht = H_laser_.transpose();
    MatrixXd PHt = P_ * Ht;
    MatrixXd S = H_laser_ * PHt + R_laser_;
    MatrixXd Si = S.inverse();
    MatrixXd K = PHt * Si;
    
    //state update
    x_ = x_ + (K * y);
    
    //angle normalization
    while (x_(3)> M_PI) x_(3)-=2.*M_PI;
    while (x_(3)<-M_PI) x_(3)+=2.*M_PI;
    
    long x_size = x_.size();
    MatrixXd I = MatrixXd::Identity(x_size, x_size);
    P_ = (I - K * H_laser_) * P_;
    
    //calculate NIS(Normalized Innovation Squared)
    NIS_laser_ = y.transpose() * S.inverse() * y;
    
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
    //receive the measurement value from radar
    VectorXd z = meas_package.raw_measurements_;

    //create matrix for sigma points in measurement space
    MatrixXd Zsig = MatrixXd(dim_radar_measurement_, 2*n_aug_+1);
    TransformToRadarMeasurementSpace(&Zsig);
  
    //create mean predicted measurement
    VectorXd z_pred = VectorXd(dim_radar_measurement_);
    z_pred.fill(0.0);
    //create measurement covariance matrix S
    MatrixXd S = MatrixXd(dim_radar_measurement_,dim_radar_measurement_);
    S.fill(0.0);
    //calculate z_pred and S
    PredictMeanAndCovariance(&z_pred, &S, Zsig, dim_radar_measurement_);
  
    //add measurement noise covariance matrix
    S = S + R_radar_;
   
    //create matrix for cross correlation Tc
    MatrixXd Tc = MatrixXd(n_x_, dim_radar_measurement_);
    //calculate cross correlation matrix
    Tc = GetCrossCorrelationMatrix(Xsig_pred_, x_, Zsig, z_pred, dim_radar_measurement_);
  
    NIS_radar_ = UpdateState(Tc, S, z, z_pred);
  
}

double UKF::UpdateState(MatrixXd Tc, MatrixXd S, VectorXd raw_measurements, VectorXd z_pred) {
  //Kalman gain K;
  MatrixXd K = Tc * S.inverse();
  
  //residual
  VectorXd z_diff = raw_measurements - z_pred;
  //angle normalization
  while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
  while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;
  
  //update state mean and covariance matrix
  x_ = x_ + K * z_diff;
  P_ = P_ - K*S*K.transpose();
  
  double NIS = z_diff.transpose() * S.inverse() * z_diff;
  return NIS;
  
}

void UKF::TransformToRadarMeasurementSpace(MatrixXd* Z_out) {
  //create matrix for sigma points in measurement space
  MatrixXd Zsig = MatrixXd(dim_radar_measurement_, 2*n_aug_+1);
  
  //transform sigma points into measurement space
  for (int i=0; i < 2*n_aug_+1; i++) {
    //2n+1 simga points
    // extract values for better readibility
    double p_x = Xsig_pred_(0,i);
    double p_y = Xsig_pred_(1,i);
    double v   = Xsig_pred_(2,i);
    double yaw = Xsig_pred_(3,i);
    
    double v1 = cos(yaw)*v;
    double v2 = sin(yaw)*v;
    
    // measurement model
    Zsig(0,i) = sqrt(p_x*p_x + p_y*p_y);                        //r
    Zsig(1,i) = atan2(p_y,p_x);                                 //phi
    Zsig(2,i) = (p_x*v1 + p_y*v2 ) / sqrt(p_x*p_x + p_y*p_y);   //r_dot
  }
  
  //write result
  *Z_out = Zsig;
}


MatrixXd UKF::GetCrossCorrelationMatrix(MatrixXd Xsig_in, VectorXd x_in, MatrixXd Zsig_in, VectorXd z_in, int dim) {
  //create matrix for cross correlation
  MatrixXd Tc = MatrixXd(n_x_, dim);
  Tc.fill(0.0);
  
  //calculate cross correlation matrix
  for (int i=0; i<2*n_aug_+1; i++) {
    //residual
    VectorXd z_diff = Zsig_in.col(i) - z_in;
    if (dim==dim_radar_measurement_) {
      //the measurement comes from radar
      //angle normalization
      while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
      while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;
    }
    
    // state difference
    VectorXd x_diff = Xsig_in.col(i) - x_in;
    //angle normalization
    while (x_diff(3)> M_PI) x_diff(3)-=2.*M_PI;
    while (x_diff(3)<-M_PI) x_diff(3)+=2.*M_PI;
    
    Tc = Tc + weights_(i) * x_diff * z_diff.transpose();
  }
  
  return Tc;
}


