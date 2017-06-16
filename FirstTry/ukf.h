#ifndef UKF_H
#define UKF_H

#include "measurement_package.h"
#include "Eigen/Dense"
#include <vector>
#include <string>
#include <fstream>

using Eigen::MatrixXd;
using Eigen::VectorXd;

class UKF {
public:

  ///* initially set to false, set to true in first call of ProcessMeasurement
  bool is_initialized_;

  ///* if this is false, laser measurements will be ignored (except for init)
  bool use_laser_;

  ///* if this is false, radar measurements will be ignored (except for init)
  bool use_radar_;

  ///* state vector: [pos1 pos2 vel_abs yaw_angle yaw_rate] in SI units and rad
  VectorXd x_;

  ///* state covariance matrix
  MatrixXd P_;

  ///* predicted sigma points matrix
  MatrixXd Xsig_pred_;

  ///* time when the state is true, in us
  long long time_us_;

  ///* Process noise standard deviation longitudinal acceleration in m/s^2
  double std_a_;

  ///* Process noise standard deviation yaw acceleration in rad/s^2
  double std_yawdd_;

  ///* Laser measurement noise standard deviation position1 in m
  double std_laspx_;

  ///* Laser measurement noise standard deviation position2 in m
  double std_laspy_;

  ///* Radar measurement noise standard deviation radius in m
  double std_radr_;

  ///* Radar measurement noise standard deviation angle in rad
  double std_radphi_;

  ///* Radar measurement noise standard deviation radius change in m/s
  double std_radrd_ ;

  ///* Weights of sigma points
  VectorXd weights_;

  ///* State dimension
  int n_x_;

  ///* Augmented state dimension
  int n_aug_;
  
  ///* Sigma point spreading parameter
  double lambda_;
  
  ///************************MY ADAPTION********************************
  
  ///* dimension of lidar measurement
  int dim_laser_measurement_;
  
  ///* dimension of radar measurement
  int dim_radar_measurement_;
  
  ///* the current NIS (Normalized Innovation Squared) for Lidar
  double NIS_laser_;
  
  ///* the current NIS (Normalized Innovation Squared) for Radar
  double NIS_radar_;
  
  ///* Linear measurement matrix for lidar updata (predict measurement)
  MatrixXd H_laser_;
  
  ///* Lidar measurement covariance matrix
  MatrixXd R_laser_;
  
  ///* Radar measurement covariance matrix
  MatrixXd R_radar_;
  
  ///*******************************************************************

  /**
   * Constructor
   */
  UKF();

  /**
   * Destructor
   */
  virtual ~UKF();

  /**
   * ProcessMeasurement
   * @param meas_package The latest measurement data of either radar or laser
   */
  void ProcessMeasurement(MeasurementPackage meas_package);

  /**
   * Prediction Predicts sigma points, the state, and the state covariance
   * matrix
   * @param delta_t Time between k and k+1 in s
   */
  MatrixXd Prediction(double delta_t);

  /**
   * Updates the state and the state covariance matrix using a laser measurement
   * @param meas_package The measurement at k+1
   */
  void UpdateLidar(MeasurementPackage meas_package);

  /**
   * Updates the state and the state covariance matrix using a radar measurement
   * @param meas_package The measurement at k+1
   */
  void UpdateRadar(MeasurementPackage meas_package);
  
  ///* A help function to generate augmented sigma points, used in the prediction step.
  void AugmentedSigmaPoints(MatrixXd* Xsig_out);
  
  ///* A help function to generate predict sigma points, used in the prediction step.
  void SigmaPointPrediction(MatrixXd* Xsig_out);
  
  ///* A help function to calculate mean and covariance of a given matrix X_in, where columns are variables.
  void PredictMeanAndCovariance(VectorXd* x_out, MatrixXd* P_out, MatrixXd X_in);
  
  ///* A help function to calculate the Cross Correlation Matrix, given predict sigma points matrix Xsig_in, mean vector of predict state x_in, measurement sigma points matrix Zsig_in, mean vector of measurement state z_in and dimension of measurement space dim.
  MatrixXd GetCrossCorrelationMatrix(MatrixXd Xsig_in, VectorXd x_in, MatrixXd Zsig_in, VectorXd z_in, int dim);
  
  ///* A help function to update state when correlation matrix Tc, measurement covariance matrix S, raw measurements and mean predict measurement are ready
  void UpdateState(MatrixXd Tc, MatrixXd S, VectorXd raw_measurements, VectorXd z_pred, int dim);
  
  ///* A help function to transform the predict sigma points to measurement space.
  void TransformToRadarMeasurementSpace(MatrixXd* Z_out);
  
  
};

#endif /* UKF_H */
