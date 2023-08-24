/******************************************************************************
 * Copyright 2019 The Apollo Authors. All Rights Reserved.
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 * http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 *****************************************************************************/

/**
 * @file
 **/

#pragma once

#include <memory>
#include <vector>

#include "Eigen/Eigen"

#ifdef ALIVE
#undef ALIVE
#endif

#include "modules/common_msgs/config_msgs/vehicle_config.pb.h"
#include "modules/common/math/vec2d.h"
#include "modules/common/vehicle_state/proto/vehicle_state.pb.h"
#include "modules/planning/common/trajectory/discretized_trajectory.h"
#include "modules/planning/open_space/coarse_trajectory_generator/hybrid_a_star.h"
#include "modules/planning/open_space/trajectory_smoother/distance_approach_problem.h"
#include "modules/planning/open_space/trajectory_smoother/dual_variable_warm_start_problem.h"
#include "modules/planning/open_space/trajectory_smoother/iterative_anchoring_smoother.h"
#include "modules/planning/proto/open_space_task_config.pb.h"
/***
 * HybridAstar search
 * path smooth, trajectory plan
 * 
*/
namespace apollo {
namespace planning {
class OpenSpaceTrajectoryOptimizer {
 public:
  OpenSpaceTrajectoryOptimizer(
      const OpenSpaceTrajectoryOptimizerConfig& config);

  virtual ~OpenSpaceTrajectoryOptimizer() = default;

  common::Status Plan(
      const std::vector<common::TrajectoryPoint>& stitching_trajectory,
      const std::vector<double>& end_pose, const std::vector<double>& XYbounds,
      double rotate_angle, const common::math::Vec2d& translate_origin,
      const Eigen::MatrixXi& obstacles_edges_num,
      const Eigen::MatrixXd& obstacles_A, const Eigen::MatrixXd& obstacles_b,
      const std::vector<std::vector<common::math::Vec2d>>&
          obstacles_vertices_vec,
      double* time_latency);

  void GetStitchingTrajectory(
      std::vector<common::TrajectoryPoint>* stitching_trajectory) {
    stitching_trajectory->clear();
    *stitching_trajectory = stitching_trajectory_;
  }
  // 输出
  void GetOptimizedTrajectory(DiscretizedTrajectory* optimized_trajectory) {
    optimized_trajectory->clear();
    *optimized_trajectory = optimized_trajectory_;
  }

  void RecordDebugInfo(
      const common::TrajectoryPoint& trajectory_stitching_point,
      const common::math::Vec2d& translate_origin, const double rotate_angle,
      const std::vector<double>& end_pose, const Eigen::MatrixXd& xWS,
      const Eigen::MatrixXd& uWs, const Eigen::MatrixXd& l_warm_up,
      const Eigen::MatrixXd& n_warm_up, const Eigen::MatrixXd& dual_l_result_ds,
      const Eigen::MatrixXd& dual_n_result_ds,
      const Eigen::MatrixXd& state_result_ds,
      const Eigen::MatrixXd& control_result_ds,
      const Eigen::MatrixXd& time_result_ds,
      const std::vector<double>& XYbounds,
      const std::vector<std::vector<common::math::Vec2d>>&
          obstacles_vertices_vec);

  void UpdateDebugInfo(
      ::apollo::planning_internal::OpenSpaceDebug* open_space_debug);

  apollo::planning_internal::OpenSpaceDebug* mutable_open_space_debug() {
    return &open_space_debug_;
  }

 private:
  bool IsInitPointNearDestination(
      const common::TrajectoryPoint& planning_init_point,
      const std::vector<double>& end_pose, double rotate_angle,
      const common::math::Vec2d& translate_origin);

  void PathPointNormalizing(double rotate_angle,
                            const common::math::Vec2d& translate_origin,
                            double* x, double* y, double* phi);

  void PathPointDeNormalizing(double rotate_angle,
                              const common::math::Vec2d& translate_origin,
                              double* x, double* y, double* phi);

  void LoadTrajectory(const Eigen::MatrixXd& state_result_ds,
                      const Eigen::MatrixXd& control_result_ds,
                      const Eigen::MatrixXd& time_result_ds);

  void LoadHybridAstarResultInEigen(HybridAStartResult* result,
                                    Eigen::MatrixXd* xWS, Eigen::MatrixXd* uWS);

  void UseWarmStartAsResult(
      const Eigen::MatrixXd& xWS, const Eigen::MatrixXd& uWS,
      const Eigen::MatrixXd& l_warm_up, const Eigen::MatrixXd& n_warm_up,
      Eigen::MatrixXd* state_result_ds, Eigen::MatrixXd* control_result_ds,
      Eigen::MatrixXd* time_result_ds, Eigen::MatrixXd* dual_l_result_ds,
      Eigen::MatrixXd* dual_n_result_ds);

  bool GenerateDistanceApproachTraj(
      const Eigen::MatrixXd& xWS, const Eigen::MatrixXd& uWS,
      const std::vector<double>& XYbounds,
      const Eigen::MatrixXi& obstacles_edges_num,
      const Eigen::MatrixXd& obstacles_A, const Eigen::MatrixXd& obstacles_b,
      const std::vector<std::vector<common::math::Vec2d>>&
          obstacles_vertices_vec,
      const Eigen::MatrixXd& last_time_u, const double init_v,
      Eigen::MatrixXd* state_result_ds, Eigen::MatrixXd* control_result_ds,
      Eigen::MatrixXd* time_result_ds, Eigen::MatrixXd* l_warm_up,
      Eigen::MatrixXd* n_warm_up, Eigen::MatrixXd* dual_l_result_ds,
      Eigen::MatrixXd* dual_n_result_ds);

  bool GenerateDecoupledTraj(
      const Eigen::MatrixXd& xWS, const double init_a, const double init_v,
      const std::vector<std::vector<common::math::Vec2d>>&
          obstacles_vertices_vec,
      Eigen::MatrixXd* state_result_dc, Eigen::MatrixXd* control_result_dc,
      Eigen::MatrixXd* time_result_dc);

  void LoadResult(const DiscretizedTrajectory& discretized_trajectory,
                  Eigen::MatrixXd* state_result_dc,
                  Eigen::MatrixXd* control_result_dc,
                  Eigen::MatrixXd* time_result_dc);

  void CombineTrajectories(
      const std::vector<Eigen::MatrixXd>& xWS_vec,
      const std::vector<Eigen::MatrixXd>& uWS_vec,
      const std::vector<Eigen::MatrixXd>& state_result_ds_vec,
      const std::vector<Eigen::MatrixXd>& control_result_ds_vec,
      const std::vector<Eigen::MatrixXd>& time_result_ds_vec,
      const std::vector<Eigen::MatrixXd>& l_warm_up_vec,
      const std::vector<Eigen::MatrixXd>& n_warm_up_vec,
      const std::vector<Eigen::MatrixXd>& dual_l_result_ds_vec,
      const std::vector<Eigen::MatrixXd>& dual_n_result_ds_vec,
      Eigen::MatrixXd* xWS, Eigen::MatrixXd* uWS,
      Eigen::MatrixXd* state_result_ds, Eigen::MatrixXd* control_result_ds,
      Eigen::MatrixXd* time_result_ds, Eigen::MatrixXd* l_warm_up,
      Eigen::MatrixXd* n_warm_up, Eigen::MatrixXd* dual_l_result_ds,
      Eigen::MatrixXd* dual_n_result_ds);

 private:
  OpenSpaceTrajectoryOptimizerConfig config_;

  std::unique_ptr<HybridAStar> warm_start_;
  std::unique_ptr<DistanceApproachProblem> distance_approach_;
  std::unique_ptr<DualVariableWarmStartProblem> dual_variable_warm_start_;
  std::unique_ptr<IterativeAnchoringSmoother> iterative_anchoring_smoother_;

  std::vector<common::TrajectoryPoint> stitching_trajectory_;
  DiscretizedTrajectory optimized_trajectory_;   // result 输出

  apollo::planning_internal::OpenSpaceDebug open_space_debug_;
};
}  // namespace planning
}  // namespace apollo



/*
message PlannerOpenSpaceConfig {
  // Open Space ROIConfig
  optional ROIConfig roi_config = 1;
  // Hybrid A Star Warm Start
  optional WarmStartConfig warm_start_config = 2;
  // Dual Variable Warm Start
  optional DualVariableWarmStartConfig dual_variable_warm_start_config = 3;
  // Distance Approach Configs
  optional DistanceApproachConfig distance_approach_config = 4;
  // Iterative Anchoring Configs
  optional IterativeAnchoringConfig iterative_anchoring_smoother_config = 5;
  // Trajectory PartitionConfig Configs
  optional TrajectoryPartitionConfig trajectory_partition_config = 6;
  optional float delta_t = 7 [default = 1.0];
  optional double is_near_destination_threshold = 8 [default = 0.001];
  optional bool enable_check_parallel_trajectory = 9 [default = false];
  optional bool enable_linear_interpolation = 10 [default = false];
  optional double is_near_destination_theta_threshold = 11 [default = 0.05];
}

message ROIConfig {
  // longitudinal range of parking roi backward
  optional double roi_longitudinal_range_start = 1 [default = 10.0];
  // longitudinal range of parking roi forward
  optional double roi_longitudinal_range_end = 2 [default = 10.0];
  // parking spot range detection threshold
  optional double parking_start_range = 3 [default = 7.0];
  // Parking orientation for reverse parking
  optional bool parking_inwards = 4 [default = false];
}

message WarmStartConfig {
  // Hybrid a star for warm start
  optional double xy_grid_resolution = 1 [default = 0.2];
  optional double phi_grid_resolution = 2 [default = 0.05];
  optional uint64 next_node_num = 3 [default = 10];
  optional double step_size = 4 [default = 0.5];
  optional double traj_forward_penalty = 5 [default = 0.0];
  optional double traj_back_penalty = 6 [default = 0.0];
  optional double traj_gear_switch_penalty = 7 [default = 10.0];
  optional double traj_steer_penalty = 8 [default = 100.0];
  optional double traj_steer_change_penalty = 9 [default = 10.0];
  // Grid a star for heuristic
  optional double grid_a_star_xy_resolution = 15 [default = 0.1];
  optional double node_radius = 16 [default = 0.5];
  optional PiecewiseJerkSpeedOptimizerConfig s_curve_config = 17;
  optional double traj_kappa_contraint_ratio = 10 [default = 0.7];
}

message IterativeAnchoringConfig {
  // Ipopt configs
  optional double interpolated_delta_s = 1 [default = 0.1];
  optional int32 reanchoring_trails_num = 2 [default = 50];
  optional double reanchoring_pos_stddev = 3 [default = 0.25];
  optional double reanchoring_length_stddev = 4 [default = 1.0];
  optional bool estimate_bound = 5 [default = false];
  optional double default_bound = 6 [default = 2.0];
  optional double vehicle_shortest_dimension = 7 [default = 1.04];
  optional FemPosDeviationSmootherConfig fem_pos_deviation_smoother_config = 8;
  optional double collision_decrease_ratio = 9 [default = 0.9];
  // TODO(QiL, Jinyun): Merge following with overall config for open space
  optional double max_forward_v = 10 [default = 2.0];
  optional double max_reverse_v = 11 [default = 2.0];
  optional double max_forward_acc = 12 [default = 3.0];
  optional double max_reverse_acc = 13 [default = 2.0];
  optional double max_acc_jerk = 14 [default = 4.0];
  optional double delta_t = 15 [default = 0.2];
  optional PiecewiseJerkSpeedOptimizerConfig s_curve_config = 16;
}

/*