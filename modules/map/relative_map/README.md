# Relative Map

## Introduction
Relative map module is a middle layer connecting HDMap/Perception module and planning module. This module generates real-time relative map in the body coordinate system (FLU) and a reference line for planning. The inputs for relative map module have two parts: offline and online. The offline part is navigation line (human driving path) and the HDMap information on and near navigation line. And the online part is the traffic sign related information from perception module, e.g., lane marker, crosswalk, traffic light etc. The generation of relative map can leverage both online and offline parts. It also works with either online or offline part only.

- 基于ego 坐标系
- relative_map 包含两部分数据
    离线部分： 人类伺机驾驶路线 navigation line, HDMap info(左中右车道)
    在线部分： 感知：交通灯、车道线， 人行道

## Inputs
  * NavigationInfo from dreamview module
  * LaneMarker from perception module
  * Localization from localization module

## Output
  * Relative map follows map format defined in modules/map/proto/map.proto

  common_msgs/map_msgs/map.proto
