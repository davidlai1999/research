#ifndef CHANNEL_H
#define CHANNEL_H

#include <iostream>
#include <fstream>
#include <string>

#include "ns3/core-module.h"
#include "ns3/network-module.h"
#include "ns3/mobility-module.h"
#include "global_environment.h"
#include "AP_Node.h"
#include "UE_Node.h"

using namespace ns3;

double RadToDeg(const double & radian);

double DegToRad(const double & degree);

double GetDistance_AP_UE(Ptr<Node> AP, Ptr<Node> UE);

double GetDistance_AP_UE(AP_Node ap, UE_Node ue); // Debug

double Get_Incidence_Angle_AP_UE(Ptr<Node> AP, Ptr<Node> UE);

double Get_Incidence_Angle_AP_UE(AP_Node ap, UE_Node ue); // Debug

void Calculate_Channel_Gain_Matrix(std::vector<AP_Node> & myAPlist, std::vector<UE_Node> & myUElist, std::vector<std::vector<double>> & Channel_Gain_Matrix);

double Estimate_one_RF_Channel_Gain(AP_Node ap, UE_Node ue);

double Estimate_one_VLC_Channel_Gain(AP_Node ap, UE_Node ue);

double Estimate_one_VLC_Channel_Gain2(AP_Node ap, UE_Node ue);

#endif
