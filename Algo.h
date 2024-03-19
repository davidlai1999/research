#ifndef ALGO_H
#define ALGO_H

#include <iostream>
#include <fstream>
#include <string>

#include "ns3/core-module.h"
#include "ns3/network-module.h"
#include "ns3/mobility-module.h"
#include "Channel.h"

#include "global_environment.h"
#include "UE_Node.h"
#include "AP_Node.h"
#include "ns3/point-to-point-module.h"

using namespace ns3;

void Initialize_UE_Node_list(NodeContainer &UE_Nodes);

void Initialize_AP_Node_list(NodeContainer & AP_Nodes);

void Update_UE_Condiction(std::vector<UE_Node> & UElist);

double set_BER_Required_SINR(const int& number_of_layer);

double get_layer_rate(int layer, double SINR);

double get_achievable_data_rate(std::vector<int> mode);

void create_modulation_mode_rate_table();

double get_layer_required_power(const int ap_id, const int ue_id, const int mod, const int layer, const double sum_prev_ue_power);

double get_VLC_minimum_required_power(const int ap_id, const int ue_id, const double sum_prev_ue_power);

double get_RF_minimum_required_power(const int ap_id, const int ue_id, const double sum_prev_ue_power);

std::vector<double> power_allocation_after_handover(int ap_id, int ue_id);

int find_Best_VLC_Channel_Gain(int ue_id);

void set_UE_Associated_AP_Based_On_Channel();

void Update_VLC_AP_power_allocation(int ap_id);

void Update_RF_AP_power_allocation(int ap_id);

void Update_All_AP_power_allocation();

double generate_VHO_overhead();

double generate_HHO_overhead();

void Do_algorithm();

#endif
