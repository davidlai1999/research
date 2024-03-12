#ifndef INSTALL_MOBILITY_H
#define INSTALL_MOBILITY_H
#include <iostream>
#include <fstream>
#include <string>

#include "ns3/core-module.h"
#include "ns3/network-module.h"
#include "ns3/mobility-module.h"


using namespace ns3;

void install_AP_mobility(NodeContainer &AP_Nodes);
void install_UE_mobility(NodeContainer &UE_Nodes);

void print_AP_position(NodeContainer &AP_Nodes);
void print_UE_position(NodeContainer &UE_Nodes);
#endif

