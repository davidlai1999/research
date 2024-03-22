/* -*- Mode:C++; c-file-style:"gnu"; indent-tabs-mode:nil; -*- */
/*
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License version 2 as
 * published by the Free Software Foundation;
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 *
 */

//
// Network topology
//
//           10Mb/s, 10ms       10Mb/s, 10ms
//       n0-----------------n1-----------------n2
//
//
// - Tracing of queues and packet receptions to file
//   "tcp-large-transfer.tr"
// - pcap traces also generated in the following files
//   "tcp-large-transfer-$n-$i.pcap" where n and i represent node and interface
// numbers respectively
//  Usage (e.g.): ./waf --run tcp-large-transfer

#include <iostream>
#include <fstream>
#include <string>

#include "ns3/core-module.h"
#include "ns3/applications-module.h"
#include "ns3/network-module.h"
#include "ns3/internet-module.h"
#include "ns3/point-to-point-module.h"
#include "ns3/ipv4-global-routing-helper.h"
#include "global_environment.h"
#include "install_mobility.h"
#include "UE_Node.h"
#include "Algo.h"
#include "Channel.h"

using namespace ns3;

NS_LOG_COMPONENT_DEFINE ("VisibleLightCommunication");

int main (int argc, char *argv[])
{
  /** Create AP Nodes, install mobility model and initialize the AP list for the algorithm**/
  NodeContainer AP_Nodes;

  AP_Nodes.Create(RF_AP_Num+VLC_AP_Num);

  install_AP_mobility(AP_Nodes);

  Initialize_AP_Node_list(AP_Nodes);

  // Create UE Nodes , install mobility model and initialize the UE list for the algorithm
  NodeContainer UE_Nodes;

  UE_Nodes.Create (UE_Num);

  install_UE_mobility(UE_Nodes);

  Initialize_UE_Node_list(UE_Nodes);

  set_BER_Required_SINR(number_of_layer);

  create_modulation_mode_rate_table();

  //schedule the algorithm
  //Do_algorithm();
  Simulator::Schedule(Seconds(0.0),&Do_algorithm);

  //set the simulation time
  Simulator::Stop (Seconds (simulation_time));

  //start simulation
  Simulator::Run ();

  //Output();
  Simulator::Destroy ();

  print_analysis_data();
}
