#include <iostream>
#include <fstream>
#include <string>

#include "ns3/core-module.h"
#include "ns3/network-module.h"
#include "ns3/mobility-module.h"
#include "UE_Node.h"

using namespace ns3;
UE_Node::UE_Node(int id, double required_rate, Ptr<Node> _node) {
    Node_ID = id;
    required_datarate = required_rate;
    associated_AP = -1;
    node = _node;
    mode_index = -1;
    modulation_mod = {0, 0, 0, 0};
}

void UE_Node::Set_UE_ID(int id){
     Node_ID = id;
}

int UE_Node::GetID(void) const{
    return Node_ID;
}

void UE_Node::SetVelocity() {
    Ptr<MobilityModel> VLC_UE_MobilityModel = node->GetObject<MobilityModel>();
    ue_velocity = VLC_UE_MobilityModel->GetVelocity();
}

double UE_Node::GetVelocity(void) {
    double result = sqrt(pow(ue_velocity.x, 2) + pow(ue_velocity.y, 2));
    return result;
}

void UE_Node::SetPosition(){
    Ptr<MobilityModel> VLC_UE_MobilityModel = node->GetObject<MobilityModel>();
    pos = VLC_UE_MobilityModel->GetPosition();
}

void UE_Node::SetPosition(double x, double y){
    pos.x = x;
    pos.y = y;
}

Vector UE_Node::GetPosition(void){
    return pos;
}

void UE_Node::UpdateCondiction() {
    SetVelocity();
    SetPosition();
}

void UE_Node::Set_Required_DataRate(double data_rate_in_Mbps){
    required_datarate = data_rate_in_Mbps;
}

double UE_Node::Get_Required_DataRate(void){
    return required_datarate;
}

void UE_Node::Set_Achievable_DataRate(double data_rate_in_Mbps) {
    achievable_datarate = data_rate_in_Mbps;
}

double UE_Node::Get_Achievable_DataRate(void) {
    return achievable_datarate;
}

void UE_Node::Set_Associated_AP(int associated_AP_index){
    associated_AP = associated_AP_index;
}

int UE_Node::Get_Associated_AP(void) const{
    return associated_AP;
}

Ptr<Node> UE_Node::Get_Node(void){
    return node;
}

void UE_Node::Set_Mode_Index(int index) {
    mode_index = index;
}

int UE_Node::Get_Mode_Index() {
    return mode_index;
}


void UE_Node::Set_Modulation_Mod(std::vector<int> mode){
    modulation_mod = mode;
}

std::vector<int> UE_Node::Get_Modulation_Mod(void){
    return modulation_mod;
}

