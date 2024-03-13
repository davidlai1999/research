#include <iostream>
#include <fstream>
#include <string>

#include "ns3/core-module.h"
#include "ns3/network-module.h"
#include "ns3/mobility-module.h"
#include "AP_Node.h"

using namespace ns3;
AP_Node::AP_Node(int id,Vector in_pos, Ptr<Node> _node){
    Node_ID = id;
    pos = in_pos;
    residual_power = VLC_Max_Power;
    node = _node;
}

void AP_Node::Set_AP_ID(int id){
     Node_ID = id;
}

int AP_Node::GetID(void){
    return Node_ID;
}

void AP_Node::SetPosition(Vector pos_from_model){
    pos = pos_from_model;
}

Vector AP_Node::GetPosition(void){
    return pos;
}

void AP_Node::Add_Associated_UE(int uid, double power){
    accociated_UEs.push_back({uid, power});
}

void AP_Node::Insert_Associated_UE(int uid, int dst, double power){
    accociated_UEs.insert(accociated_UEs.begin()+dst, {uid, power});
}

int AP_Node::Associated_UE_Num() {
    return accociated_UEs.size();
}

void AP_Node::Remove_Associated_UE(int uid){
    for (int i = 0 ; i < accociated_UEs.size() ; i++) {
        if (accociated_UEs[i].first == uid)
            accociated_UEs.erase(accociated_UEs.begin()+i);
    }
}

void AP_Node::Clear_Associated_UE() {
    ;//accociated_UEs.clear();
}

std::vector<int> AP_Node::Get_UEs(void){
    std::vector<int> UEs;
    for (int i = 0 ; i < accociated_UEs.size() ; i++) {
        UEs.push_back(accociated_UEs[i].first);
    }

    return UEs;
}

Ptr<Node> AP_Node::Get_Node(void){
    return node;
}

bool AP_Node::Find_UE(int uid){
    for (int i = 0 ; i < accociated_UEs.size() ; i++) {
        if (accociated_UEs[i].first == uid)
            return true;
    }

    return false;
}

void AP_Node::Set_UE_Power(int uid, double required_power){
    for (int i = 0 ; i < accociated_UEs.size() ; i++) {
        if (accociated_UEs[i].first == uid)
            accociated_UEs[i].second = required_power;
    }
}

double AP_Node::Get_UE_Power(int uid){
    for (int i = 0 ; i < accociated_UEs.size() ; i++) {
        if (accociated_UEs[i].first == uid)
            return accociated_UEs[uid].second;
    }
}

double AP_Node::Get_Residual_Power(void){
    return residual_power;
}

void AP_Node::Set_Residual_Power(double Power){
    residual_power = Power;
}
