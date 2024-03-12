#include <iostream>
#include <fstream>
#include <string>
#include <cmath>

#include "ns3/core-module.h"
#include "ns3/network-module.h"
#include "ns3/mobility-module.h"
#include "install_mobility.h"
#include "global_environment.h"

using namespace ns3;
//安裝AP的移動模型
void install_AP_mobility(NodeContainer & AP_Nodes){

    MobilityHelper AP_Mobility;
    Ptr < ListPositionAllocator > AP_Pos_list = CreateObject<ListPositionAllocator>();
    /** 設定RF_AP的初始位置 **/
    for (double y = RF_AP_Y ; y > -(room_size/2) ; y-=RF_AP_SPACE) {
        for(double x = -RF_AP_X ; x < (room_size/2) ; x+=RF_AP_SPACE) {
            AP_Pos_list->Add(Vector(x,y,VLC_AP_height));
        }
    }

    /** 設定VLC_AP的初始位置 **/
    for (double y = VLC_AP_Y ; y > -(room_size/2) ; y-=VLC_AP_SPACE) {
        for(double x = -VLC_AP_Y ; x < (room_size/2) ; x+=VLC_AP_SPACE) {
            AP_Pos_list->Add(Vector(x,y,VLC_AP_height));
        }
    }

    //用剛剛的postion list 設定 PositionAllocator
    AP_Mobility.SetPositionAllocator(AP_Pos_list);

    //AP靜止不動
    AP_Mobility.SetMobilityModel("ns3::ConstantPositionMobilityModel");
    AP_Mobility.Install(AP_Nodes);
}

//安裝UE的移動模型
void install_UE_mobility(NodeContainer & UE_Nodes){
    // random
    SeedManager::SetSeed(time(NULL));

    MobilityHelper UE_Mobility;

    // the position allocator for initial position and RWP model
    ObjectFactory pos;
    pos.SetTypeId("ns3::RandomRectanglePositionAllocator");
    std::stringstream ssPos;
    ssPos << "ns3::UniformRandomVariable[Min="<<-room_size / 2<<"|Max=" << room_size / 2 << "]";
    pos.Set("X", StringValue(ssPos.str()));
    pos.Set("Y", StringValue(ssPos.str()));
    Ptr<PositionAllocator>position_allocator = pos.Create()->GetObject<PositionAllocator>();
    UE_Mobility.SetPositionAllocator(position_allocator);

    // the speed
    std::stringstream ssSpeed;
    ssSpeed << "ns3::UniformRandomVariable[Min="<< 0.1 <<"|Max=" << avg_speed*2 << "]";
    
    // RandomWaypointMobilityModel
    // the pause time for RWP model
    std::stringstream ssPause;
    ssPause << "ns3::UniformRandomVariable[Min=0.0|Max=" << pause_time << "]";

    UE_Mobility.SetMobilityModel("ns3::RandomWaypointMobilityModel",
                                 "Speed", StringValue(ssSpeed.str()),
                                 "Pause", StringValue(ssPause.str()),
                                 "PositionAllocator", PointerValue(position_allocator));

    //install the normal RWP model for the rest of the UEs
    for(int i = ceil(UE_Num*crowded_portion) ; i < UE_Num ; i++)
        UE_Mobility.Install(UE_Nodes.Get(i));
}
//print the position of the APs(for debug)
void print_AP_position(NodeContainer &AP_Nodes){
    int AP_Index = 0;
    for (NodeContainer::Iterator it = AP_Nodes.Begin() ; it != AP_Nodes.End() ; ++it){
        Ptr<MobilityModel> AP_MobilityModel = (*it)->GetObject<MobilityModel>();
        Vector pos = AP_MobilityModel->GetPosition();
        std::cout << "Position of AP " << AP_Index ++ << " =("<<pos.x<<","<<pos.y<<","<<pos.z<<")" << std::endl;
    }
    std::cout << std::endl;
}

//print the position of the UEs(for debug)
void print_UE_position(NodeContainer &UE_Nodes){
    int UE_Index = 1;
    for (NodeContainer::Iterator it = UE_Nodes.Begin() ; it != UE_Nodes.End() ; ++it){
        Ptr<MobilityModel>UE_MobilityModel = (*it)->GetObject<MobilityModel>();
        Vector pos = UE_MobilityModel->GetPosition();
        std::cout << "Position of UE " << UE_Index ++ << " =(" << pos.x << "," << pos.y << "," << pos.z << ")" << std::endl;
    }
    std::cout<<std::endl;
}

