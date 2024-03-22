#ifndef UE_NODE_H
#define UE_NODE_H

#include <iostream>
#include <fstream>
#include <string>
#include <vector>

#include "ns3/core-module.h"
#include "ns3/network-module.h"
#include "ns3/mobility-module.h"
#include "global_environment.h"

using namespace ns3;

class UE_Node {
    public :
        UE_Node(int id, double required_rate, Ptr<Node> _node);

        void Set_UE_ID(int id);

        int GetID(void) const;

        void UpdateCondiction();

        void SetVelocity();

        double GetVelocity(void);

        void SetPosition();

        void SetPosition(double x, double y);

        Vector GetPosition(void);

        void Set_Required_DataRate(double data_rate_in_Mbps); // set the required rate of the UE

        double Get_Required_DataRate(void); // get the required rate of the UE

        void Set_Achievable_DataRate(double data_rate_in_Mbps); // set the achievable rate of the UE

        double Get_Achievable_DataRate(void); // get the achievable rate of the UE

        void Set_Associated_AP(int associated_AP_index); // set the associated AP of the UE

        int Get_Associated_AP(void) const; // get the associated AP of the UE

        Ptr<Node> Get_Node(void); // get the ns3 node of the UE

        void Set_Best_Channel();

        double get_Best_Channel();

        void Set_Mode_Index(int index);

        int Get_Mode_Index();

        void Set_Modulation_Mod(std::vector<int> mode);

        std::vector<int> Get_Modulation_Mod(void);

    private :
        int Node_ID;                                 //user id
        Vector ue_velocity;                          //velocity
        Vector pos;                                  //user position
        Ptr<Node> node;                              //ns3 node
        double required_datarate;                    //required data rate
        int mode_index;
        std::vector<int> modulation_mod;             //modulation mod of the UE


        double achievable_datarate;                  //achievable data rate
        int associated_AP;                           //associated AP of the UE
};
#endif

