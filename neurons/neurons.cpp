#include <vector>
#include <iostream>
struct NeuronLayer;
//template<class Derived>
template<class Derived>
struct Neurons {
    void connect(Derived &other);
    //std::vector<struct NeuronLayer> layersOut;
};

struct Neuron: Neurons <Neuron>{
    std::vector<Neuron*> in;
    std::vector<Neuron*> out;
    unsigned int ID;
    Neuron() {
        static int ID = 1;
        this->ID = ID++;
    }    
    Neuron* begin() {return this;}
    Neuron* end() {return this+1;}
};

struct NeuronLayer: std::vector<Neuron>, Neurons<struct NeuronLayer> {
    NeuronLayer(int number_of_neurons) {
        while (number_of_neurons-- >0) {
            emplace_back(Neuron{});
        }
    }
};
template<class Derived>
void Neurons<Derived>::connect(Derived &other) {
    //auto me = static_cast<Derived*>(this);
    //auto he = static_cast<Derived*>(&other);
    //for(auto it = me->begin(); it< me->end(); ++it) {
    //    for(auto it2 = he->begin(); it2< he->end(); ++it2) {
    //        (it)->out.push_back(*it2);
    //        (it2)->in.push_back(*it);
    //    }
    //}
}

std::ostream &operator<<(std::ostream &console, Neurons<Neuron> &neurons) {
    //???
    //auto me = static_cast<struct NeuronLayer&>(neurons);
    //for(auto it = me.begin(); it< me.end(); ++it) {
    //    for (Neuron *n : (*it).in)
    //        console << n->ID << "\t>\t" << (*it).ID << "*" << "\n";
    //    for (Neuron *n : (*it).out)
    //        console << (*it).ID << "*\t>\t" << n->ID << "\n";
    //}
    return console;
}

int main() {
   // NeuronLayer layer_1{1}, layer_2{2};
   // layer_1.connect(layer_2);
    return 0;
}