
#include <vector>


template <class Derived>
class Neurons {
    public:
        void Connect(Neurons<Derived> &other) {
            auto me = static_cast<Derived*>(this);
            auto he = static_cast<Derived*>(other);
            for(auto member: me) {
                for(auto h_member: he) {
                    for(auto it = member->begin() ; it< member->end() ; ++it) {
                        for(auto it2 = h_member->begin() ; it2 < h_member->end(); ++it2) {
                            it->out.push_back(it2);
                            it2->in.push_back(it);
                        }
                    }
                }
            }
        }
};




class Neuron: public Neurons<Neuron> {
    public:
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


class NeuronLayer: public std::vector<Neuron>, Neurons<NeuronLayer> {
    public:
        NeuronLayer(int number_of_neurons) {
            while (number_of_neurons-- >0) {
            emplace_back(Neuron{});
        }
        }
};

int main() {
    NeuronLayer layer1(5), layer2(3);
    static_cast<Neurons<NeuronLayer>>(layer1).Connect(static_cast<Neurons<NeuronLayer>>(layer2));
    
    return 0;
}