
#include <iostream>
#include <vector>
class NeuronLayer;
class Neuron;
template <class Derived>
class Neurons {
    public:
        template<class T>
        void Connect(T &other) {
            auto me = static_cast<Derived*>(this);
            for(auto it = me->begin() ; it != me->end() ; ++it) {
                for(auto it2 = other.begin() ; it2  != other.end() ;++it2) {
                    (*it).out.push_back(&(*it2));
                    (*it2).in.push_back(&(*it));
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

class NeuronLayer : public std::vector<Neuron>, public Neurons<NeuronLayer> {
    //std::vector<Neuron*> neurons;
    public:
    NeuronLayer(int number_of_neurons) {
        while (number_of_neurons-- >0) {
            emplace_back(Neuron{});
        }
    }
};

template<class T>
std::ostream& operator<<(std::ostream& console, Neurons<T> &neuronn) {

    if constexpr (std::is_same_v<T, Neuron>) {
        auto neuron = static_cast<Neuron&>(neuronn);
        for (Neuron *n : neuron.in)
            console << n->ID << "\t>\t" << neuron.ID << "*\n";
        for (Neuron *n : neuron.out)
            console << neuron.ID << "*\t>\t" << n->ID << "\n";
    } else if constexpr (std::is_same_v<T, NeuronLayer>) {
        auto neuron = static_cast<NeuronLayer&>(neuronn);
        for(auto& n: neuron) {
            console<<n;
        }
    }

    return console;
}


int main() {

    Neuron single_neuron_1, single_neuron_2;
    NeuronLayer layer1{1}, layer2{2};
    single_neuron_1.Connect(single_neuron_2);
    single_neuron_1.Connect(layer1);
    layer1.Connect(layer2);
    layer2.Connect(single_neuron_2);
    std::cout << single_neuron_1 << "\n";
    std::cout << single_neuron_2 << "\n";
    std::cout << layer1 << "\n";
    std::cout << layer2 << "\n";
    return 0;
}