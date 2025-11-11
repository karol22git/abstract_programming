#include <vector>
class Neuron;
class NeuronLayer;

class NeuronLayer : public std::vector<Neuron> {
    //std::vector<Neuron*> neurons;
    public:
    NeuronLayer(int number_of_neurons) {
        while (number_of_neurons-- >0) {
            emplace_back(Neuron{});
        }
    }
};

template <class Derived>
class Neurons {
    public:
    void Connect(NeuronLayer &other) {
        for(Neuron* member: nLayer) {
            for(Neuron* other_member: other) {
                member->out.push_back(other_member);
                other_member->in.push_back(member);
            }
        }
    }
    Neurons(int n): nLayer(n){}
    std::vector<Neurons*> in;
    std::vector<Neurons*> out;
    NeuronLayer nLayer;
};

class Neuron: public Neurons<Neuron> {
    std::vector<Neuron*> in;
    std::vector<Neuron*> out;
    unsigned int ID;
    Neuron(): Neurons(0) {
        static int ID = 1;
        this->ID = ID++;
    }
    Neuron* begin() {return this;}
    Neuron* end() {return this+1;}
};

int main() {
    Neurons<Neuron> layer1(5);
    Neurons<Neuron> layer2(4);
    //layer1.Connect(layer2);
    return 0;
}