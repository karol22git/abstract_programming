#include <vector>
class Neuron;
class NeuronLayer;

template <class Derived>
class Neurons {
    public:
    void Connect(Derived &other) {
        auto me = static_cast<Derived*>(this);
        auto me_as_vector = static_cast<std::vector<Neuron>>(me);
        auto other_as_vector = static_cast<std::vector<Neuron>>(other);
        //auto he = static_cast<Derived*>(other);
        for(Neuron neuron: me_as_vector) {
            for(Neuron other_neuron: other_as_vector) {
                neuron.out.push_back(other_neuron);
                other_neuron.in.push_back(neuron);
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

int main() {

    NeuronLayer layer1(5);
    NeuronLayer layer2(4);
    auto layer1_cast = static_cast<Neurons<NeuronLayer>*>(&layer1);
    //auto layer2_cast = static_cast<Neurons<NeuronLayer>*>(&layer2);
    //layer1.Connect(layer2);
    layer1_cast->Connect(layer2);//(layer2_cast);
    return 0;
}