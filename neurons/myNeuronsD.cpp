#include <iostream>
#include <vector>
//Derived to layer, derived2 to neuron
class NeuronLayer;
template <class Derived, class Derived2>
class Neurons {
    public:
        void Connect(Derived &other) {
            auto me = static_cast<Derived*>(this);
            //auto me_as_vector = static_cast<std::vector<Derived2>>(*me);
            //auto other_as_vector = static_cast<std::vector<Derived2>>(other);
            for(Derived2& neuron: *me/*me_as_vector*/) {
                for(Derived2& other_neuron: other) {
                    neuron.out.push_back(&other_neuron);
                    other_neuron.in.push_back(&neuron);
                }
            }
        }
        void Connect(Derived2 &other) {
            auto me = static_cast<Derived2*>(this);
            me->out.push_back(&other);
            other.in.push_back(me);
        }
};

class Neuron: public Neurons<NeuronLayer, Neuron> {
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

class NeuronLayer : public std::vector<Neuron>, public Neurons<NeuronLayer, Neuron> {
    //std::vector<Neuron*> neurons;
    public:
    NeuronLayer(int number_of_neurons) {
        while (number_of_neurons-- >0) {
            emplace_back(Neuron{});
        }
    }
};

//template<typename T>
//std::ostream &operator<<(std::ostream &console, /*Neurons<NeuronLayer, Neuron>*/T &neurons) {
//    //???
//    //auto me = static_cast<struct NeuronLayer&>(neurons);
//    //for(auto it = me.begin(); it< me.end(); ++it) {
//    //    for (Neuron *n : (*it).in)
//    //        console << n->ID << "\t>\t" << (*it).ID << "*" << "\n";
//    //    for (Neuron *n : (*it).out)
//    //        console << (*it).ID << "*\t>\t" << n->ID << "\n";
//    //}
//    if constexpr (std::is_same_v<T, Neuron>) {
//        auto me = static_cast<Neuron*>(&neurons);
//        for (Neuron *n : me->in)
//            console << n->ID << "\t>\t" << me->ID << "*" << "\n";
//        for (Neuron *n : me->out)
//            console << me->ID << "*\t>\t" << n->ID << "\n";
//        //return console;
//    } else if constexpr (std::is_same_v<T, NeuronLayer>) {
//        auto me = static_cast<NeuronLayer*>(&neurons);
//    for(Neuron& neuron: *me/*me_as_vector*/) {
//        for(Neuron* m: neuron.out) {
//            console<<neuron.ID<<"* > "<< m->ID<<"\n";
//        }
//        for(Neuron* m: neuron.in) {
//            console<<m->ID<<" > *"<< neuron.ID<<"\n";
//        }
//    }
//    }
//    //auto me = static_cast<NeuronLayer*>(&neurons);
//    //for(Neuron& neuron: *me/*me_as_vector*/) {
//    //    for(Neuron* m: neuron.out) {
//    //        std::cout<<neuron.ID<<"* > "<< m->ID<<std::endl;
//    //    }
//    //    for(Neuron* m: neuron.in) {
//    //        std::cout<<m->ID<<" > *"<< neuron.ID<<std::endl;
//    //    }
//    //}
//    return console;
//}

//template<typename T>
//std::ostream &operator<<(std::ostream &console, T &neurons) {
//    if constexpr (std::is_same_v<T, Neuron>) {
//        // neurons JUŻ JEST Neuron& - nie trzeba rzutować!
//        for (Neuron *n : neurons.in)
//            console << n->ID << "\t>\t" << neurons.ID << "*\n";
//        for (Neuron *n : neurons.out)
//            console << neurons.ID << "*\t>\t" << n->ID << "\n";
//    } else if constexpr (std::is_same_v<T, NeuronLayer>) {
//        // neurons JUŻ JEST NeuronLayer&
//        for(Neuron& neuron : neurons) {
//            for(Neuron* m : neuron.out)
//                console << neuron.ID << "* > " << m->ID << "\n";
//            for(Neuron* m : neuron.in)
//                console << m->ID << " > *" << neuron.ID << "\n";
//        }
//    }
//    return console;
//}

std::ostream& operator<<(std::ostream& console, Neuron& neuron) {
    for (Neuron *n : neuron.in)
        console << n->ID << "\t>\t" << neuron.ID << "*\n";
    for (Neuron *n : neuron.out)
        console << neuron.ID << "*\t>\t" << n->ID << "\n";
    return console;
}

std::ostream& operator<<(std::ostream& console, NeuronLayer& layer) {
    for(Neuron& neuron : layer) {
        console << neuron;  // używa operatora dla Neuron
    }
    return console;
}
int main() {
    //NeuronLayer layer1(5);
    //NeuronLayer layer2(4);
    //auto layer1_cast = static_cast<Neurons<NeuronLayer,Neuron>*>(&layer1);
    //layer1_cast->Connect(layer2);//(layer2_cast);
    ////auto layer2_cast = static_cast<std::vector<Neuron>>(layer2);
    ////for(Neuron n: layer2_cast) {
    ////    for(Neuron* m: n.out) {
    ////        std::cout<<n.ID<<"* > "<< m->ID<<std::endl;
    ////    }
    ////    for(Neuron* m: n.in) {
    ////        std::cout<<m->ID<<" > *"<< n.ID<<std::endl;
    ////    }
    ////}
    //std::cout<<layer2;
    //std::cout<<"sie"<<std::endl;
    Neuron single_neuron_1, single_neuron_2;
    //NeuronLayer layer_1{1}, layer_2{2};
    NeuronLayer layer1{1}, layer2{2};
    single_neuron_1.Connect(single_neuron_2);
    layer1.Connect(layer2);
    std::cout << single_neuron_1 << "\n";
    std::cout << single_neuron_2 << "\n";
    std::cout << layer1 << "\n";
    std::cout << layer2 << "\n";
    return 0;
}