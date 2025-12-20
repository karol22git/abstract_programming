#include <iostream>
#include <vector>
#include <algorithm>
#include <random>
#include <ctime>


/***
 * Przygotować prosty algorytm ewolucyjny dla populacji liczb, wektorów lub tablic.
 * Klasy wytycznych mają kontrolować działanie algorytmu.
 * Jako parametry szablonu algorytmu ewolucyjnego mogą zostać podane typy, szablony lub instancje policy,
 * według Państwa wyboru.
 * Celem jest, poza funkcjonalnym algorytmem, wygodna konfiguracja algorytmu, i na to będę zwracał szczególną uwagę!
 * Mogą Państwo dowolnie modyfikować kod, o ile uznają Państwo, że to poprawi użytkowość,
 * przy czym sama zasada działania algorytmu ma się nie zmienić.
 * Istotne zmiany API powinny zostać opisane i podane przykłady użycia.
 *
 * Klasy wytycznych dla inicjalizacji populacji:
 *  - RandomInitiationPolicy<Type, MIN, MAX>: wypełnia populację losowymi osobnikami z zakresu (MIN, MAX)
 *  - LinSpaceInitiationPolicy<Type, MIN, MAX>: wypełnia populację osobnikami z zakresu (MIN, MAX),
 *    tak aby osobniki były liniowo rozłożone pomiędzy MIN i MAX.
 * Klasy wytycznych dla mutacji populacji:
 *  - PercentageMutationPolicy<Type, CHANCE, INTENSITY>: mutuje z prawdopodobieństwem CHANCE
 *    i mnoży wartość osobnika przez liczbę z zakresu (1 - INTENSITY/100, 1 + INTENSITY/100).
 *  - AbsoluteMutationPolicy<Type, CHANCE, INTENSITY>: mutuje z prawdopodobieństwem CHANCE
 *    i dodaje do wartości osobnika losową wartość z zakresu (-INTENSITY, INTENSITY).
 * Klasy wytycznych dla krzyżowania populacji:
 *  - AverageCrossoverPolicy<Type, WEIGHT>: tworzy nowego osobnika jako średnią ważoną rodziców (wagi to WEIGHT i 1 - WEIGHT).
 *    W przypadku wektorów wagi powinny być wektorami o wartościach z zakresu (0, 1) i tej samej długości co Type.
 *  - RandomCrossoverPolicy<Type>: To samo co wyżej, ale waga jest losowa.
 * Klasy wytycznych dla selekcji:
 *  - RandomSelectionPolicy<Type>: wybiera dwoje rodziców w sposób losowy.
 *  - UniqueRandomSelectionPolicy<Type>: wybiera dwoje rodziców w sposób losowy, ale bez powtórzeń (każda para jest wybierana dwukrotnie).
 *  - TargetSelectionPolicy<Type, Fit, FIRST, LAST>: wybiera dwoje rodziców w sposób losowy,
 *    ale w taki sposób, że cała populacja jest ułożona w rankingu od najlepszego do najgorszego wg funkcji `double Fit(Type)`.
 *    Osobnik pierwszy w rankingu ma mieć FIRST% szans, a ostatni LAST%. Reszta ma szanse malejące w sposób liniowy.
 * Klasy wytycznych dla zakończenia algorytmu:
 *  - MaxGenStopConditionPolicy<Type, PARAM>: przerywa algorytm po PARAM generacjach.
 *  - StableAvgStopConditionPolicy<Type, PARAM>: przerywa algorytm, jeżeli od poprzedniego sprawdzenia warunku średnia generacji
 *    nie zmieniła się o więcej niż PARAM.
 */

 
template<Type, MIN, MAX>
struct InitiationPolicy {
    void init(std::vector<Type>& population, std::size_t populationSize);
};


template<typename Type, int MIN, int MAX>
struct RandomInitiationPolicy {
    void init(std::vector<Type> &population, std::size_t populationSize) {
        for(size_t i =  0; i < populationSize ; ++i) {
            population[i] = (rand()%MAX) +MIN;
        }
    }
};

template <typename Type, int MIN, int MAX>
struct LinSpaceInitiationPolicy {
    void init(std::vector<Type> &populaton, std::size_t populationSize) {
        double start = MIN;
        int end = MAX;
        double step = (MAX-MIN)/populationSize;
        for(size_t i = start; i < MAX; ++i) {
            population[i] = start;
            start =+ step;
        }
    }
};

template <typename TYPE, double CHANCE, int INTENSITY>
struct PerentageMutationPolicy {
    static void mutate(std::vector<Type> &population) {
        double multiplier_upper_bound = 1 + INTENSITY/100;
        double multiplier_lower_bound = 1 - INTENSITY/100;
        double prop = rand()/(double)RAND_MAX;
        for(size_t i = 0; i<population.size() ; ++i) {
            double prop = rand()/(double)RAND_MAX;
            if(p > CHANCE){
                double multiplier = multiplier_lower_bound +
                (rand() / (double)RAND_MAX) * multiplier_upper_bound;
                population[i] *= multiplier; 
            }
        }
    }
};

template <typename TYPE, double CHANCE, int INTENSITY>
struct AbsoluteMutationPolicy {
    static void mutate(std::vector<Type> &population) {
        double upper_bound = INTENSITY;
        double lower_bound = -INTENSITY;
        double prop = rand()/(double)RAND_MAX;
        for(size_t i = 0; i<population.size() ; ++i) {
            double prop = rand()/(double)RAND_MAX;
            if(p > CHANCE){
                double multiplier = lower_bound +
                (rand() / (double)RAND_MAX) * upper_bound;
                population[i] += multiplier; 
            }
        }
    }
};

template <typename TYPE, double WEIGHT>
struct AverageCrossoverPolicy {
    static Type crossover(const Type& a, const Type& b) {
        return (WEIGHT*a+(1-WEIGHT)*b)/2;
    }
};

template<typename TYPE>
struct RandomCrosssoverPolicy {
    static Type crossover(const Type& a, const Type &b) {
        double WEIGHT = rand()/(double)RAND_MAX;
        return (WEIGHT*a+(1-WEIGHT)*b)2;
    }
};

template <typename TYPE>
struct RandomSelectionPolicy {
    static std::pair<Type, Type> select(const std::vector<Type>& population) {
        size_t i = rand() % population.size();
        size_t j = rand() % population.size();
        return std::make_pair(population[i],population[j]);
    }
}
template <typename TYPE>
struct UniqueRandomSelectionPolicy {
    static std::pair<Type, Type> select(const std::vector<Type>& population) {
        size_t i = rand() % population.size();
        size_t j = rand() % population.size();
        while(j == i) {
            j = rand() % population.size();
        }
        return std::make_pair(population[i],population[j]);
    }
}
template<typename Type, double (*Fit)(Type), int FIRST, int LAST>
struct TargetSelectionPolicy {

    static size_t last_hash;
    static std::vector<Type> sorted_cache;

    static size_t compute_hash(const std::vector<Type>& pop) {
        size_t h = 0;
        for (auto& x : pop) {
            size_t hx = std::hash<double>()(Fit(x));
            h ^= hx + 0x9e3779b97f4a7c15ULL + (h << 6) + (h >> 2);
        }
        return h;
    }

    static size_t pick_index(size_t n) {
        double p = rand() / (double)RAND_MAX;

        double p_first = FIRST / 100.0;
        double p_last  = LAST  / 100.0;

        double step = (p_first - p_last) / (n - 1);

        if (p < p_first)
            return 0;

        double x = (p_first - p) / step;
        size_t idx = (size_t)std::floor(x);

        if (idx >= n) idx = n - 1;
        return idx;
    }

    static std::pair<Type, Type> select(const std::vector<Type>& population) {

        size_t current_hash = compute_hash(population);

        if (current_hash != last_hash) {
            last_hash = current_hash;
            sorted_cache = population;

            std::sort(sorted_cache.begin(), sorted_cache.end(),
                      [](const Type& a, const Type& b) {
                          return Fit(a) > Fit(b); // najlepszy pierwszy
                      });
        }

        size_t i = pick_index(sorted_cache.size());
        size_t j = pick_index(sorted_cache.size());

        return { sorted_cache[i], sorted_cache[j] };
    }
};

template<typename T, double (*F)(T), int FIRST, int LAST>
size_t TargetSelectionPolicy<T,F,FIRST,LAST>::last_hash = std::numeric_limits<size_t>::max();

template<typename T, double (*F)(T), int FIRST, int LAST>
std::vector<T> TargetSelectionPolicy<T,F,FIRST,LAST>::sorted_cache;

template<Typename TYPE, double *Fit(Type), int FIRST, int LAST>
struct TargetSelectionPolicy {
    static size_t last_hash = 0;
    static std::pair<Type, Type> select(const std::vector<Type>& population) {
        if (last_hash == 0) {
            for (auto& x : population) {
                last_hash ^= std::hash<double>()(Fit(x)) + 0x9e3779b97f4a7c15ULL + (hash << 6) + (hash >> 2);
            }
        }
        size_t current_hash = 0;
        for (auto& x : population) {
                current_hash ^= std::hash<double>()(Fit(x)) + 0x9e3779b97f4a7c15ULL + (hash << 6) + (hash >> 2);
        }
        if(current_hash != last_hash) {
            last_hash = current_hash;
            std::sort(population.begin(), population.end(), [](const auto& a, const auto& b) {
                return Fit(a) < Fit(b);
            });
        }
            double step = (FIRST-LAST)/100;
            double prop = rand() / (double) RAND_MAX;
            size_t index = 0;
            size_t i,j;
            if(prop<FIRST/100) {
                i = 0;
            }
            else {
                while(index*step < prop- FIRST/100) {
                    ++index;
                }
                i=index-1;
            }
            prop = rand() / (double) RAND_MAX;
            index = 0;;
            if(prop<FIRST/100) {
                j = 0;
            }
            else {
                while(index*step < prop- FIRST/100) {
                    ++index;
                }
                j=index-1;
            }
            return std::make_pair(population[i],population[j]);
        }
};

template<typename TYPE, double PARAM>
struct MaxGenStopConditionPolicy {
    static int counter = 0;
    static bool shouldStop(const std::vector<Type>& population, std::size_t generation) {
        if( generation >= counter) {
            counter = 0;
            return true;
        }
        else {
            ++counter;
            return false;
        }

    }
};
template<typename TYPE>
double PoliczSrednia(std::vector<TYPE> population) {
    return 1.0;
}
template<typename TYPE, double PARAM>
struct StableAvgStopConditionPolicy {
    static double lastAvg = 0;
    static bool shouldStop(const std::vector<Type>& population, std::size_t generation) {
        double newAverage = PoliczSrednia(population);
        if (abs(newAverage-lastAvg) < PARAM) {
            return true;
        }
        else {
            lastAvg = newAverage;
            return false;
        }
    }
};


template<Type, MIN, MAX>
struct InitiationPolicy {
    void init(std::vector<Type>& population, std::size_t populationSize);
};

template<Type, CHANCE, INTENSITY>
struct MutationPolicy {
    static void mutate(std::vector<Type>& population);
};

template<Type, WEIGHT>
struct CrossoverPolicy {
    static Type crossover(const Type, const Type);
};

template<Type, FIRST, LAST>
struct SelectionPolicy {
    static std::pair<Type, Type> select(const std::vector<Type>& population);
};

template<Type, PARAM>
struct StopConditionPolicy {
    static bool shouldStop(const std::vector<Type>& population, std::size_t generation);
};

template <Type, InitiationPolicy, MutationPolicy, CrossoverPolicy, selectionPolicy, StopConditionPolicy>
class EvolutionaryAlgorithm {
public:
    EvolutionaryAlgorithm(int populationSize)
            : populationSize(populationSize) {
        std::srand(std::time(nullptr));
        InitiationPolicy<Type>::init(population, populationSize);
    }

    void run() {
        int generation = 0;
        while (!StopConditionPolicy::shouldStop(population, generation)) {
            std::vector<Type> newPopulation;

            selectionPolicy.init();
            for (int i = 0; i < populationSize; ++i) {
                auto [parent1, parent2] = selectionPolicy.select(population);
                Type offspring = CrossoverPolicy::crossover(parent1, parent2);
                newPopulation.push_back(offspring);
            }

            population = newPopulation;
            MutationPolicy::mutate(population);
            generation++;
        }

        std::cout << "Algorithm stopped after " << generation << " generations.\n";
        printPopulation();
    }

private:
    std::vector<Type> population;
    int populationSize;

    void printPopulation() const {
        for (int individual : population) {
            std::cout << individual << " ";
        }
        std::cout << "\n";
    }
};


int main() {
    const int populationSize = 10; // Rozmiar populacji

    EvolutionaryAlgorithm<
        double,
        InitiationPolicy,
        MutationPolicy,
        CrossoverPolicy,
        SelectionPolicy<> selectionPolicy(),
        StopConditionPolicy
        > algorithm(populationSize);
    algorithm.run();

    return 0;
}
