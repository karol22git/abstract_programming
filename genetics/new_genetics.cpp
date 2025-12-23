#include <iostream>
#include <vector>
#include <algorithm>
#include <random>
#include <ctime>
#include <array>
#include <type_traits>
#include <cmath>
#include <functional>

template<typename U>
struct is_std_array : std::false_type {};

template<typename U, size_t N>
struct is_std_array<std::array<U, N>> : std::true_type {};


template<typename T>
struct EvolutionTrait {
    static T addition(const T&, const T&);
    static T addition(const T&, double);
    static T multiplication(const T&, double);
};

template<typename T>
struct PrintTrait {
    static void print(const T& x) {
        std::cout << x;
    }
};

template<>
struct PrintTrait<std::vector<double>> {
    static void print(const std::vector<double>& v) {
        std::cout << "[ ";
        for (double x : v)
            std::cout << x << " ";
        std::cout << "]";
    }
};

template<size_t N>
struct PrintTrait<std::array<double, N>> {
    static void print(const std::array<double, N>& a) {
        std::cout << "[ ";
        for (double x : a)
            std::cout << x << " ";
        std::cout << "]";
    }
};

template<typename Type, int MIN, int MAX>
struct InitiationPolicy {
    void init(std::vector<Type>& population, std::size_t populationSize);
};

template<typename Type,double CHANCE, int INTENSITY>
struct MutationPolicy {
    static void mutate(std::vector<Type>& population);
};

template<typename Type, double WEIGHT>
struct CrossoverPolicy {
    static Type crossover(const Type, const Type);
};
template<typename Type, int FIRST, int LAST>
struct SelectionPolicy {
    static std::pair<Type, Type> select(const std::vector<Type>& population);
};

template<typename Type, int PARAM>
struct StopConditionPolicy {
    static bool shouldStop(const std::vector<Type>& population, std::size_t generation);
};
 
template<>
struct EvolutionTrait<double> {
    static double random(int MIN, int MAX) {
        return MIN+(rand()/(double)RAND_MAX)*(MAX-MIN);
    }

    static double constant(double v) {
        return v;
    }

    static double addition(double individual, double b) {
        return individual+b;
    }

    static double multiplication(double individual, double b) {
        return individual*b;
    }
};

template<>
struct EvolutionTrait<std::vector<double>> {

    static std::vector<double> random(int MIN, int MAX, size_t size) {
        std::vector<double> result(size);
        for (double &x: result)
            x = EvolutionTrait<double>::random(MIN, MAX);
        return result;
    }

    static std::vector<double> constant(double v, size_t size) {
        return std::vector<double>(size, v);
    }

    static std::vector<double> addition(const std::vector<double>& a, const std::vector<double>& b)
    {
        std::vector<double> result(a.size());
        for (size_t i = 0; i < a.size(); ++i)
            result[i] = a[i]+b[i];
        return result;
    }

    static std::vector<double> addition(const std::vector<double>& individual, double b) {
        std::vector<double> result(individual.size());
        for (size_t i = 0; i < individual.size(); ++i)
            result[i] = individual[i]+b;
        return result;
    }

    static std::vector<double> multiplication(const std::vector<double>& individual, double b) {
        std::vector<double> result(individual.size());
        for (size_t i = 0; i < individual.size(); ++i)
            result[i] = individual[i]*b;
        return result;
    }
};

template<size_t N>
struct EvolutionTrait<std::array<double, N>> {

    static std::array<double, N> random(int MIN, int MAX) {
        std::array<double, N> a;
        for (double &x : a)
            x = EvolutionTrait<double>::random(MIN, MAX);
        return a;
    }

    static std::array<double, N> constant(double v) {
        std::array<double, N> a;
        a.fill(v);
        return a;
    }

    static std::array<double, N> addition(const std::array<double, N>& a, const std::array<double, N>& b) {
        std::array<double, N> result;
        for (size_t i = 0; i < N; ++i)
            result[i] = a[i]+b[i];
        return result;
    }

    static std::array<double, N> addition(const std::array<double, N>& individual, double b) {
        std::array<double, N> result;
        for (size_t i = 0; i < N; ++i)
            result[i] = individual[i]+b;
        return result;
    }

    static std::array<double, N> multiplication(const std::array<double, N>& individual, double b)
    {
        std::array<double, N> result;
        for (size_t i = 0; i < N; ++i)
            result[i] = individual[i]*b;
        return result;
    }
};



template<typename T, int MIN, int MAX>
struct RandomInitiationPolicy {
    static void init(std::vector<T>& population, size_t size) {
        for (auto& x: population)
            x = EvolutionTrait<T>::random(MIN, MAX);
    }
};


template<typename T, int MIN, int MAX>
struct LinSpaceInitiationPolicy {
    static void init(std::vector<T>& population, size_t size) {
        double step = (MAX - MIN) / double(size - 1);
        double value = MIN;

        for (size_t i = 0; i < size; ++i) {
            if constexpr (std::is_same_v<T, double>) {
                population[i] = value;
            }
            else if constexpr (std::is_same_v<T, std::vector<double>>) {
                population[i] = EvolutionTrait<T>::constant(value, population[i].size());
            }
            else if constexpr (is_std_array<T>::value) {
                population[i] = EvolutionTrait<T>::constant(value);
            }
            value += step;
        }
    }
};


template<typename T, double CHANCE, int INTENSITY>
struct PercentageMutationPolicy {
    static void mutate(std::vector<T>& population) {
        double multiplier_upper_bound = 1.0 + INTENSITY/100.0;
        double multiplier_lower_bound = 1.0 - INTENSITY/100.0;
        for (auto& x: population) {
            double prob = (rand()/(double)RAND_MAX);
            if (prob < CHANCE) {
                double multiplier = multiplier_lower_bound +
                (rand() / (double)RAND_MAX) * multiplier_upper_bound;
                x = EvolutionTrait<T>::multiplication(x,multiplier);
            }
        }
    }
};

template<typename T, double CHANCE, int INTENSITY>
struct AbsoluteMutationPolicy {
    static void mutate(std::vector<T>& population) {
        double upper_bound = INTENSITY;
        double lower_bound = -INTENSITY;
        for (T& x: population) {
            double p = rand() / (double)RAND_MAX;
            if (p < CHANCE) {
                double multiplier = lower_bound +(rand() / (double)RAND_MAX) * upper_bound;
                x = EvolutionTrait<T>::addition(x, multiplier);
            }
        }
    }
};
template<typename T>
double individualAverage(const T& x) {
    double sum = 0.0;
    for (double v: x) sum += v;
    return sum/x.size();
}

template<>
double individualAverage(const double& x) {
    return x;
}
template<typename T>
double populationAverage(const std::vector<T>& population) {
    double sum = 0.0;
    for (const T& x: population)
        sum += individualAverage(x);
    return sum / population.size();
}

template<typename T, int MAXGEN>
struct MaxGenStopConditionPolicy {
    static bool shouldStop(const std::vector<T>&, size_t gen) {
        return gen >= MAXGEN;
    }
};

template<typename T, double PARAM>
struct StableAvgStopConditionPolicy {
    static double lastAvg;
    static bool firstCall;
    static bool shouldStop(const std::vector<T>& population, size_t) {
        double avg = populationAverage(population);
        if (firstCall) {
            firstCall = false;
            lastAvg = avg;
            return false;
        }
        if (std::abs(avg - lastAvg) < PARAM) {
            return true;
        }
        lastAvg = avg;
        return false;
    }
};

template<typename T, double PARAM>
double StableAvgStopConditionPolicy<T, PARAM>::lastAvg = 0.0;

template<typename T, double PARAM>
bool StableAvgStopConditionPolicy<T, PARAM>::firstCall = true;

template<typename T>
struct RandomSelectionPolicy {
    static std::pair<T,T> select(const std::vector<T>& population) {
        size_t i = rand() % population.size();
        size_t j = rand() % population.size();
        return std::make_pair(population[i],population[j]);
    }
};

template<typename T>
struct UniqueRandomSelectionPolicy {
    static std::pair<T, T> select(const std::vector<T>& population) {
        size_t i = rand() % population.size();
        size_t j = rand() % population.size();

        while (j == i) {
            j = rand() % population.size();
        }
        return std::make_pair(population[i],population[j]);
    }
};
//nie moje
template<typename T, double (*Fit)(const T&), int FIRST, int LAST>
struct TargetSelectionPolicy {

    static size_t lastHash;
    static std::vector<T> sortedCache;

    static size_t computeHash(const std::vector<T>& pop) {
        size_t h = 0;
        for (const T& x: pop) {
            size_t fx = std::hash<double>()(Fit(x));
            h ^= fx + 0x9e3779b97f4a7c15ULL + (h << 6) + (h >> 2);
        }
        return h;
    }

    static size_t pickIndex2(size_t n) {
    double p = rand() / (double)RAND_MAX;

    double p_first = FIRST / 100.0;
    double p_last  = LAST  / 100.0;

    double step = (p_first - p_last) / (n - 1);

    double cumulative = 0.0;

    for (size_t i = 0; i < n; i++) {
        double prob_i = p_first - i * step;
        cumulative += prob_i;

        if (p < cumulative)
            return i;
    }

    return n - 1;
}

static size_t pickIndex(size_t n) {
    if (n == 0) return 0;
    if (n == 1) return 0;

    double p = rand() / (double)RAND_MAX;
    double p_first = FIRST / 100.0;
    double p_last = LAST / 100.0;
    double step = (p_first - p_last) / (n - 1);
    double total_sum = n * (p_first + p_last) / 2.0;
    double cumulative = 0.0;
    
    for (size_t i = 0; i < n; i++) {
        double prob_i = (p_first - i * step) / total_sum;
        cumulative += prob_i;
        
        if (p < cumulative || (i == n-1 && p <= 1.0)) {
            return i;
        }
    }
    
    return n - 1;
}

    static std::pair<T, T> select(const std::vector<T>& population) {

        size_t h = computeHash(population);

        if (h != lastHash) {
            lastHash = h;
            sortedCache = population;

            std::sort(sortedCache.begin(), sortedCache.end(),
                      [](const T& a, const T& b) {
                          return Fit(a) > Fit(b);
                      });
        }

        size_t i = pickIndex(sortedCache.size());
        size_t j = pickIndex(sortedCache.size());

         return std::make_pair(population[i],population[j]);
    }
};

template<typename T, double (*Fit)(const T&), int FIRST, int LAST>
size_t TargetSelectionPolicy<T, Fit, FIRST, LAST>::lastHash = 0;

template<typename T, double (*Fit)(const T&), int FIRST, int LAST>
std::vector<T> TargetSelectionPolicy<T, Fit, FIRST, LAST>::sortedCache;


template<typename T>
struct RandomCrossoverPolicy {
    static T crossover(const T& a, const T& b) {
        double weight = rand() / (double)RAND_MAX;
        T mom = EvolutionTrait<T>::multiplication(a, weight);
        T dad = EvolutionTrait<T>::multiplication(b, 1.0 - weight);
        return EvolutionTrait<T>::addition(mom, dad);
    }
};

template<typename T, double WEIGHT>
struct AverageCrossoverPolicy {
    static T crossover(const T& a, const T& b) {
        return EvolutionTrait<T>::addition(
            EvolutionTrait<T>::multiplication(a, WEIGHT),
            EvolutionTrait<T>::multiplication(b, 1.0 - WEIGHT)
        );
    }
};


template <typename Type, class InitiationPolicy, class MutationPolicy, class CrossoverPolicy, class  SelectionPolicy, class StopConditionPolicy>
class EvolutionaryAlgorithm {
public:
    EvolutionaryAlgorithm(int populationSize)
            : populationSize(populationSize) {
        std::srand(std::time(nullptr));
        population.resize(populationSize);
        InitiationPolicy::init(population, populationSize);
    }

    void run() {
        int generation = 0;
        while (!StopConditionPolicy::shouldStop(population, generation)) {
            std::vector<Type> newPopulation;

            for (int i = 0; i < populationSize; ++i) {
                auto [parent1, parent2] = SelectionPolicy::select(population);
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
    for (const auto& individual : population) {
        PrintTrait<Type>::print(individual);
        std::cout << " ";
    }
    std::cout << "\n";
}

};

void testTraits() {
    std::cout << "=== TRAITS TEST ===\n";

    double a = 2.0;
    double b = 3.0;

    std::cout << "double add: " << EvolutionTrait<double>::addition(a,b) << "\n";
    std::cout << "double mul: " << EvolutionTrait<double>::multiplication(a,2.0) << "\n";

    std::vector<double> v1 = {1,2,3};
    std::vector<double> v2 = {4,5,6};

    auto vadd = EvolutionTrait<std::vector<double>>::addition(v1,v2);
    auto vmul = EvolutionTrait<std::vector<double>>::multiplication(v1,2.0);

    PrintTrait<std::vector<double>>::print(vadd); std::cout << "\n";
    PrintTrait<std::vector<double>>::print(vmul); std::cout << "\n";

    std::array<double,3> a1 = {1,2,3};
    std::array<double,3> a2 = {4,5,6};

    auto aadd = EvolutionTrait<std::array<double,3>>::addition(a1,a2);
    auto amul = EvolutionTrait<std::array<double,3>>::multiplication(a1,2.0);

    PrintTrait<std::array<double,3>>::print(aadd); std::cout << "\n";
    PrintTrait<std::array<double,3>>::print(amul); std::cout << "\n";

    std::cout << "=== END TRAITS TEST ===\n\n";
}
double FitDouble(const double& x) {
    return -std::abs(x); // im bliżej 0, tym lepiej
}

void testTargetSelection() {
    std::cout << "\n=== TARGET SELECTION POLICY TEST ===\n";

    using T = double;

    // Tworzymy populację z wartościami o różnej "jakości"
    std::vector<T> population = { 5.0, 3.0, 1.0, -2.0, -0.5, 4.0 };

    // Polityka selekcji: FIRST=60%, LAST=5%
    using Sel = TargetSelectionPolicy<T, FitDouble, 60, 5>;

    // 1. Pierwsze wywołanie — powinno posortować populację
    auto [p1, p2] = Sel::select(population);

    std::cout << "First selection: " << p1 << " , " << p2 << "\n";

    // 2. Sprawdzamy czy cache działa — drugie wywołanie NIE powinno sortować ponownie
    size_t oldHash = Sel::lastHash;
    auto [p3, p4] = Sel::select(population);

    std::cout << "Second selection: " << p3 << " , " << p4 << "\n";

    if (Sel::lastHash == oldHash)
        std::cout << "Cache OK (no re-sort)\n";
    else
        std::cout << "Cache ERROR (re-sort happened)\n";

    // 3. Sprawdzamy czy pickIndex wybiera sensowne indeksy
    std::cout << "Testing pickIndex distribution...\n";
    std::vector<int> counts(population.size(), 0);

    for (int i = 0; i < 5000; i++) {
        size_t idx = Sel::pickIndex(population.size());
        counts[idx]++;
    }

    std::cout << "Index distribution:\n";
    for (size_t i = 0; i < counts.size(); i++) {
        std::cout << "  idx " << i << ": " << counts[i] << "\n";
    }

    std::cout << "=== END TARGET SELECTION POLICY TEST ===\n\n";
}

int main() {
    using T = double;
    EvolutionaryAlgorithm<
        T,
        RandomInitiationPolicy<T, -5, 5>,
        PercentageMutationPolicy<T, 0.5, 10>,
        RandomCrossoverPolicy<T>,
        RandomSelectionPolicy<T>,
        MaxGenStopConditionPolicy<T, 10>
    > alg(5);
    alg.run();
return 0;
}