#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <numeric>
#include <cmath>
#include <cstdlib>
#include <ctime>
#include <chrono>
#include <execution>
#include <random>
#include <string>
#include <filesystem>
#include <iomanip>

template <typename T>
std::vector<size_t> argsort(const std::vector<T>& array) {
    // Create a vector of indices from 0 to N-1
    std::vector<size_t> indices(array.size());
    std::iota(indices.begin(), indices.end(), 0); // Fills indices with 0, 1, 2, ...

    // Sort the indices vector based on the values in the original 'array'
    std::sort(indices.begin(), indices.end(),
        [&array](size_t i, size_t j) {
            return array[i] < array[j];
        });

    return indices;
}

struct Vec3 {
    double x, y, z;
    
    Vec3() : x(0), y(0), z(0) {}
    Vec3(double x, double y, double z) : x(x), y(y), z(z) {}

    double min() const {
        return std::min(x, std::min(y, z));
    }
    
    double max() const {
        return std::max(x, std::max(y, z));
    }

    bool operator<(Vec3 other) const {
        return (x < other.x) && (y < other.y) && (z < other.z);
    }
       
    Vec3 operator+(Vec3 other) const {
        return Vec3(x + other.x, y + other.y, z + other.z);
    }

    Vec3 operator-(Vec3 other) const {
        return Vec3(x - other.x, y - other.y, z - other.z);
    }
};

struct Box3D {
    Vec3 min;
    Vec3 max;

    Box3D() : min(0, 0, 0), max(0, 0, 0) {}
    
    Box3D(Vec3 minCoord, Vec3 maxCoord) 
        : min(minCoord), max(maxCoord) {}
    
    Vec3 getSize() const {
        return max - min;
    }
    
    double getVolume() const {
        Vec3 size = getSize();
        return size.x * size.y * size.z;
    }

    bool intersects(Box3D other) const {
        return (min.x < other.max.x && other.min.x < max.x) &&
               (min.y < other.max.y && other.min.y < max.y) &&
               (min.z < other.max.z && other.min.z < max.z);
    }
    
    bool inscribed(Box3D other) const {
        return (min.x >= other.min.x && max.x <= other.max.x) &&
               (min.y >= other.min.y && max.y <= other.max.y) &&
               (min.z >= other.min.z && max.z <= other.max.z);
    }

    bool fits(Box3D other) const {
        return getSize() < (other.getSize());
    }
       
    bool operator<(Box3D other) const {
        return min < other.min && max < other.max;
    }

    bool operator==(Box3D other) const {
        return (min.x == other.min.x) && (min.y == other.min.y) && (min.z == other.min.z) &&
               (max.x == other.max.x) && (max.y == other.max.y) && (max.z == other.max.z);
    }
};

struct Bin {
    Vec3 dimensions = Vec3();
    double empty_volume = 0;
    double max_ems_volume = 0;
    std::vector<Box3D> EMSs = std::vector<Box3D>();
    std::vector<Box3D> boxes = std::vector<Box3D>();

    Bin(Vec3 dim) : dimensions(dim) {
        boxes.reserve(15);
        EMSs.reserve(30);
        EMSs.emplace_back(Vec3(), dimensions);
        empty_volume = dimensions.x * dimensions.y * dimensions.z;
        max_ems_volume = empty_volume;
    }

    void update(Box3D boxToPlace, Box3D selectedEMS, double minVol, double minDim) {
        auto box = Box3D(selectedEMS.min, selectedEMS.min + boxToPlace.max);
        boxes.push_back(box);
        empty_volume -= box.getVolume();
        auto EMSsCopy = EMSs;
        for (auto ems : EMSsCopy) {
            if (box.intersects(ems)) {
                auto it = std::find(EMSs.begin(), EMSs.end(), ems);
                if (it != EMSs.end()) {
                    EMSs.erase(it);
                }

                Box3D newEMSs[] = {
                    Box3D(Vec3(box.max.x, ems.min.y, ems.min.z), Vec3(ems.max.x, ems.max.y, ems.max.z)),
                    Box3D(Vec3(ems.min.x, box.max.y, ems.min.z), Vec3(ems.max.x, ems.max.y, ems.max.z)),
                    Box3D(Vec3(ems.min.x, ems.min.y, box.max.z), Vec3(ems.max.x, ems.max.y, ems.max.z)),
                };

                for (auto new_ems : newEMSs) {
                    bool isValid = true;

                    // Do not add new EMS having smaller dimensions of the smallest dimension of remaining boxes
                    if (new_ems.getSize().min() < minDim) {
                        isValid = false;
                    }

                    // Do not add new EMS smaller than the volume of remaining boxes
                    else if (new_ems.getVolume() < minVol) {
                        isValid = false;
                    }

                    // Eliminate new EMSs which are totally inscribed by other EMSs
                    else {
                        for (auto otherEms : EMSs) {
                            if (new_ems.inscribed(otherEms)) {
                                isValid = false;
                            }
                        }
                    }

                    if (isValid) {
                        EMSs.push_back(new_ems);
                    }
                }
            }
        }
        max_ems_volume = 0;
        for (auto ems : EMSs) {
            if (ems.getVolume() > max_ems_volume) {
                max_ems_volume = ems.getVolume();
            }
        }
    }

    double load() {
        double total_volume = dimensions.x * dimensions.y * dimensions.z;
        return (total_volume - empty_volume) / total_volume;
    }
};

struct PlacementProcedure {
    Vec3 bin_dimensions;
    std::vector<Bin> bins;
    const std::vector<Box3D> *boxes;
    std::vector<size_t> bps;
    std::vector<double> vbo;
    size_t infeasible_bins = 200;
    bool infeasible = false;

    PlacementProcedure(const std::vector<Box3D> *boxes, Vec3 bin_dimensions, size_t bin_amount, const std::vector<double> *solution)
        : boxes(boxes), bin_dimensions(bin_dimensions) {
        bins = std::vector<Bin>();
        bins.reserve(bin_amount);
        bins.emplace_back(bin_dimensions);

        bps = std::vector<size_t>(solution->size()/2);
        std::iota(bps.begin(), bps.end(), 0);
        std::sort(bps.begin(), bps.end(),
            [solution](size_t i, size_t j) {
                return (*solution)[i] < (*solution)[j];
            });

        vbo = std::vector<double>(solution->begin() + boxes->size(), solution->end());
    }

    void place(bool can_rotate) {
        for (size_t i = 0; i < bps.size(); i++) {
            auto box = (*boxes)[bps[i]];
            bool bin_found = false;
            size_t selected_bin_index = 0;
            Box3D selected_ems;
            for (int k = 0; k < bins.size(); k++) {
                if (bins[k].max_ems_volume < box.getVolume()) {
                    continue;
                }
                auto ems = DFTRC_2(box, bins[k], can_rotate);
                if (ems.getSize().x != 0) {
                    bin_found = true;
                    selected_bin_index = k;
                    selected_ems = ems;
                    break;
                }
            }
            if (!bin_found) {
                bins.emplace_back(bin_dimensions);
                selected_bin_index = bins.size() - 1;
                if (bins.size() > infeasible_bins) {
                    infeasible = true;
                    return;
                }
                selected_ems = *(bins[selected_bin_index].EMSs.begin());
            }

            auto box_orientation = box.getSize();
            if (can_rotate) {
                box_orientation = selectBoxOrientation(vbo[i], box, selected_ems);
            }
            auto mins = eliminationRule(i+1);
            bins[selected_bin_index].update(Box3D(Vec3(0,0,0), box_orientation), selected_ems, mins.first, mins.second);
        }
    }

    Box3D DFTRC_2(Box3D box, const Bin &bin, bool can_rotate) {
        double max_dist = -1;
        Box3D selected_ems;
        auto dim = box.getSize();
        Vec3 orientations[] = {
            Vec3(dim.x, dim.y, dim.z),
            Vec3(dim.x, dim.z, dim.y),
            Vec3(dim.y, dim.x, dim.z),
            Vec3(dim.y, dim.z, dim.x),
            Vec3(dim.z, dim.x, dim.y),
            Vec3(dim.z, dim.y, dim.x),
        };

        if (!can_rotate) {
            for (auto ems : bin.EMSs) {
                if (Box3D(Vec3(), dim).fits(ems)) {
                    return ems;
                }
            }
            return selected_ems;
        }

        for (auto ems : bin.EMSs) {
            if (!can_rotate) {
                if (Box3D(Vec3(), dim).fits(ems)) {
                    selected_ems = ems;
                    break;
                }
            }
            for (int i = 0; i < 6; i++) {
                if (Box3D(Vec3(0,0,0), orientations[i]).fits(ems)) {
                    auto distance = std::pow(bin.dimensions.x - ems.min.x - orientations[i].x, 2) +
                                    std::pow(bin.dimensions.y - ems.min.y - orientations[i].y, 2) +
                                    std::pow(bin.dimensions.z - ems.min.z - orientations[i].z, 2);
                    if (distance > max_dist) {
                        max_dist = distance;
                        selected_ems = ems;
                    }
                }
            }
        }
        return selected_ems;
    }

    Vec3 selectBoxOrientation(double bo, Box3D box, Box3D ems) {
        auto dim = box.getSize();
        Vec3 orientations[] = {
            Vec3(dim.x, dim.y, dim.z),
            Vec3(dim.x, dim.z, dim.y),
            Vec3(dim.y, dim.x, dim.z),
            Vec3(dim.y, dim.z, dim.x),
            Vec3(dim.z, dim.x, dim.y),
            Vec3(dim.z, dim.y, dim.x),
        };
        std::vector<size_t> indexes;
        indexes.reserve(6);
        for (size_t i = 0; i < 6; i++) {
            if (Box3D(Vec3(0,0,0), orientations[i]).fits(ems)) {
                indexes.push_back(i);
            }
        }
        int c = (int) (std::ceil(bo*(indexes.size()))-1) % indexes.size();
        Vec3 selected_bo = orientations[indexes[c]];
        return selected_bo;
    }

    std::pair<double, double> eliminationRule(size_t remaining_start) {
        if (remaining_start == bps.size()) {
            return std::pair<double, double>(0, 0);
        }
        double minVol = 1e10;
        double minDim = 1e10;
        for (size_t i = remaining_start; i < bps.size(); i++) {
            auto box = (*boxes)[bps[i]];
            auto dim = box.max.min();
            if (dim < minDim) {
                minDim = dim;
            }
            auto vol = box.getVolume();
            if (vol < minVol) {
                minVol = vol;
            }
        }
        return std::pair<double, double>(minVol, minDim);
    }

    double evaluate() {
        if (infeasible) {
            return 1e10;
        }
        double leastLoad = 1;
        for (auto &bin : bins) {
            auto load = bin.load();
            if (load < leastLoad) {
                leastLoad = load;
            }
        }
        return bins.size() + leastLoad;
    }
};

struct Result {
    double best_fitness;
    int best_fitness_gen;
};

struct BRKGA {
    Vec3 bin_dimensions;
    const std::vector<Box3D> *boxes;
    size_t num_boxes;

    size_t num_generations;
    size_t num_individuals;
    size_t num_gene;

    size_t num_elites;
    size_t num_mutants;
    double eliteCProb;

    bool can_rotate;

    size_t used_bins;
    //std::vector<double> solution;
    double best_fitness;
    int best_iter;

    std::vector<double> mean;
    std::vector<double> min;

    BRKGA(Vec3 bin_dimensions, const std::vector<Box3D> *boxes)
        : bin_dimensions(bin_dimensions), boxes(boxes) {
            num_boxes = boxes->size();
            num_generations = 200;
            num_individuals = 30*num_boxes;
            num_gene = 2*num_boxes;
            num_elites = static_cast<size_t>(std::floor(0.1*num_individuals));
            num_mutants = static_cast<size_t>(std::floor(0.15*num_individuals));
            eliteCProb = 0.7;
            can_rotate = true;
            mean = std::vector<double>();
            mean.reserve(num_generations);
            min = std::vector<double>();
            min.reserve(num_generations);
    }

/*     std::vector<double> calculate_fitness(const std::vector<std::vector<double>> *population) {
        std::vector<double> fitness_list;
        fitness_list.reserve(population->size());
        for (auto solution : *population) {
            auto decoder = PlacementProcedure(boxes, bin_dimensions, 200, &solution);
            decoder.place();
            fitness_list.push_back(decoder.evaluate());
        }
        return fitness_list;
    } */

    std::vector<double> calculate_fitness(const std::vector<std::vector<double>> *population) {
        std::vector<double> fitness_list(population->size());

        std::transform(
            std::execution::par,
            population->begin(),
            population->end(),
            fitness_list.begin(),
            [&](const std::vector<double>& solution) {
                auto decoder = PlacementProcedure(boxes, bin_dimensions, 200, &solution);
                decoder.place(can_rotate);
                return decoder.evaluate();
            }
        );

        return fitness_list;
    }

    Result fit(bool verbose) {
        std::srand(static_cast<unsigned int>(std::time({})));
        std::vector<std::vector<double>> population;
        population.reserve(num_individuals);
        
        for (int i = 0; i < num_individuals; i++) {
            std::vector<double> gene;
            gene.reserve(num_gene);
            for (int k = 0; k < num_gene; k++) {
                gene.push_back((double)rand()/RAND_MAX);
            }
            population.push_back(gene);
        }

        std::vector<double> fitness_list = calculate_fitness(&population);

        if (verbose) {
            std::cout << "\nInitial Population:\n";
            std::cout << "   ->   population size: " << num_individuals << std::endl;
            std::cout << "   ->   Best Fitness: " <<  *(std::max_element(fitness_list.begin(), fitness_list.end())) << std::endl;
        }

        auto best_fitness_it = std::min_element(fitness_list.begin(), fitness_list.end());
        best_fitness = *(best_fitness_it);
        auto best_solution = population[std::distance(fitness_list.begin(), best_fitness_it)];
        min.push_back(best_fitness);
        mean.push_back(std::accumulate(fitness_list.begin(), fitness_list.end(), 0.0) / fitness_list.size());

        best_iter = 0;
        for (int g = 0; g < num_generations; g++) {
            std::vector<size_t> sorted_indexes = argsort(fitness_list);
            std::vector<std::vector<double>> elites;
            std::vector<double> elite_fitness_list;
            elites.reserve(num_elites);
            elite_fitness_list.reserve(num_elites);
            for (int i = 0; i < num_elites; i++) {
                elites.push_back(population[sorted_indexes[i]]);
                elite_fitness_list.push_back(fitness_list[sorted_indexes[i]]);
            }
            std::vector<std::vector<double>> non_elites;
            non_elites.reserve(num_individuals-num_elites);
            for (size_t i = num_elites; i < num_individuals; i++) {
                non_elites.push_back(population[sorted_indexes[i]]);
            }
            size_t num_offspring = num_individuals - num_elites - num_mutants;
            std::vector<std::vector<double>> offsprings;
            offsprings.reserve(num_offspring);
            for (int i = 0; i < num_offspring; i++) {
                int elite_index = rand() % num_elites;
                int non_elite_index = rand() % (num_individuals - num_elites);
                std::vector<double> offspring;
                offspring.reserve(num_gene);
                for (int gene = 0; gene < num_gene; gene++) {
                    if (((double)rand() / (double)RAND_MAX) < eliteCProb) {
                        offspring.push_back(elites[elite_index][gene]);
                    }
                    else {
                        offspring.push_back(non_elites[non_elite_index][gene]);
                    }
                }
                offsprings.push_back(offspring);
            }
            std::vector<std::vector<double>> mutants;
            mutants.reserve(num_mutants);
            for (int i = 0; i < num_mutants; i++) {
                std::vector<double> mutant;
                mutant.reserve(num_gene);
                for (int k = 0; k < num_gene; k++) {
                    mutant.push_back((double)((double)rand() / (double)RAND_MAX));
                }
                mutants.push_back(mutant);
            }
            std::vector<std::vector<double>> offspring;
            offspring.reserve(mutants.size() + offsprings.size());
            offspring.insert(offspring.end(), mutants.begin(), mutants.end());
            offspring.insert(offspring.end(), offsprings.begin(), offsprings.end());
            auto offspring_fitness_list = calculate_fitness(&offspring);

            std::vector<std::vector<double>> new_population;
            new_population.reserve(num_elites + offsprings.size());
            new_population.insert(new_population.end(), elites.begin(), elites.end());
            new_population.insert(new_population.end(), offsprings.begin(), offsprings.end());
            
            std::vector<double> new_fitness_list;
            new_fitness_list.reserve(elite_fitness_list.size() + offspring_fitness_list.size());
            new_fitness_list.insert(new_fitness_list.end(), elite_fitness_list.begin(), elite_fitness_list.end());
            new_fitness_list.insert(new_fitness_list.end(), offspring_fitness_list.begin(), offspring_fitness_list.end());

            double min_fitness = 999999;
            double fitness_sum = 0;
            for (double fitness : new_fitness_list) {
                if (fitness < best_fitness) {
                    best_iter = g;
                    best_fitness = fitness;
                    best_solution = population[sorted_indexes[0]];
                }
                if (fitness < min_fitness) {
                    min_fitness = fitness;
                }
                fitness_sum += fitness;
            }
            min.push_back(min_fitness);
            mean.push_back(fitness_sum / new_fitness_list.size());

            if (verbose) {
                std::cout << "Generation: " << g << " \t(Best Fitness: " << best_fitness << ")\n";
            }
        }
        used_bins = static_cast<size_t>(std::floor(best_fitness));
        std::cout << "Best Fitness: " << best_fitness << std::endl;
        std::cout << "Generation for Best Fitness: " << best_iter << std::endl;
        return Result{best_fitness, best_iter};
    }
};

double rand_uniform(double a, double b, std::mt19937& gen) {
    std::uniform_real_distribution<double> dist(a, b);
    return dist(gen);
}

struct Instance {
    Vec3 dimensions;
    std::vector<Box3D> boxes;
};

Instance generateInstance(size_t n, int instanceClass) {
    std::mt19937 gen(std::random_device{}());

    double W, H, D;
    W = H = D = 100;
    Instance instance;
    instance.boxes = std::vector<Box3D>();
    instance.boxes.reserve(n);

    if (instanceClass < 6) {
        W = H = D = 100.0f;

        // probabilities
        std::vector<double> probs{0.1, 0.1, 0.1, 0.1, 0.1};
        probs[instanceClass - 1] = 0.6;

        std::discrete_distribution<int> type_dist(
            probs.begin(), probs.end()
        );

        // low / high bounds
        std::vector<std::vector<double>> low = {{
            {1.0f, 2.0f/3.0f * H, 2.0f/3.0f * D},
            {2.0f/3.0f * W, 1.0f, 2.0f/3.0f * D},
            {2.0f/3.0f * W, 2.0f/3.0f * H, 1.0f},
            {0.5f * W, 0.5f * H, 0.5f * D},
            {1.0f, 1.0f, 1.0f}
        }};

        std::vector<std::vector<double>> high = {{
            {0.5f * W, H, D},
            {W, 0.5f * H, D},
            {W, H, 0.5f * D},
            {W, H, D},
            {0.5f * W, 0.5f * H, 0.5f * D}
        }};

        for (int i = 0; i < n; ++i) {
            int type = type_dist(gen);
            instance.boxes.push_back(Box3D(Vec3(0,0,0), Vec3(
                rand_uniform(low[type][0], high[type][0], gen),
                rand_uniform(low[type][1], high[type][1], gen),
                rand_uniform(low[type][2], high[type][2], gen)
            )));
        }
    }
    else if (instanceClass == 6) {
        W = H = D = 10.0f;
        for (int i = 0; i < n; ++i)
            instance.boxes.push_back(Box3D(Vec3(0,0,0), Vec3(
                rand_uniform(1.0f, 10.0f, gen),
                rand_uniform(1.0f, 10.0f, gen),
                rand_uniform(1.0f, 10.0f, gen)
            )));
    }
    else if (instanceClass == 7) {
        W = H = D = 40.0f;
        for (int i = 0; i < n; ++i)
            instance.boxes.push_back(Box3D(Vec3(0,0,0), Vec3(
                rand_uniform(1.0f, 35.0f, gen),
                rand_uniform(1.0f, 35.0f, gen),
                rand_uniform(1.0f, 35.0f, gen)
            )));
    }
    else if (instanceClass == 8) {
        W = H = D = 100.0f;
        for (int i = 0; i < n; ++i)
            instance.boxes.push_back(Box3D(Vec3(0,0,0), Vec3(
                rand_uniform(1.0f, 100.0f, gen),
                rand_uniform(1.0f, 100.0f, gen),
                rand_uniform(1.0f, 100.0f, gen)
            )));
    }

    instance.dimensions = Vec3(W, H, D);
    return instance;
}

int main() {
    size_t n = 50;

    Instance instances[8][10];

    for (int i = 0; i < 8; i++) {
        std::string filename = "data/class-" + std::to_string(i+1) + "-" + std::to_string(n) + ".txt";

        if (std::filesystem::exists(filename)) {
            std::ifstream file(filename);
            double W, H, D;
            file >> W >> H >> D;

            for (int k = 0; k < 10; k++) {
                instances[i][k].dimensions = Vec3(W, H, D);
                instances[i][k].boxes.reserve(n);
                for (int m = 0; m < n; m++) {
                    double x, y, z;
                    file >> x >> y >> z;
                    instances[i][k].boxes.emplace_back(Vec3(), Vec3(x, y, z));
                }
            }
            file.close();
        }
        else {
            std::ofstream file(filename);
            file << std::fixed << std::setprecision(5);
            for (int k = 0; k < 10; k++) {
                instances[i][k] = generateInstance(n, i+1);
            }
            file << instances[i][0].dimensions.x << " " << instances[i][0].dimensions.y << " " << instances[i][0].dimensions.z << std::endl;
            for (int k = 0; k < 10; k++) {
                for (int m = 0; m < n; m++) {
                    file << instances[i][k].boxes[m].max.x << " " << instances[i][k].boxes[m].max.y << " " << instances[i][k].boxes[m].max.z << std::endl;
                }
            }
            file.close();
        }
    }

    for (int i = 0; i < 8; i++) {
        std::cout << "\n@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n";
        std::cout << "Class " << i+1 << " (" << n << " boxes)" << std::endl;
        std::cout << "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n\n";

        double fitness_sum = 0;
        int generation_sum = 0;
        long long time_sum = 0;

        for (int k = 0; k < 10; k++) {
            std::cout << "\n--------------------------------\n";
            std::cout << "Run " << std::to_string(k+1) << std::endl;
            std::cout << "----------------------------------\n\n";

            BRKGA ga(instances[i][k].dimensions, &instances[i][k].boxes);

            std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();

            Result r = ga.fit(false);
            
            std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();

            long long time_diff = std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count();
            std::cout << "Time difference = " << time_diff << "[ms]" << std::endl;
            time_sum += time_diff;
            fitness_sum += r.best_fitness;
            generation_sum += r.best_fitness_gen;
        }
        std::cout << "\nAverage best fitness: " << fitness_sum / 10 << std::endl;
        std::cout << "Average generation for best fitness: " << generation_sum / 10 << std::endl;
        std::cout << "Average runtime: " << time_sum / 10 << "[ms]" << std::endl;
    }


    return 0;
}