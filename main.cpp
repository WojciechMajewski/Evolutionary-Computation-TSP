#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <string>
#include <cstdlib>
#include <algorithm>
#include <chrono>

std::vector <std::vector <int>> read_dataset(std::string path) {
    std::ifstream file(path);

    std::vector <std::vector <int>> rows;
    std::vector <int> row;
    std::string row_string;

    while (getline(file, row_string)) {
        row.clear();
        size_t pos_start = 0, pos_end;
        std::string token;
        while ((pos_end = row_string.find(";", pos_start)) != std::string::npos) {
            token = row_string.substr(pos_start, pos_end - pos_start);
            pos_start = pos_end + 1;
            row.push_back(std::stoi(token));
        }
        row.push_back(std::stoi(row_string.substr(pos_start)));
        rows.push_back(row);
    }
    file.close();
    return rows;
} 

int euclidean_distance(int xa, int ya, int xb, int yb){
    return std::round(std::sqrt(std::pow(xa - xb, 2) + std::pow(ya - yb, 2)));
}

std::vector <std::vector <int>> create_distance_matrix(std::vector <std::vector <int>> dataset){
    std::vector <std::vector <int>> distance_matrix;
    std::vector <int> row;
    for(int i = 0; i < dataset.size(); i++){
        row.push_back(0);
    }
    for(int i = 0; i < dataset.size(); i++){
        distance_matrix.push_back(row);
    }

    for(int i = 0; i < dataset.size(); i++){
        for(int j = 0; j < i; j++){
            int distance = euclidean_distance(dataset[i][0], dataset[i][1], dataset[j][0], dataset[j][1]);
            distance_matrix[i][j] = distance;
            distance_matrix[j][i] = distance;
        }
    }

    return distance_matrix;
}

int get_cost(int start_node, int end_node, std::vector <std::vector <int>> & dataset, std::vector <std::vector <int>> & distance_matrix){
    return distance_matrix[start_node][end_node] + dataset[end_node][2];
}

int get_path_cost(std::vector <int> & path, std::vector <std::vector <int>> & dataset, std::vector <std::vector <int>> & distance_matrix){
    int cost = 0;

    for(int i = 0; i < path.size(); i++){
        cost += get_cost(path[i], path[(i+1) % path.size()], dataset, distance_matrix);
    }

    return cost;
}

std::vector <int> random_solution(std::vector <std::vector <int>> & dataset, std::vector <std::vector <int>> & distance_matrix){
    std::vector <int> path;

    std::vector <int> nodes;
    for(int i = 0; i < 200; i++){
        nodes.push_back(i);
    }

    std::random_shuffle(nodes.begin(), nodes.end());

    for(int i = 0; i < 100; i++){
        path.push_back(nodes[i]);
    }

    return path;
}

std::vector <int> nearest_solution(int starting_node, std::vector <std::vector <int>> & dataset, std::vector <std::vector <int>> & distance_matrix){
    std::vector <int> path;
    std::vector <int> available;
    int path_target_length = std::floor(dataset.size() / 2);

    path.push_back(starting_node);
    
    for(int i = 0; i < dataset.size(); i++){
        if(i != starting_node){
            available.push_back(i);
        }
    }

    for(int i = path.size(); i < path_target_length; i++){
        int shortest_distance = 1000000;
        int chosen_node_index = -1;

        for(int j = 0; j < available.size(); j++){
            int cost = get_cost(path.back(), available[j], dataset, distance_matrix);
            if(cost < shortest_distance){
                shortest_distance = cost;
                chosen_node_index = j;
            }
        }

        path.push_back(available[chosen_node_index]);
        available.erase(available.begin() + chosen_node_index);
    }

    return path;
}

std::vector <int> greedy_cycle_solution(int starting_node, std::vector <std::vector <int>> & dataset, std::vector <std::vector <int>> & distance_matrix){
    std::vector <int> path;
    std::vector <int> available;
    int path_target_length = std::floor(dataset.size() / 2);

    path.push_back(starting_node);

    for(int i = 0; i < dataset.size(); i++){
        if(i != starting_node){
            available.push_back(i);
        }
    }

    // Get second node by nearest neighbour
    int shortest_distance = 1000000;
    int chosen_node_index = -1;
    for(int j = 0; j < available.size(); j++){
        int cost = get_cost(path.back(), available[j], dataset, distance_matrix);
        // + get_cost(available[j], path.back(), dataset, distance_matrix)
        // This imitates the behavior of greedy cycle with only 1 chosen node
        if(cost < shortest_distance){
            shortest_distance = cost;
            chosen_node_index = j;
        }
    }
    path.push_back(available[chosen_node_index]);
    available.erase(available.begin() + chosen_node_index);

    // Get all subsequent nodes by greedy cycle
    for(int i = path.size(); i < path_target_length; i++){
        int shortest_extra_distance = 1000000;
        int chosen_node_index = -1;
        int chosen_path_index = -1;

        for(int j = 0; j < available.size(); j++){

            for(int k = 0; k < path.size(); k++){
                int extra_distance = get_cost(path[k], available[j], dataset, distance_matrix) + \
                                    get_cost(available[j], path[(k+1) % path.size()], dataset, distance_matrix) - \
                                    get_cost(path[k], path[(k+1) % path.size()], dataset, distance_matrix);
                if(extra_distance < shortest_extra_distance){
                    shortest_extra_distance = extra_distance;
                    chosen_node_index = j;
                    chosen_path_index = k+1;
                }
            }
        }

        path.insert(path.begin() + chosen_path_index, available[chosen_node_index]);
        available.erase(available.begin() + chosen_node_index);
    }

    return path;
}


std::vector <int> greedy_2_regret_solution(int starting_node, std::vector <std::vector <int>> & dataset, std::vector <std::vector <int>> & distance_matrix){
    std::vector <int> path;
    std::vector <int> available;
    int path_target_length = std::floor(dataset.size() / 2);

    path.push_back(starting_node);

    for(int i = 0; i < dataset.size(); i++){
        if(i != starting_node){
            available.push_back(i);
        }
    }

    int shortest_distance = 1000000;
    int chosen_node_index = -1;
    for(int j = 0; j < available.size(); j++){
        int cost = get_cost(path.back(), available[j], dataset, distance_matrix);
        if(cost < shortest_distance){
            shortest_distance = cost;
            chosen_node_index = j;
        }
    }
    path.push_back(available[chosen_node_index]);
    available.erase(available.begin() + chosen_node_index);

    for(int i = path.size(); i < path_target_length; i++){
        int highest_regret = -1000000;
        int chosen_node_index = -1;
        int chosen_path_index = -1;

        for(int j = 0; j < available.size(); j++){

            int best_location_for_this_node = -1;
            std::vector <int> top_2_scores;
            top_2_scores.push_back(1000000);
            top_2_scores.push_back(1000000);

            for(int k = 0; k < path.size(); k++){
                int extra_distance = get_cost(path[k], available[j], dataset, distance_matrix) + \
                                    get_cost(available[j], path[(k+1) % path.size()], dataset, distance_matrix) - \
                                    get_cost(path[k], path[(k+1) % path.size()], dataset, distance_matrix);

                for(int m = 0; m < 2; m++){
                    if(extra_distance < top_2_scores[m]){
                        top_2_scores.insert(top_2_scores.begin() + m, extra_distance);
                        top_2_scores.pop_back();
                        if(m == 0){
                            best_location_for_this_node = k+1;
                        }
                        break;
                    }
                }
            }

            int regret_2 = top_2_scores[1] - top_2_scores[0];

            if(regret_2 > highest_regret){
                highest_regret = regret_2;
                chosen_node_index = j;
                chosen_path_index = best_location_for_this_node;
            }
        }

        path.insert(path.begin() + chosen_path_index, available[chosen_node_index]);
        available.erase(available.begin() + chosen_node_index);
    }

    return path;
}


std::vector <int> greedy_weighted_solution(int starting_node, std::vector <std::vector <int>> & dataset, std::vector <std::vector <int>> & distance_matrix){
    std::vector <int> path;
    std::vector <int> available;
    int path_target_length = std::floor(dataset.size() / 2);

    path.push_back(starting_node);

    for(int i = 0; i < dataset.size(); i++){
        if(i != starting_node){
            available.push_back(i);
        }
    }

    int shortest_distance = 1000000;
    int chosen_node_index = -1;
    for(int j = 0; j < available.size(); j++){
        int cost = get_cost(path.back(), available[j], dataset, distance_matrix);
        if(cost < shortest_distance){
            shortest_distance = cost;
            chosen_node_index = j;
        }
    }
    path.push_back(available[chosen_node_index]);
    available.erase(available.begin() + chosen_node_index);

    for(int i = path.size(); i < path_target_length; i++){
        int highest_regret = -1000000;
        int chosen_node_index = -1;
        int chosen_path_index = -1;

        for(int j = 0; j < available.size(); j++){

            int best_location_for_this_node = -1;
            std::vector <int> top_2_scores;
            top_2_scores.push_back(1000000);
            top_2_scores.push_back(1000000);

            for(int k = 0; k < path.size(); k++){
                int extra_distance = get_cost(path[k], available[j], dataset, distance_matrix) + \
                                    get_cost(available[j], path[(k+1) % path.size()], dataset, distance_matrix) - \
                                    get_cost(path[k], path[(k+1) % path.size()], dataset, distance_matrix);

                for(int m = 0; m < 2; m++){
                    if(extra_distance < top_2_scores[m]){
                        top_2_scores.insert(top_2_scores.begin() + m, extra_distance);
                        top_2_scores.pop_back();
                        if(m == 0){
                            best_location_for_this_node = k+1;
                        }
                        break;
                    }
                }
            }

            int regret_2 = top_2_scores[1] - top_2_scores[0];
            // The only change from 2-regret:
            regret_2 -= top_2_scores[0];

            if(regret_2 > highest_regret){
                highest_regret = regret_2;
                chosen_node_index = j;
                chosen_path_index = best_location_for_this_node;
            }
        }

        path.insert(path.begin() + chosen_path_index, available[chosen_node_index]);
        available.erase(available.begin() + chosen_node_index);
    }

    return path;
}


std::vector <std::vector <int>> calculate_best_paths(std::vector <std::vector <int>> & dataset, std::vector <std::vector <int>> & distance_matrix, std::string filename = "", std::string dataset_name = "example"){
    int iterations = 200;
    std::vector <std::vector <int>> best_paths;
    std::vector <int> best_scores;
    std::vector <int> worst_scores;
    std::vector <float> average_scores;
    std::vector <int> times;
    
    std::vector <std::string> algorithm_names;
    algorithm_names.push_back("Nearest");
    algorithm_names.push_back("Greedy Cycle");
    algorithm_names.push_back("Greedy 2-Regret");
    algorithm_names.push_back("Greedy Weighted");

    best_paths.push_back(nearest_solution(0, dataset, distance_matrix));
    best_paths.push_back(greedy_cycle_solution(0, dataset, distance_matrix));
    best_paths.push_back(greedy_2_regret_solution(0, dataset, distance_matrix));
    best_paths.push_back(greedy_weighted_solution(0, dataset, distance_matrix));

    for(int i = 0; i < best_paths.size(); i++){
        best_scores.push_back(get_path_cost(best_paths[i], dataset, distance_matrix));
        worst_scores.push_back(get_path_cost(best_paths[i], dataset, distance_matrix));
        average_scores.push_back(get_path_cost(best_paths[i], dataset, distance_matrix));
    }

    for(int j = 0; j < algorithm_names.size(); j++){
        std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();

        for(int i = 1; i < iterations; i++){
        //if(i%20 == 0) std::cout << "iteration " << i << "\n";
            std::vector <int> solution;
            int cost;

            if(j == 0){
                solution = nearest_solution(i, dataset, distance_matrix);
            }
            if(j == 1){
                solution = greedy_cycle_solution(i, dataset, distance_matrix);
            }
            if(j == 2){
                solution = greedy_2_regret_solution(i, dataset, distance_matrix);
            }
            if(j == 3){
                solution = greedy_weighted_solution(i, dataset, distance_matrix);
            }

            cost = get_path_cost(solution, dataset, distance_matrix);
            average_scores[j] += cost;
            if(cost < best_scores[j]){
                best_scores[j] = cost;
                best_paths[j] = solution;
            }
            if(cost > worst_scores[j]){
                worst_scores[j] = cost;
            }
        }
        average_scores[j] /= iterations;
        std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
        times.push_back(std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count());

        if(filename == ""){
            std::cout << "Best score:\n" << best_scores[j] << "\n";
            std::cout << "Worst score:\n" << worst_scores[j] << "\n";
            std::cout << "Average score:\n" << average_scores[j] << "\n";
            std::cout << "Time difference:\n" << std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count() << "ms" << std::endl;
        }
    }

    
    if(filename != ""){
        std::ofstream ofs;
        ofs.open(filename, std::ios_base::app);
        ofs << "\nDATASET " << dataset_name << "\n";
        for(int j = 0; j < algorithm_names.size(); j++){
            ofs << "\n" << algorithm_names[j] << "\n";
            ofs << "Best score:\n" << best_scores[j] << "\n";
            for(int k = 0; k < best_paths[j].size() - 1; k++){
                ofs << best_paths[j][k] << ", ";
            }
            ofs << best_paths[j][best_paths.size() - 1] << "\n";
            ofs << "Worst score:\n" << worst_scores[j] << "\n";
            ofs << "Average score:\n" << average_scores[j] << "\n";
            ofs << "Time:\n" << times[j] << "ms\n";
        }
        ofs.close();
    }
    return best_paths;
}


int main(){
    std::srand(148253);

    std::string filename = "cpp_results.txt";


    if(filename != ""){
        std::ofstream ofs;
        ofs.open(filename, std::ofstream::out | std::ofstream::trunc);
        ofs.close();
    }

    std::vector <std::string> dataset_paths;
    dataset_paths.push_back("Data/TSPA.csv");
    dataset_paths.push_back("Data/TSPB.csv");
    dataset_paths.push_back("Data/TSPC.csv");
    dataset_paths.push_back("Data/TSPD.csv");
    
    std::vector <std::string> algorithm_names;
    algorithm_names.push_back("Nearest");
    algorithm_names.push_back("Greedy Cycle");
    algorithm_names.push_back("Greedy 2-Regret");
    algorithm_names.push_back("Greedy Weighted");


    for(int i = 0; i < dataset_paths.size(); i++){
        std::vector <std::vector <int>> dataset = read_dataset(dataset_paths[i]);
        std::vector <std::vector <int>> distance_matrix = create_distance_matrix(dataset);
        std::string dataset_name;
        dataset_name += (char(65 + i));
        std::vector <std::vector <int>> best_paths = calculate_best_paths(dataset, distance_matrix, filename, dataset_name);

        
    }
    return 0;
}


