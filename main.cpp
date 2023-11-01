#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <string>
#include <cstdlib>
#include <algorithm>

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


int main(){
    std::srand(148253);

    std::vector <std::vector <int>> dataset_A = read_dataset("Data/TSPA.csv");
    std::vector <std::vector <int>> distance_matrix_A = create_distance_matrix(dataset_A);

    std::vector <int> solution_random = random_solution(dataset_A, distance_matrix_A);
    std::cout << get_path_cost(solution_random, dataset_A, distance_matrix_A) << "\n";

    std::vector <int> solution_nearest = nearest_solution(0, dataset_A, distance_matrix_A);
    std::cout << get_path_cost(solution_nearest, dataset_A, distance_matrix_A) << "\n";

    std::vector <int> solution_greedy_cycle = greedy_cycle_solution(0, dataset_A, distance_matrix_A);
    std::cout << get_path_cost(solution_greedy_cycle, dataset_A, distance_matrix_A) << "\n";

    std::vector <int> solution_greedy_2_regret = greedy_2_regret_solution(0, dataset_A, distance_matrix_A);
    std::cout << get_path_cost(solution_greedy_2_regret, dataset_A, distance_matrix_A) << "\n";

    std::vector <int> solution_greedy_weighted = greedy_weighted_solution(0, dataset_A, distance_matrix_A);
    std::cout << get_path_cost(solution_greedy_weighted, dataset_A, distance_matrix_A) << "\n";
    //for(int i = 0; i < distance_matrix_A.size(); i++){
    //    for(int j = 0; j < distance_matrix_A[i].size(); j++){
    //        std::cout << distance_matrix_A[i][j] << "\n";
    //    }
    //}

    // Create extra cost list

    // Greedy heuristic algorithms


    return 0;
}


