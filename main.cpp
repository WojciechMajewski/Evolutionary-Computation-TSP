#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <string>
#include <cstdlib>
#include <algorithm>
#include <chrono>
#include <algorithm>
#include <queue>
#include <numeric>
#include <random>

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

    std::random_device rd;
    std::mt19937 g(rd());

    std::shuffle(nodes.begin(), nodes.end(), g);

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

    // Get second node by nearest neighbor
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

std::vector <int> local_search(bool steepest_neighborhood, bool edges_exchange, std::vector <int> solution, std::vector <std::vector <int>> & dataset, std::vector <std::vector <int>> & distance_matrix){
    // if not steepest_neighborhood, then greedy
    // if not edges_exchange, then nodes exchange // But intra- and inter- should both be used
    // if starting_node == -1 then random starting solution

    // sort the starting solution, then for each of 200 check if in sorted, if yes advance sorted, if no add to available
    std::vector <int> available_nodes;
    std::vector <int> unavailable_nodes = solution;
    std::sort(unavailable_nodes.begin(), unavailable_nodes.end());

    //std::cout << "\n\nLOCAL SEARCH  \n";
    int index = 0;
    for(int i = 0; i < dataset.size(); i++){
        if(index < unavailable_nodes.size() && unavailable_nodes[index] == i){
            index++;
        }
        else{
            available_nodes.push_back(i);
        }
    }

    // 100 iterations instead of while(true)
    for(int t = 0; t < 1000; t++){
        int best_improvement = 0;
        int new_node_index = -1;
        int replacing_node_index = -1;

        int first_edge_start_index = -1; // By definition the end node will be (i+1 mod size)
        int second_edge_start_index = -1;

        int first_to_switch_index = -1;
        int second_to_switch_index = -1;

        if(steepest_neighborhood){
            // First new nodes
            for(int i = 0; i < solution.size(); i++){
                int old_cost = 0;
                old_cost += get_cost(solution[(i-1 + solution.size()) % solution.size()], solution[i], dataset, distance_matrix);
                old_cost += get_cost(solution[i], solution[(i+1) % solution.size()], dataset, distance_matrix);
                for(int j = 0; j < available_nodes.size(); j++){
                    int cost_improvement = 0;
                    cost_improvement -= get_cost(solution[(i-1 + solution.size()) % solution.size()], available_nodes[j], dataset, distance_matrix);
                    cost_improvement -= get_cost(available_nodes[j], solution[(i+1) % solution.size()], dataset, distance_matrix);
                    cost_improvement += old_cost;
                    if(cost_improvement > best_improvement){
                        best_improvement = cost_improvement;
                        replacing_node_index = i;
                        new_node_index = j;
                    }
                }
            }

            // Then edge/node reordering
            if(edges_exchange){
                // Two pairs of consecutive nodes
                for(int i = 0; i < solution.size(); i++){
                    for(int j = 0; (j+1) % solution.size() < i; j++){
                        if(j == i || j == (i+1) % solution.size() || (j+1) % solution.size() == i) continue; // catch intersections
                        int old_cost = 0;
                        int cost_improvement = 0;

                        old_cost += get_cost(solution[i], solution[(i+1) % solution.size()], dataset, distance_matrix);
                        old_cost += get_cost(solution[j], solution[(j+1) % solution.size()], dataset, distance_matrix);

                        old_cost += dataset[solution[i]][2];
                        //old_cost += dataset[solution[j]][2];

                        //cost_improvement -= get_cost(solution[i], solution[(j+1) % solution.size()], dataset, distance_matrix);
                        //cost_improvement -= get_cost(solution[j], solution[(i+1) % solution.size()], dataset, distance_matrix);

                        cost_improvement -= get_cost(solution[(j+1) % solution.size()], solution[(i+1) % solution.size()], dataset, distance_matrix);
                        cost_improvement -= get_cost(solution[j], solution[i], dataset, distance_matrix);

                        cost_improvement -= dataset[solution[(j+1) % solution.size()]][2];
                        //cost_improvement -= dataset[solution[j]][2];

                        cost_improvement += old_cost;

                        if(cost_improvement > best_improvement){
                            best_improvement = cost_improvement;

                            first_edge_start_index = j;
                            second_edge_start_index = i; // As j is always smaller than i

                            new_node_index = -1;
                        }
                    }
                }

            }
            else{ // nodes exchange
                for(int i = 0; i < solution.size(); i++){
                    for(int j = 0; j < i; j++){
                        int old_cost = 0;
                        int cost_improvement = 0;
                        if(i == solution.size() - 1 && j == 0){ // Then i is just before j as the path is a cycle
                            //continue;
                            // TODO: FIX
                            old_cost += get_cost(solution[(i-1 + solution.size()) % solution.size()], solution[i], dataset, distance_matrix);
                            old_cost += get_cost(solution[i], solution[j], dataset, distance_matrix);
                            old_cost += get_cost(solution[j], solution[(j+1) % solution.size()], dataset, distance_matrix);

                            cost_improvement -= get_cost(solution[(i-1 + solution.size()) % solution.size()], solution[j], dataset, distance_matrix);
                            cost_improvement -= get_cost(solution[j], solution[i], dataset, distance_matrix);
                            cost_improvement -= get_cost(solution[i], solution[(j+1) % solution.size()], dataset, distance_matrix);
                        }
                        else{
                            if(j == i-1){
                            //continue;
                            old_cost += get_cost(solution[(j-1 + solution.size()) % solution.size()], solution[j], dataset, distance_matrix);
                            old_cost += get_cost(solution[j], solution[i], dataset, distance_matrix);
                            old_cost += get_cost(solution[i], solution[(i+1) % solution.size()], dataset, distance_matrix);

                            cost_improvement -= get_cost(solution[(j-1 + solution.size()) % solution.size()], solution[i], dataset, distance_matrix);
                            cost_improvement -= get_cost(solution[i], solution[j], dataset, distance_matrix);
                            cost_improvement -= get_cost(solution[j], solution[(i+1) % solution.size()], dataset, distance_matrix);  
                            }
                            else{
                                //continue;
                                old_cost += get_cost(solution[(i-1 + solution.size()) % solution.size()], solution[i], dataset, distance_matrix);
                                old_cost += get_cost(solution[i], solution[(i+1) % solution.size()], dataset, distance_matrix);
                                old_cost += get_cost(solution[(j-1 + solution.size()) % solution.size()], solution[j], dataset, distance_matrix);
                                old_cost += get_cost(solution[j], solution[(j+1) % solution.size()], dataset, distance_matrix);

                                cost_improvement -= get_cost(solution[(i-1 + solution.size()) % solution.size()], solution[j], dataset, distance_matrix);
                                cost_improvement -= get_cost(solution[j], solution[(i+1) % solution.size()], dataset, distance_matrix);
                                cost_improvement -= get_cost(solution[(j-1 + solution.size()) % solution.size()], solution[i], dataset, distance_matrix);
                                cost_improvement -= get_cost(solution[i], solution[(j+1) % solution.size()], dataset, distance_matrix);
                            }
                        }

                        cost_improvement += old_cost;

                        if(cost_improvement > best_improvement){
                            best_improvement = cost_improvement;
                            first_to_switch_index = i;
                            second_to_switch_index = j;
                            new_node_index = -1;
                        }

                    }

                }
            }

            if(best_improvement > 0){
                if(new_node_index != -1){ // Then new node is best
                
                    //std::cout << "[imp " << best_improvement << " node " << solution[replacing_node_index] << " into " << available_nodes[new_node_index] << " in " << replacing_node_index << "], ";

                    int temp = solution[replacing_node_index];
                    solution[replacing_node_index] = available_nodes[new_node_index];
                    available_nodes[new_node_index] = temp;
                    
                }
                else{
                    if(edges_exchange){ // Then edge exchange worked
                        if(((first_edge_start_index + 1) % solution.size()) > ((second_edge_start_index + 1) % solution.size())){
                            int temp = first_edge_start_index;
                            first_edge_start_index = second_edge_start_index;
                            second_edge_start_index = temp;
                        }
                        // Edge exchange through flipping a subpath
                        std::reverse(solution.begin() + ((first_edge_start_index + 1) % solution.size()), solution.begin() + ((second_edge_start_index + 1) % solution.size()));
                        
                        //std::cout << "[imp " << best_improvement << " edge " << solution[first_edge_start_index] << "->" << solution[(first_edge_start_index + 1) % solution.size()];
                        //std::cout << " into " << solution[second_edge_start_index] << "->" << solution[(second_edge_start_index + 1) % solution.size()];
                        //std::cout << " in " << first_edge_start_index << " and " << second_edge_start_index << "], ";
                    }
                    else{ // Then node exchange worked
                    
                        int temp = solution[first_to_switch_index];
                        solution[first_to_switch_index] = solution[second_to_switch_index];
                        solution[second_to_switch_index] = temp;
                    }
                }
            }
            else{ // No improvement, break the local search loop, local optimum found!
                return solution;
            }
        }
        else{
            int random_offset = rand() % solution.size(); //TODO FIX
            //int io = (i + random_offset) % solution.size(); // that's for greedy

            int i_old_node_iterating = 0;
            int i_old_node = (i_old_node_iterating + random_offset) % solution.size();
            int j_new_node = 0;
            int old_cost_new_node = 0;
            bool new_node_finished = false;

            int i_intra_iter = 0; // This calculates when the loop should end
            int i_intra = (i_intra_iter + random_offset) % solution.size(); // This is used as an index in calculations, and is always offset from i_intra_iter by set random amount
            int j_intra = 0;
            bool intra_finished = false;

            // Similar to steepest, but it needs to randomize whether it tries to add a new node or do an internal swap
            while(true){
                int choice = rand() % 2; // TODO FIX

                // New node advance
                if(choice == 0 && !new_node_finished){ 
                    j_new_node++;
                    if(j_new_node >= available_nodes.size()){
                        j_new_node = 0;
                        i_old_node_iterating++;
                        if(i_old_node_iterating >= solution.size()){
                            new_node_finished = true;
                            continue; // No node replacement will improve the score
                        }

                        i_old_node = (i_old_node_iterating + random_offset) % solution.size();

                        old_cost_new_node = 0;
                        old_cost_new_node += get_cost(solution[(i_old_node-1 + solution.size()) % solution.size()], solution[i_old_node], dataset, distance_matrix);
                        old_cost_new_node += get_cost(solution[i_old_node], solution[(i_old_node+1) % solution.size()], dataset, distance_matrix);
                    }
                    int cost_improvement = old_cost_new_node;
                    cost_improvement -= get_cost(solution[(i_old_node-1 + solution.size()) % solution.size()], available_nodes[j_new_node], dataset, distance_matrix); 
                    cost_improvement -= get_cost(available_nodes[j_new_node], solution[(i_old_node+1) % solution.size()], dataset, distance_matrix);
                    if(cost_improvement > 0){
                        int temp = solution[i_old_node];
                        solution[i_old_node] = available_nodes[j_new_node];
                        available_nodes[j_new_node] = temp;
                        break; // Found improvement, end iteration, go find the next
                    }
                }
                else if (choice != 0 && !intra_finished){
                    // Edge rearrangement advance
                    if(edges_exchange){
                        
                        j_intra++;
                        if(j_intra+1 >= i_intra){
                            j_intra = 0;
                            i_intra_iter++;
                            if(i_intra_iter >= solution.size()){
                                intra_finished = true;
                                continue; // No edge rearrangement will improve the score
                            }

                            i_intra = (i_intra_iter + random_offset) % solution.size();

                        }

                        if(j_intra == i_intra || j_intra == (i_intra+1) % solution.size() || (j_intra+1) % solution.size() == i_intra) continue;
                        int old_cost = 0;
                        int cost_improvement = 0;

                        old_cost += get_cost(solution[i_intra], solution[(i_intra+1) % solution.size()], dataset, distance_matrix);
                        old_cost += get_cost(solution[j_intra], solution[(j_intra+1) % solution.size()], dataset, distance_matrix);

                        old_cost += dataset[solution[i_intra]][2];
                        //old_cost += dataset[solution[j_intra]][2];

                        cost_improvement -= get_cost(solution[(j_intra+1) % solution.size()], solution[(i_intra+1) % solution.size()], dataset, distance_matrix);
                        cost_improvement -= get_cost(solution[j_intra], solution[i_intra], dataset, distance_matrix);

                        cost_improvement -= dataset[solution[(j_intra+1) % solution.size()]][2];
                        //cost_improvement -= dataset[solution[j_intra]][2];

                        cost_improvement += old_cost;


                        if(cost_improvement > 0){
                            std::reverse(solution.begin() + ((j_intra + 1) % solution.size()), solution.begin() + ((i_intra + 1) % solution.size()));
                            break;
                        }
                        
                    }

                    // Node rearrangement advance
                    else{
                        j_intra++;
                        if(j_intra >= i_intra){
                            j_intra = 0;
                            i_intra_iter++;
                            if(i_intra_iter >= solution.size()){
                                intra_finished = true;
                                continue; // No edge rearrangement will improve the score
                            }

                            i_intra = (i_intra_iter + random_offset) % solution.size();
                        }

                        int old_cost = 0;
                        int cost_improvement = 0;
                        if(i_intra == solution.size() - 1 && j_intra == 0){ // Then i is just before j as the path is a cycle
                            old_cost += get_cost(solution[(i_intra-1 + solution.size()) % solution.size()], solution[i_intra], dataset, distance_matrix);
                            old_cost += get_cost(solution[i_intra], solution[j_intra], dataset, distance_matrix);
                            old_cost += get_cost(solution[j_intra], solution[(j_intra+1) % solution.size()], dataset, distance_matrix);
                            cost_improvement -= get_cost(solution[(i_intra-1 + solution.size()) % solution.size()], solution[j_intra], dataset, distance_matrix);
                            cost_improvement -= get_cost(solution[j_intra], solution[i_intra], dataset, distance_matrix);
                            cost_improvement -= get_cost(solution[i_intra], solution[(j_intra+1) % solution.size()], dataset, distance_matrix);
                        }
                        else if(j_intra == i_intra-1){
                            old_cost += get_cost(solution[(j_intra-1 + solution.size()) % solution.size()], solution[j_intra], dataset, distance_matrix);
                            old_cost += get_cost(solution[j_intra], solution[i_intra], dataset, distance_matrix);
                            old_cost += get_cost(solution[i_intra], solution[(i_intra+1) % solution.size()], dataset, distance_matrix);
                            cost_improvement -= get_cost(solution[(j_intra-1 + solution.size()) % solution.size()], solution[i_intra], dataset, distance_matrix);
                            cost_improvement -= get_cost(solution[i_intra], solution[j_intra], dataset, distance_matrix);
                            cost_improvement -= get_cost(solution[j_intra], solution[(i_intra+1) % solution.size()], dataset, distance_matrix);
                        }
                        else{
                            old_cost += get_cost(solution[(i_intra-1 + solution.size()) % solution.size()], solution[i_intra], dataset, distance_matrix);
                            old_cost += get_cost(solution[i_intra], solution[(i_intra+1) % solution.size()], dataset, distance_matrix);
                            old_cost += get_cost(solution[(j_intra-1 + solution.size()) % solution.size()], solution[j_intra], dataset, distance_matrix);
                            old_cost += get_cost(solution[j_intra], solution[(j_intra+1) % solution.size()], dataset, distance_matrix);
                            cost_improvement -= get_cost(solution[(i_intra-1 + solution.size()) % solution.size()], solution[j_intra], dataset, distance_matrix);
                            cost_improvement -= get_cost(solution[j_intra], solution[(i_intra+1) % solution.size()], dataset, distance_matrix);
                            cost_improvement -= get_cost(solution[(j_intra-1 + solution.size()) % solution.size()], solution[i_intra], dataset, distance_matrix);
                            cost_improvement -= get_cost(solution[i_intra], solution[(j_intra+1) % solution.size()], dataset, distance_matrix);
                        }

                        cost_improvement += old_cost;

                        if(cost_improvement > best_improvement){
                            int temp = solution[i_intra];
                            solution[i_intra] = solution[j_intra];
                            solution[j_intra] = temp;
                            break;
                        }
                    }
                }
                else{
                    
                    return solution;
                    // No improvement can be found by greedy, so end
                }
            }
        }
    }

    return solution;
}

std::vector <int> local_candidate_moves(std::vector <int> solution, std::vector <std::vector <std::tuple <int, int>>> closest_nodes, std::vector <std::vector <int>> & dataset, std::vector <std::vector <int>> & distance_matrix){
    // sort the starting solution, then for each of 200 check if in sorted, if yes advance sorted, if no add to available
    std::vector <int> available_nodes;
    std::vector <int> unavailable_nodes = solution;
    std::sort(unavailable_nodes.begin(), unavailable_nodes.end());

    // Creating a set of nodes that are not in the chosen path but could be placed there
    int index = 0;
    for(int i = 0; i < dataset.size(); i++){
        if(index < unavailable_nodes.size() && unavailable_nodes[index] == i){
            index++;
        }
        else{
            available_nodes.push_back(i);
        }
    }

    for(int t = 0; t < 1000; t++){
        int best_improvement = 0;
        int new_node_index = -1;
        int replaced_node_index = -1;

        int first_edge_start_index = -1; // By definition edge end will be ((edge_start + 1) % size)
        int second_edge_start_index = -1;

        int first_to_switch_index = -1;
        int second_to_switch_index = -1;

        // Consider replacing nodes
        
        // Using candidate moves:
        // For each node in solution, check if the next node could be replaced by any one of 10 best nodes.
        // By adding distance to one of the 10 plus distance from it to the next node in solution
        
        // This is only calculated for nodes inside the solution
        for(int i = 0; i < solution.size(); i++){
            for(int j = 0; j < 10; j++){
                int neighbor_node = std::get<0>(closest_nodes[solution[i]][j]);
                int neighbor_distance = std::get<1>(closest_nodes[solution[i]][j]);

                int next_node = solution[(i+1) % solution.size()];
                int further_next_node = solution[(i+2) % solution.size()];
                
                if(next_node == neighbor_node){
                    continue;
                }
                

                bool neighbor_inside = false;
                int neighbor_index = 0;
                for(int k = 0; k < solution.size(); k++){
                    if(solution[k] == neighbor_node){
                        neighbor_inside = true;
                        neighbor_index = k;
                        break;
                    }
                }
                if(!neighbor_inside){
                    for(int k = 0; k < available_nodes.size(); k++){
                        if(available_nodes[k] == neighbor_node){
                            neighbor_index = k;
                            break;
                        }
                    }
                }

                // If neighbor_node on the path, consider an edge exchange
                if(neighbor_inside){
                    // Four nodes: node, next node and neighbor, next neighbor
                    // Ends with: node, neighbor and next node, next neighbor

                    // All need to be different!
                    
                    // If the neighbor is the previous node, edges can't be exchanged
                    if((neighbor_index + 1) % solution.size() == i){
                        continue;
                    }

                    int old_cost = 0;
                    old_cost += get_cost(solution[i], next_node, dataset, distance_matrix);
                    old_cost += get_cost(neighbor_node, solution[(neighbor_index + 1) % solution.size()], dataset, distance_matrix);
                    old_cost += dataset[neighbor_node][2];

                    int cost_improvement = 0;
                    cost_improvement -= neighbor_distance;
                    cost_improvement -= get_cost(next_node, solution[(neighbor_index + 1) % solution.size()], dataset, distance_matrix);
                    cost_improvement -= dataset[next_node][2];
                    
                    cost_improvement += old_cost;

                    if(cost_improvement > best_improvement){
                        best_improvement = cost_improvement;

                        first_edge_start_index = neighbor_index;
                        second_edge_start_index = i; // As j is always smaller than i

                        //new_node_index = -1;
                    }

                }
                // If neighbor_node not on the path, consider replacing next_node with neighbor node (node replacement)
                else{
                    int old_cost = 0;
                    old_cost += get_cost(solution[i], next_node, dataset, distance_matrix);
                    old_cost += get_cost(next_node, further_next_node, dataset, distance_matrix);

                    int cost_improvement = 0;
                    cost_improvement -= neighbor_distance;
                    cost_improvement -= get_cost(neighbor_node, further_next_node, dataset, distance_matrix);

                    cost_improvement += old_cost;

                    if(cost_improvement > best_improvement){
                        best_improvement = cost_improvement;
                        replaced_node_index = (i+1) % solution.size();
                        new_node_index = neighbor_index; // index from available nodes;
                    }
                }
            }
        }
        
        
        for(int i = 0; i < available_nodes.size(); i++){
            for(int j = 0; j < 10; j++){
                int neighbor_node = std::get<0>(closest_nodes[available_nodes[i]][j]);
                int neighbor_distance = std::get<1>(closest_nodes[available_nodes[i]][j]);
                
                int neighbor_index = -1;
                for(int k = 0; k < solution.size(); k++){
                    if(solution[k] == neighbor_node){
                        neighbor_index = k;
                        break;
                    }
                }

                if(neighbor_index == -1){ // Both the node and the neighbor are outside the path
                    continue;
                }

                int neighbor_prev_index = (neighbor_index - 1 + solution.size()) % solution.size();
                int neighbor_more_prev_index = (neighbor_index - 2 + solution.size()) % solution.size();

                int old_cost = 0;
                old_cost += get_cost(solution[neighbor_more_prev_index], solution[neighbor_prev_index], dataset, distance_matrix);
                old_cost += get_cost(solution[neighbor_prev_index], solution[neighbor_index], dataset, distance_matrix);

                int cost_improvement = 0;
                cost_improvement -= neighbor_distance;
                cost_improvement -= get_cost(solution[neighbor_more_prev_index], available_nodes[i], dataset, distance_matrix);

                cost_improvement += old_cost;

                if(cost_improvement > best_improvement){
                    best_improvement = cost_improvement;
                    replaced_node_index = neighbor_prev_index;
                    new_node_index = i; // index from available nodes;
                }

                
            }
        }
        

        if(best_improvement > 0){
            // New node is best
            if(new_node_index != -1){
                int temp = solution[replaced_node_index];
                solution[replaced_node_index] = available_nodes[new_node_index];
                available_nodes[new_node_index] = temp;
            }
            // Then edge exchange is best
            else{
                if(((first_edge_start_index + 1) % solution.size()) > ((second_edge_start_index + 1) % solution.size())){
                    int temp = first_edge_start_index;
                    first_edge_start_index = second_edge_start_index;
                    second_edge_start_index = temp;
                }

                std::reverse(solution.begin() + ((first_edge_start_index + 1) % solution.size()), solution.begin() + ((second_edge_start_index + 1) % solution.size()));
            }
        }
        else{ // No improvement, break the local search loop, local optimum found!
            return solution;
        }

    }

    return solution;
}

// Descending order
bool sort_by_first(const std::vector<int>& v1, const std::vector<int>& v2){
    return v1[0] > v2[0];
}

int get_index_to_insert_descending(const std::vector<std::vector<int>>& array, std::vector<int>& item){
    int compared_index = 0;
    int insert_index = array.size();

    if (array.size() == 0 || (array[array.size() - 1][compared_index] > item[compared_index]))
        return insert_index;
    
    int low = 0;
    int high = array.size() - 1;
    while (low <= high){
        int mid = low + (high - low) / 2;
        if (array[mid][compared_index] > item[compared_index]){
            low = mid + 1;
        }
        else{
            high = mid - 1;
        }
    }
    insert_index = high + 1;
    return insert_index;
}


std::vector <int> local_delta_working(std::vector <int> solution, std::vector <std::vector <int>> & dataset, std::vector <std::vector <int>> & distance_matrix){
    std::vector <int> available_nodes;
    std::vector <int> unavailable_nodes = solution;
    std::sort(unavailable_nodes.begin(), unavailable_nodes.end());

    //std::cout << "\n\nDELTA \n";
    int index = 0;
    for(int i = 0; i < dataset.size(); i++){
        if(index < unavailable_nodes.size() && unavailable_nodes[index] == i){
            index++;
        }
        else{
            available_nodes.push_back(i);
        }
    }
    
    std::vector <std::vector <int>> improvement_list; // <improvement, (0 - node replacement / 1 - edge exchange), (cycle node, available node) / (edge1 start, edge1 end, edge2 start, edge2 end)>
    // <improvement, 0 (node replacement), node, node before, node after, replacement node>
    // <improvement, 1 (edge exchange), edge1 start node, edge1 end node, edge2 start node, edge2 end node>

    // Initiate improvement list
    for(int i = 0; i < solution.size(); i++){
        int old_cost = 0;
        old_cost += get_cost(solution[(i-1 + solution.size()) % solution.size()], solution[i], dataset, distance_matrix);
        old_cost += get_cost(solution[i], solution[(i+1) % solution.size()], dataset, distance_matrix);
        for(int j = 0; j < available_nodes.size(); j++){

            int cost_improvement = 0;
            cost_improvement -= get_cost(solution[(i-1 + solution.size()) % solution.size()], available_nodes[j], dataset, distance_matrix);
            cost_improvement -= get_cost(available_nodes[j], solution[(i+1) % solution.size()], dataset, distance_matrix);
            cost_improvement += old_cost;
            if(cost_improvement > 0){
                std::vector <int> move{cost_improvement, 0, solution[i], solution[(i-1 + solution.size()) % solution.size()], solution[(i+1) % solution.size()], available_nodes[j]};
                improvement_list.push_back(move);
            }
        }
    }
    
    // Edges
    for(int i = 0; i < solution.size(); i++){
        for(int j = 0; (j+1) % solution.size() < i; j++){
            if(j == i || j == (i+1) % solution.size() || (j+1) % solution.size() == i) continue; // catch intersections
            int old_cost = 0;
            int cost_improvement = 0;

            old_cost += get_cost(solution[i], solution[(i+1) % solution.size()], dataset, distance_matrix);
            old_cost += get_cost(solution[j], solution[(j+1) % solution.size()], dataset, distance_matrix);
            old_cost += dataset[solution[i]][2];

            cost_improvement -= get_cost(solution[(j+1) % solution.size()], solution[(i+1) % solution.size()], dataset, distance_matrix);
            cost_improvement -= get_cost(solution[j], solution[i], dataset, distance_matrix);
            cost_improvement -= dataset[solution[(j+1) % solution.size()]][2];

            cost_improvement += old_cost;

            if(cost_improvement > 0){
                std::vector <int> move{cost_improvement, 1, solution[i], solution[(i+1) % solution.size()], solution[j], solution[(j+1) % solution.size()]};
                improvement_list.push_back(move);
            }
        }
    }

    ///*
    // Reversed edges (one of them)
    for(int i = 0; i < solution.size(); i++){
        for(int j = 0; (j+1) % solution.size() < i; j++){
            if(j == i || j == (i+1) % solution.size() || (j+1) % solution.size() == i) continue; // catch intersections
            int old_cost = 0;
            int cost_improvement = 0;

            old_cost += get_cost(solution[i], solution[(i+1) % solution.size()], dataset, distance_matrix);
            old_cost += get_cost(solution[(j+1) % solution.size()], solution[j], dataset, distance_matrix); // Switched
            old_cost += dataset[solution[i]][2];

            cost_improvement -= get_cost(solution[j], solution[(i+1) % solution.size()], dataset, distance_matrix); // j switched
            cost_improvement -= get_cost(solution[(j+1) % solution.size()], solution[i], dataset, distance_matrix); // j switched
            cost_improvement -= dataset[solution[j]][2]; // j switched

            cost_improvement += old_cost;

            if(cost_improvement > 0){
                std::vector <int> move{cost_improvement, 1, solution[i], solution[(i+1) % solution.size()], solution[(j+1) % solution.size()], solution[j]};
                improvement_list.push_back(move);
            }
        }
    }
    //*/

    // Sort the available moves
    std::sort(improvement_list.begin(), improvement_list.end(), sort_by_first);

    // 100 iterations instead of while(true)
    for(int t = 0; t < 1000; t++){
        // Evaluate new moves
        bool found_improvement = false;
        std::vector <std::vector <int>> to_insert;

        for(int index = 0; index < improvement_list.size(); index++){

            if(improvement_list[index][1] == 0){
                
                int node_index = 0;
                for(;node_index < solution.size() && solution[node_index] != improvement_list[index][2]; node_index++);
                if(node_index == solution.size()){ // No longer in solution
                    continue;
                }

                // <improvement, 0 (node replacement), node, node before, node after, replacement node>

                //if(solution[improvement_list[index][2]] != improvement_list[index][3]) continue; // Node no longer there

                if(!((solution[(node_index-1 + solution.size()) % solution.size()] == improvement_list[index][4] &&
                    solution[(node_index+1) % solution.size()] == improvement_list[index][3]) || 
                    (solution[(node_index-1 + solution.size()) % solution.size()] == improvement_list[index][3] &&
                    solution[(node_index+1) % solution.size()] == improvement_list[index][4]))){ // If no longer 3 nodes in a row (either way)
                    continue;
                }
                
                int replacement_index = 0;
                for(;replacement_index < available_nodes.size() && available_nodes[replacement_index] != improvement_list[index][5]; replacement_index++);
                if(replacement_index == available_nodes.size()){ // No longer in available
                    continue;
                }
                //if(available_nodes[improvement_list[index][6]] != improvement_list[index][7]) continue;

                // FOUND legal node replacement
                //std::cout << "[imp " << improvement_list[index][0] << " node " << solution[node_index] << " into " << available_nodes[replacement_index] << " in " << node_index << "], ";
                
                found_improvement = true;
                // update the solution

                int temp = solution[node_index];
                solution[node_index] = available_nodes[replacement_index];
                available_nodes[replacement_index] = temp;

                // add new moves for the replaced node, node before and after
                for(int offset = -1; offset < 2; offset++){
                    int i = (node_index + offset + solution.size()) % solution.size();
                    int old_cost = 0;
                    old_cost += get_cost(solution[(i-1 + solution.size()) % solution.size()], solution[i], dataset, distance_matrix);
                    old_cost += get_cost(solution[i], solution[(i+1) % solution.size()], dataset, distance_matrix);
                    for(int j = 0; j < available_nodes.size(); j++){
                        int cost_improvement = 0;
                        cost_improvement -= get_cost(solution[(i-1 + solution.size()) % solution.size()], available_nodes[j], dataset, distance_matrix);
                        cost_improvement -= get_cost(available_nodes[j], solution[(i+1) % solution.size()], dataset, distance_matrix);
                        cost_improvement += old_cost;
                        if(cost_improvement > 0){
                            std::vector <int> move{cost_improvement, 0, solution[i], solution[(i-1 + solution.size()) % solution.size()], solution[(i+1) % solution.size()], available_nodes[j]};

                            to_insert.push_back(move);
                            
                            //improvement_list.insert(improvement_list.begin() + get_index_to_insert_descending(improvement_list, move), move);
                        }
                    }
                }
                // New moves replacing with the new available node
                for(int i = 0; i < solution.size(); i++){
                    int j = replacement_index;

                    int old_cost = 0;
                    old_cost += get_cost(solution[(i-1 + solution.size()) % solution.size()], solution[i], dataset, distance_matrix);
                    old_cost += get_cost(solution[i], solution[(i+1) % solution.size()], dataset, distance_matrix);

                    int cost_improvement = 0;

                    cost_improvement -= get_cost(solution[(i-1 + solution.size()) % solution.size()], available_nodes[j], dataset, distance_matrix);

                    cost_improvement -= get_cost(available_nodes[j], solution[(i+1) % solution.size()], dataset, distance_matrix);
                    cost_improvement += old_cost;

                    if(cost_improvement > 0){
                        std::vector <int> move{cost_improvement, 0, solution[i], solution[(i-1 + solution.size()) % solution.size()], solution[(i+1) % solution.size()], available_nodes[j]};
                        
                        to_insert.push_back(move);
                    }
                    
                }

                // add new edge moves for this new node
                for(int i = 0; i < solution.size(); i++){
                    int j = node_index;
                    if(j == i || j == (i+1) % solution.size() || (j+1) % solution.size() == i) continue; // catch intersections
                    int old_cost = 0;
                    int cost_improvement = 0;

                    old_cost += get_cost(solution[i], solution[(i+1) % solution.size()], dataset, distance_matrix);
                    old_cost += get_cost(solution[j], solution[(j+1) % solution.size()], dataset, distance_matrix);
                    old_cost += dataset[solution[i]][2];

                    cost_improvement -= get_cost(solution[(j+1) % solution.size()], solution[(i+1) % solution.size()], dataset, distance_matrix);
                    cost_improvement -= get_cost(solution[j], solution[i], dataset, distance_matrix);
                    cost_improvement -= dataset[solution[(j+1) % solution.size()]][2];

                    cost_improvement += old_cost;

                    if(cost_improvement > 0){
                        std::vector <int> move{cost_improvement, 1, solution[i], solution[(i+1) % solution.size()], solution[j], solution[(j+1) % solution.size()]};
                        
                        to_insert.push_back(move);
                    }
                
                }
                for(int i = 0; i < solution.size(); i++){
                    int j = (node_index-1 + solution.size()) % solution.size();
                    if(j == i || j == (i+1) % solution.size() || (j+1) % solution.size() == i) continue; // catch intersections
                    int old_cost = 0;
                    int cost_improvement = 0;

                    old_cost += get_cost(solution[i], solution[(i+1) % solution.size()], dataset, distance_matrix);  // This could be calculated outside
                    old_cost += get_cost(solution[j], solution[(j+1) % solution.size()], dataset, distance_matrix);
                    old_cost += dataset[solution[i]][2];

                    cost_improvement -= get_cost(solution[(j+1) % solution.size()], solution[(i+1) % solution.size()], dataset, distance_matrix);
                    cost_improvement -= get_cost(solution[j], solution[i], dataset, distance_matrix);
                    cost_improvement -= dataset[solution[(j+1) % solution.size()]][2];

                    cost_improvement += old_cost;

                    if(cost_improvement > 0){
                        std::vector <int> move{cost_improvement, 1, solution[i], solution[(i+1) % solution.size()], solution[j], solution[(j+1) % solution.size()]};
                        
                        to_insert.push_back(move);
                    }
                
                }
                ///*
                // Switched edges
                for(int i = 0; i < solution.size(); i++){
                    int j = node_index;
                    if(j == i || j == (i+1) % solution.size() || (j+1) % solution.size() == i) continue; // catch intersections
                    int old_cost = 0;
                    int cost_improvement = 0;

                    old_cost += get_cost(solution[i], solution[(i+1) % solution.size()], dataset, distance_matrix);
                    old_cost += get_cost(solution[(j+1) % solution.size()], solution[j], dataset, distance_matrix); // Switched
                    old_cost += dataset[solution[i]][2];

                    cost_improvement -= get_cost(solution[j], solution[(i+1) % solution.size()], dataset, distance_matrix);// Switched
                    cost_improvement -= get_cost(solution[(j+1) % solution.size()], solution[i], dataset, distance_matrix);// Switched
                    cost_improvement -= dataset[solution[j]][2];// Switched

                    cost_improvement += old_cost;

                    if(cost_improvement > 0){
                        std::vector <int> move{cost_improvement, 1, solution[i], solution[(i+1) % solution.size()], solution[(j+1) % solution.size()], solution[j]};
                        
                        to_insert.push_back(move);
                    }
                
                }
                for(int i = 0; i < solution.size(); i++){
                    int j = (node_index-1 + solution.size()) % solution.size();
                    if(j == i || j == (i+1) % solution.size() || (j+1) % solution.size() == i) continue; // catch intersections
                    int old_cost = 0;
                    int cost_improvement = 0;

                    old_cost += get_cost(solution[i], solution[(i+1) % solution.size()], dataset, distance_matrix);
                    old_cost += get_cost(solution[(j+1) % solution.size()], solution[j], dataset, distance_matrix);
                    old_cost += dataset[solution[i]][2];

                    cost_improvement -= get_cost(solution[j], solution[(i+1) % solution.size()], dataset, distance_matrix);
                    cost_improvement -= get_cost(solution[(j+1) % solution.size()], solution[i], dataset, distance_matrix);
                    cost_improvement -= dataset[solution[j]][2];

                    cost_improvement += old_cost;

                    if(cost_improvement > 0){
                        std::vector <int> move{cost_improvement, 1, solution[i], solution[(i+1) % solution.size()], solution[(j+1) % solution.size()], solution[j]};
                        
                        to_insert.push_back(move);
                    }
                
                }
                //*/
                
                // Remove moves which were not applicable
                improvement_list.erase(improvement_list.begin(), improvement_list.begin() + index+1);
                break;
            }
            else{
                // <improvement, 1 (edge exchange), edge1 start node, edge1 end node, edge2 start node, edge2 end node

                int edge1_start_index = 0;
                for(;edge1_start_index < solution.size() && solution[edge1_start_index] != improvement_list[index][2]; edge1_start_index++);
                if(edge1_start_index == solution.size()){ // No longer in solution, so remove
                    continue;
                }
                int edge2_start_index = 0;
                for(;edge2_start_index < solution.size() && solution[edge2_start_index] != improvement_list[index][4]; edge2_start_index++);
                if(edge2_start_index == solution.size()){ // No longer in solution, so remove
                    continue;
                }
                
                if(((solution[(edge1_start_index+1) % solution.size()] == improvement_list[index][3] &&
                    solution[(edge2_start_index-1 + solution.size()) % solution.size()] == improvement_list[index][5]) || 
                    (solution[(edge1_start_index-1 + solution.size()) % solution.size()] == improvement_list[index][3] &&
                    solution[(edge2_start_index+1) % solution.size()] == improvement_list[index][5]))){ // Twisted, 1 edge is in a different direction than the other, save for later
                    to_insert.push_back(improvement_list[index]);
                    continue;
                }
                if(!((solution[(edge1_start_index+1) % solution.size()] == improvement_list[index][3] &&
                    solution[(edge2_start_index+1) % solution.size()] == improvement_list[index][5]) || 
                    (solution[(edge1_start_index-1 + solution.size()) % solution.size()] == improvement_list[index][3] &&
                    solution[(edge2_start_index-1 + solution.size()) % solution.size()] == improvement_list[index][5]))){ // If not twisted and not placed correctly, then move no longer exists
                    continue;
                }
                
                // Perform operation in reverse if both edges are reversed
                if((solution[(edge1_start_index-1 + solution.size()) % solution.size()] == improvement_list[index][3] &&
                    solution[(edge2_start_index-1 + solution.size()) % solution.size()] == improvement_list[index][5])){ // Both reversed, so just take a step back
                    edge1_start_index = (edge1_start_index-1 + solution.size()) % solution.size();
                    edge2_start_index = (edge2_start_index-1 + solution.size()) % solution.size();
                }

                // FOUND legal edge exchange
                found_improvement = true;
                // update the solution
                if(((edge1_start_index + 1) % solution.size()) > ((edge2_start_index + 1) % solution.size())){
                    int temp = edge1_start_index;
                    edge1_start_index = edge2_start_index;
                    edge2_start_index = temp;
                }

                std::reverse(solution.begin() + ((edge1_start_index + 1) % solution.size()), solution.begin() + ((edge2_start_index + 1) % solution.size()));
                
                //std::cout << "[imp " << improvement_list[index][0] << " edge " << solution[edge1_start_index] << "->" << solution[(edge1_start_index + 1) % solution.size()];
                //std::cout << " into " << solution[edge2_start_index] << "->" << solution[(edge2_start_index + 1) % solution.size()];
                //std::cout << " in " << edge1_start_index << " and " << edge2_start_index << "], ";    

                // Add new moves for all nodes of the exchanged edges
                std::vector <int> nodes_to_update{edge1_start_index, int((edge1_start_index + 1) % solution.size()), edge2_start_index, int((edge2_start_index + 1) % solution.size())};
                for(int i = 0; i < solution.size(); i++){
                    for(int k = 0; k < 2; k++){
                        int j = nodes_to_update[k*2];
                        if(j == i || j == (i+1) % solution.size() || (j+1) % solution.size() == i) continue; // catch intersections
                        int old_cost = 0;
                        int cost_improvement = 0;

                        old_cost += get_cost(solution[i], solution[(i+1) % solution.size()], dataset, distance_matrix);
                        old_cost += get_cost(solution[j], solution[(j+1) % solution.size()], dataset, distance_matrix);
                        old_cost += dataset[solution[i]][2];

                        cost_improvement -= get_cost(solution[(j+1) % solution.size()], solution[(i+1) % solution.size()], dataset, distance_matrix);
                        cost_improvement -= get_cost(solution[j], solution[i], dataset, distance_matrix);
                        cost_improvement -= dataset[solution[(j+1) % solution.size()]][2];

                        cost_improvement += old_cost;

                        if(cost_improvement > 0){
                            std::vector <int> move{cost_improvement, 1, solution[i], solution[(i+1) % solution.size()], solution[j], solution[(j+1) % solution.size()]};
                            to_insert.push_back(move);

                        }
                    }
                }
                for(int i = 0; i < solution.size(); i++){
                    for(int k = 0; k < 2; k++){
                        int j = nodes_to_update[k*2];
                        if(j == i || j == (i+1) % solution.size() || (j+1) % solution.size() == i) continue; // catch intersections
                        int old_cost = 0;
                        int cost_improvement = 0;

                        old_cost += get_cost(solution[i], solution[(i+1) % solution.size()], dataset, distance_matrix);
                        old_cost += get_cost(solution[(j+1) % solution.size()], solution[j], dataset, distance_matrix); // Switch
                        old_cost += dataset[solution[i]][2];

                        cost_improvement -= get_cost(solution[j], solution[(i+1) % solution.size()], dataset, distance_matrix);
                        cost_improvement -= get_cost(solution[(j+1) % solution.size()], solution[i], dataset, distance_matrix);
                        cost_improvement -= dataset[solution[j]][2];

                        cost_improvement += old_cost;

                        if(cost_improvement > 0){
                            std::vector <int> move{cost_improvement, 1, solution[i], solution[(i+1) % solution.size()], solution[(j+1) % solution.size()], solution[j]};
                            to_insert.push_back(move);

                        }
                    }
                }
                
                // No reversed
                /*
                for(int i = 0; i < solution.size(); i++){
                    for(int j = int((edge1_start_index + 1) % solution.size()); j <= int((edge2_start_index + 1) % solution.size()); j++){
                        //int j = nodes_to_update[k];
                        if(j == i || j == (i+1) % solution.size() || (j+1) % solution.size() == i) continue; // catch intersections
                        int old_cost = 0;
                        int cost_improvement = 0;

                        old_cost += get_cost(solution[i], solution[(i+1) % solution.size()], dataset, distance_matrix);
                        old_cost += get_cost(solution[j], solution[(j+1) % solution.size()], dataset, distance_matrix);
                        old_cost += dataset[solution[i]][2];

                        cost_improvement -= get_cost(solution[(j+1) % solution.size()], solution[(i+1) % solution.size()], dataset, distance_matrix);
                        cost_improvement -= get_cost(solution[j], solution[i], dataset, distance_matrix);
                        cost_improvement -= dataset[solution[(j+1) % solution.size()]][2];

                        cost_improvement += old_cost;

                        if(cost_improvement > 0){
                            std::vector <int> move{cost_improvement, 1, solution[i], solution[(i+1) % solution.size()], solution[j], solution[(j+1) % solution.size()]};
                            to_insert.push_back(move);
                        }
                    }

                    int j = edge1_start_index;
                    if(j == i || j == (i+1) % solution.size() || (j+1) % solution.size() == i) continue; // catch intersections
                    int old_cost = 0;
                    int cost_improvement = 0;

                    old_cost += get_cost(solution[i], solution[(i+1) % solution.size()], dataset, distance_matrix);
                    old_cost += get_cost(solution[j], solution[(j+1) % solution.size()], dataset, distance_matrix);
                    old_cost += dataset[solution[i]][2];

                    cost_improvement -= get_cost(solution[(j+1) % solution.size()], solution[(i+1) % solution.size()], dataset, distance_matrix);
                    cost_improvement -= get_cost(solution[j], solution[i], dataset, distance_matrix);
                    cost_improvement -= dataset[solution[(j+1) % solution.size()]][2];

                    cost_improvement += old_cost;

                    if(cost_improvement > 0){
                        std::vector <int> move{cost_improvement, 1, solution[i], solution[(i+1) % solution.size()], solution[j], solution[(j+1) % solution.size()]};
                        to_insert.push_back(move);
                    }
                }
                */
                // And nodes
                for(int k = 0; k < nodes_to_update.size(); k++){
                    int i = nodes_to_update[k];
                    int old_cost = 0;
                    old_cost += get_cost(solution[(i-1 + solution.size()) % solution.size()], solution[i], dataset, distance_matrix);
                    old_cost += get_cost(solution[i], solution[(i+1) % solution.size()], dataset, distance_matrix);
                    for(int j = 0; j < available_nodes.size(); j++){
                        int cost_improvement = 0;
                        cost_improvement -= get_cost(solution[(i-1 + solution.size()) % solution.size()], available_nodes[j], dataset, distance_matrix);
                        cost_improvement -= get_cost(available_nodes[j], solution[(i+1) % solution.size()], dataset, distance_matrix);
                        cost_improvement += old_cost;
                        if(cost_improvement > 0){
                            std::vector <int> move{cost_improvement, 0, solution[i], solution[(i-1 + solution.size()) % solution.size()], solution[(i+1) % solution.size()], available_nodes[j]};
                            to_insert.push_back(move);
                        }
                    }
                }

                // Remove moves which were not applicable
                improvement_list.erase(improvement_list.begin(), improvement_list.begin() + index+1);
                break;
            }
        }

        if(!found_improvement){
            return solution;
        }

        for(int i = 0; i < to_insert.size(); i++){
            improvement_list.insert(improvement_list.begin() + get_index_to_insert_descending(improvement_list, to_insert[i]), to_insert[i]);
        }
        
    }

    return solution;
}


std::vector <int> local_delta_indexes(std::vector <int> solution, std::vector <std::vector <int>> & dataset, std::vector <std::vector <int>> & distance_matrix){
    std::vector <int> available_nodes;
    std::vector <int> unavailable_nodes = solution;
    std::sort(unavailable_nodes.begin(), unavailable_nodes.end());

    int index = 0;
    for(int i = 0; i < dataset.size(); i++){
        if(index < unavailable_nodes.size() && unavailable_nodes[index] == i){
            index++;
        }
        else{
            available_nodes.push_back(i);
        }
    }
    
    std::vector <std::vector <int>> improvement_list; // <improvement, (0 - node replacement / 1 - edge exchange), (cycle node, available node) / (edge1 start, edge1 end, edge2 start, edge2 end)>
    // <improvement, 0 (node replacement), node index, node, node before, node after, replacement index, replacement node>
    // <improvement, 1 (edge exchange), edge1 start index, edge1 start node, edge1 end node, edge2 start index, edge2 start node, edge2 end node

    // Initiate improvement list
    for(int i = 0; i < solution.size(); i++){
        int old_cost = 0;
        old_cost += get_cost(solution[(i-1 + solution.size()) % solution.size()], solution[i], dataset, distance_matrix);
        old_cost += get_cost(solution[i], solution[(i+1) % solution.size()], dataset, distance_matrix);
        for(int j = 0; j < available_nodes.size(); j++){
            int cost_improvement = 0;
            cost_improvement -= get_cost(solution[(i-1 + solution.size()) % solution.size()], available_nodes[j], dataset, distance_matrix);
            cost_improvement -= get_cost(available_nodes[j], solution[(i+1) % solution.size()], dataset, distance_matrix);
            cost_improvement += old_cost;
            if(cost_improvement > 0){
                std::vector <int> move{cost_improvement, 0, i, solution[i], solution[(i-1 + solution.size()) % solution.size()], solution[(i+1) % solution.size()], j, available_nodes[j]};
                improvement_list.push_back(move);
            }
        }
    }

    for(int i = 0; i < solution.size(); i++){
        for(int j = 0; (j+1) % solution.size() < i; j++){
            if(j == i || j == (i+1) % solution.size() || (j+1) % solution.size() == i) continue; // catch intersections
            int old_cost = 0;
            int cost_improvement = 0;

            old_cost += get_cost(solution[i], solution[(i+1) % solution.size()], dataset, distance_matrix);
            old_cost += get_cost(solution[j], solution[(j+1) % solution.size()], dataset, distance_matrix);
            old_cost += dataset[solution[i]][2];

            cost_improvement -= get_cost(solution[(j+1) % solution.size()], solution[(i+1) % solution.size()], dataset, distance_matrix);
            cost_improvement -= get_cost(solution[j], solution[i], dataset, distance_matrix);
            cost_improvement -= dataset[solution[(j+1) % solution.size()]][2];

            cost_improvement += old_cost;

            if(cost_improvement > 0){
                std::vector <int> move{cost_improvement, 1, i, solution[i], solution[(i+1) % solution.size()], j, solution[j], solution[(j+1) % solution.size()]};
                improvement_list.push_back(move);
            }
        }
    }

    // Sort the available moves
    std::sort(improvement_list.begin(), improvement_list.end(), sort_by_first);

    // 100 iterations instead of while(true)
    for(int t = 0; t < 1000; t++){
        // Evaluate new moves
        bool found_improvement = false;
        std::vector <std::vector <int>> to_insert;

        for(int index = 0; index < improvement_list.size(); index++){

            if(improvement_list[index][1] == 0){
                
                /*int node_index = 0;
                for(;node_index < solution.size() && solution[node_index] != improvement_list[index][2]; node_index++);
                if(node_index == solution.size()){ // No longer in solution
                    continue;
                }*/

                // <improvement, 0 (node replacement), node index, node, node before, node after, replacement index, replacement node>

                if(solution[improvement_list[index][2]] != improvement_list[index][3]) continue; // Node no longer there

                if(!((solution[(improvement_list[index][2]-1 + solution.size()) % solution.size()] == improvement_list[index][4] &&
                    solution[(improvement_list[index][2]+1) % solution.size()] == improvement_list[index][5]) || 
                    (solution[(improvement_list[index][2]-1 + solution.size()) % solution.size()] == improvement_list[index][5] &&
                    solution[(improvement_list[index][2]+1) % solution.size()] == improvement_list[index][4]))){ // If no longer 3 nodes in a row (either way)
                    continue;
                }
                
                /*int replacement_index = 0;
                for(;replacement_index < available_nodes.size() && available_nodes[replacement_index] != improvement_list[index][5]; replacement_index++);
                if(replacement_index == available_nodes.size()){ // No longer in available
                    continue;
                }*/
                if(available_nodes[improvement_list[index][6]] != improvement_list[index][7]) continue;

                // FOUND legal node replacement
                
                found_improvement = true;
                // update the solution

                int temp = solution[improvement_list[index][2]];
                solution[improvement_list[index][2]] = available_nodes[improvement_list[index][6]];
                available_nodes[improvement_list[index][6]] = temp;

                // add new moves for the replaced node, node before and after
                for(int offset = -1; offset < 2; offset++){
                    int i = (improvement_list[index][2] + offset + solution.size()) % solution.size();
                    int old_cost = 0;
                    old_cost += get_cost(solution[(i-1 + solution.size()) % solution.size()], solution[i], dataset, distance_matrix);
                    old_cost += get_cost(solution[i], solution[(i+1) % solution.size()], dataset, distance_matrix);
                    for(int j = 0; j < available_nodes.size(); j++){
                        int cost_improvement = 0;
                        cost_improvement -= get_cost(solution[(i-1 + solution.size()) % solution.size()], available_nodes[j], dataset, distance_matrix);
                        cost_improvement -= get_cost(available_nodes[j], solution[(i+1) % solution.size()], dataset, distance_matrix);
                        cost_improvement += old_cost;
                        if(cost_improvement > 0){
                            std::vector <int> move{cost_improvement, 0, i, solution[i], solution[(i-1 + solution.size()) % solution.size()], solution[(i+1) % solution.size()], j, available_nodes[j]};

                            to_insert.push_back(move);
                            
                            //improvement_list.insert(improvement_list.begin() + get_index_to_insert_descending(improvement_list, move), move);
                        }
                    }
                }
                // New moves replacing with the new available node
                for(int i = 0; i < solution.size(); i++){
                    int j = improvement_list[index][6];

                    int old_cost = 0;
                    old_cost += get_cost(solution[(i-1 + solution.size()) % solution.size()], solution[i], dataset, distance_matrix);
                    old_cost += get_cost(solution[i], solution[(i+1) % solution.size()], dataset, distance_matrix);

                    int cost_improvement = 0;

                    cost_improvement -= get_cost(solution[(i-1 + solution.size()) % solution.size()], available_nodes[j], dataset, distance_matrix);

                    cost_improvement -= get_cost(available_nodes[j], solution[(i+1) % solution.size()], dataset, distance_matrix);
                    cost_improvement += old_cost;

                    if(cost_improvement > 0){
                        std::vector <int> move{cost_improvement, 0, i, solution[i], solution[(i-1 + solution.size()) % solution.size()], solution[(i+1) % solution.size()], j, available_nodes[j]};
                        
                        to_insert.push_back(move);
                    }
                    
                }

                // add new edge moves for this new node
                for(int i = 0; i < solution.size(); i++){
                    int j = improvement_list[index][2];
                    if(j == i || j == (i+1) % solution.size() || (j+1) % solution.size() == i) continue; // catch intersections
                    int old_cost = 0;
                    int cost_improvement = 0;

                    old_cost += get_cost(solution[i], solution[(i+1) % solution.size()], dataset, distance_matrix);
                    old_cost += get_cost(solution[j], solution[(j+1) % solution.size()], dataset, distance_matrix);
                    old_cost += dataset[solution[i]][2];

                    cost_improvement -= get_cost(solution[(j+1) % solution.size()], solution[(i+1) % solution.size()], dataset, distance_matrix);
                    cost_improvement -= get_cost(solution[j], solution[i], dataset, distance_matrix);
                    cost_improvement -= dataset[solution[(j+1) % solution.size()]][2];

                    cost_improvement += old_cost;

                    if(cost_improvement > 0){
                        std::vector <int> move{cost_improvement, 1, i, solution[i], solution[(i+1) % solution.size()], j, solution[j], solution[(j+1) % solution.size()]};
                        
                        to_insert.push_back(move);
                    }
                
                }
                for(int i = 0; i < solution.size(); i++){
                    int j = (improvement_list[index][2]-1 + solution.size()) % solution.size();
                    if(j == i || j == (i+1) % solution.size() || (j+1) % solution.size() == i) continue; // catch intersections
                    int old_cost = 0;
                    int cost_improvement = 0;

                    old_cost += get_cost(solution[i], solution[(i+1) % solution.size()], dataset, distance_matrix);
                    old_cost += get_cost(solution[j], solution[(j+1) % solution.size()], dataset, distance_matrix);
                    old_cost += dataset[solution[i]][2];

                    cost_improvement -= get_cost(solution[(j+1) % solution.size()], solution[(i+1) % solution.size()], dataset, distance_matrix);
                    cost_improvement -= get_cost(solution[j], solution[i], dataset, distance_matrix);
                    cost_improvement -= dataset[solution[(j+1) % solution.size()]][2];

                    cost_improvement += old_cost;

                    if(cost_improvement > 0){
                        std::vector <int> move{cost_improvement, 1, i, solution[i], solution[(i+1) % solution.size()], j, solution[j], solution[(j+1) % solution.size()]};
                        
                        to_insert.push_back(move);
                    }
                
                }
                
                // Remove moves which were not applicable
                improvement_list.erase(improvement_list.begin(), improvement_list.begin() + index+1);
                break;
            }
            else{
                // <improvement, 1 (edge exchange), edge1 start index, edge1 start node, edge1 end node, edge2 start index, edge2 start node, edge2 end node
                if(solution[improvement_list[index][2]] != improvement_list[index][3]) continue;
                if(solution[improvement_list[index][5]] != improvement_list[index][6]) continue;

                if(solution[(improvement_list[index][2] + 1) % solution.size()] != improvement_list[index][4]) continue;
                if(solution[(improvement_list[index][5] + 1) % solution.size()] != improvement_list[index][7]) continue;


                /*int edge1_start_index = 0;
                for(;edge1_start_index < solution.size() && solution[edge1_start_index] != improvement_list[index][2]; edge1_start_index++);
                if(edge1_start_index == solution.size()){ // No longer in solution, so remove
                    continue;
                }
                int edge2_start_index = 0;
                for(;edge2_start_index < solution.size() && solution[edge2_start_index] != improvement_list[index][4]; edge2_start_index++);
                if(edge2_start_index == solution.size()){ // No longer in solution, so remove
                    continue;
                }*/
                /*
                if(((solution[(edge1_start_index+1) % solution.size()] == improvement_list[index][3] &&
                    solution[(edge2_start_index-1 + solution.size()) % solution.size()] == improvement_list[index][5]) || 
                    (solution[(edge1_start_index-1 + solution.size()) % solution.size()] == improvement_list[index][3] &&
                    solution[(edge2_start_index+1) % solution.size()] == improvement_list[index][5]))){ // Twisted, 1 edge i sin a different direction than the other, save for later
                    saved_edges.push_back(improvement_list[index]);
                    continue;
                }
                if(!((solution[(edge1_start_index+1) % solution.size()] == improvement_list[index][3] &&
                    solution[(edge2_start_index+1) % solution.size()] == improvement_list[index][5]) || 
                    (solution[(edge1_start_index-1 + solution.size()) % solution.size()] == improvement_list[index][3] &&
                    solution[(edge2_start_index-1 + solution.size()) % solution.size()] == improvement_list[index][5]))){ // If not twisted and not placed correctly, then move no longer exists
                    continue; // !!! THIS IS BUGGED, AS LOOKING THE OTHER WAY MEANS THE REVERSE OPERATION SHOULD BE HANDLED DIFFERENTLY
                }*/

                // FOUND legal edge exchange
                found_improvement = true;
                // update the solution
                if(((improvement_list[index][2] + 1) % solution.size()) > ((improvement_list[index][5] + 1) % solution.size())){
                    int temp = improvement_list[index][2];
                    improvement_list[index][2] = improvement_list[index][5];
                    improvement_list[index][5] = temp;
                }

                std::reverse(solution.begin() + ((improvement_list[index][2] + 1) % solution.size()), solution.begin() + ((improvement_list[index][5] + 1) % solution.size()));
                

                // Add new moves for all nodes which were reversed (including edge1 start, end, edge2 start and end and all inbetween)
                for(int i = (improvement_list[index][2] + 1) % solution.size(); i <= (improvement_list[index][5] + 1) % solution.size(); i++){
                    int old_cost = 0;
                    old_cost += get_cost(solution[(i-1 + solution.size()) % solution.size()], solution[i], dataset, distance_matrix);
                    old_cost += get_cost(solution[i], solution[(i+1) % solution.size()], dataset, distance_matrix);
                    for(int j = 0; j < available_nodes.size(); j++){
                        int cost_improvement = 0;
                        cost_improvement -= get_cost(solution[(i-1 + solution.size()) % solution.size()], available_nodes[j], dataset, distance_matrix);
                        cost_improvement -= get_cost(available_nodes[j], solution[(i+1) % solution.size()], dataset, distance_matrix);
                        cost_improvement += old_cost;
                        if(cost_improvement > 0){
                            std::vector <int> move{cost_improvement, 0, i, solution[i], solution[(i-1 + solution.size()) % solution.size()], solution[(i+1) % solution.size()], j, available_nodes[j]};
                            to_insert.push_back(move);
                        }
                    }
                }
                int i = improvement_list[index][2];
                int old_cost = 0;
                old_cost += get_cost(solution[(i-1 + solution.size()) % solution.size()], solution[i], dataset, distance_matrix);
                old_cost += get_cost(solution[i], solution[(i+1) % solution.size()], dataset, distance_matrix);
                for(int j = 0; j < available_nodes.size(); j++){
                    int cost_improvement = 0;
                    cost_improvement -= get_cost(solution[(i-1 + solution.size()) % solution.size()], available_nodes[j], dataset, distance_matrix);
                    cost_improvement -= get_cost(available_nodes[j], solution[(i+1) % solution.size()], dataset, distance_matrix);
                    cost_improvement += old_cost;
                    if(cost_improvement > 0){
                        std::vector <int> move{cost_improvement, 0, i, solution[i], solution[(i-1 + solution.size()) % solution.size()], solution[(i+1) % solution.size()], j, available_nodes[j]};
                        to_insert.push_back(move);
                    }
                }

                // add new edge moves for this new node
                // New edges exchanges for all nodes in the reversed part into all nodes in the part that did not change
                
                // Experiment with evaluating all edges in this step

                for(int i = (improvement_list[index][2] + 1) % solution.size(); i <= (improvement_list[index][5] + 1) % solution.size(); i++){
                    for(int j = 0; j < solution.size(); j++){
                        //if((improvement_list[index][2] + 1) % solution.size() <= j || j < (improvement_list[index][5] + 1) % solution.size()) continue;
                        if(j == i || j == (i+1) % solution.size() || (j+1) % solution.size() == i) continue; // catch intersections
                        int old_cost = 0;
                        int cost_improvement = 0;

                        old_cost += get_cost(solution[i], solution[(i+1) % solution.size()], dataset, distance_matrix);
                        old_cost += get_cost(solution[j], solution[(j+1) % solution.size()], dataset, distance_matrix);
                        old_cost += dataset[solution[i]][2];

                        cost_improvement -= get_cost(solution[(j+1) % solution.size()], solution[(i+1) % solution.size()], dataset, distance_matrix);
                        cost_improvement -= get_cost(solution[j], solution[i], dataset, distance_matrix);
                        cost_improvement -= dataset[solution[(j+1) % solution.size()]][2];

                        cost_improvement += old_cost;

                        if(cost_improvement > 0){
                            std::vector <int> move{cost_improvement, 1, i, solution[i], solution[(i+1) % solution.size()], j, solution[j], solution[(j+1) % solution.size()]};
                            to_insert.push_back(move);
                        }
                    }
                }
                i = improvement_list[index][2];
                for(int j = 0; j < solution.size(); j++){
                    //if((improvement_list[index][2] + 1) % solution.size() <= j || j < (improvement_list[index][5] + 1) % solution.size()) continue;
                    if(j == i || j == (i+1) % solution.size() || (j+1) % solution.size() == i) continue; // catch intersections
                    int old_cost = 0;
                    int cost_improvement = 0;

                    old_cost += get_cost(solution[i], solution[(i+1) % solution.size()], dataset, distance_matrix);
                    old_cost += get_cost(solution[j], solution[(j+1) % solution.size()], dataset, distance_matrix);
                    old_cost += dataset[solution[i]][2];

                    cost_improvement -= get_cost(solution[(j+1) % solution.size()], solution[(i+1) % solution.size()], dataset, distance_matrix);
                    cost_improvement -= get_cost(solution[j], solution[i], dataset, distance_matrix);
                    cost_improvement -= dataset[solution[(j+1) % solution.size()]][2];

                    cost_improvement += old_cost;

                    if(cost_improvement > 0){
                        std::vector <int> move{cost_improvement, 1, i, solution[i], solution[(i+1) % solution.size()], j, solution[j], solution[(j+1) % solution.size()]};
                        to_insert.push_back(move);
                    }
                }

                // Remove moves which were not applicable
                improvement_list.erase(improvement_list.begin(), improvement_list.begin() + index+1);
                break;
            }
        }

        if(!found_improvement){
            return solution;
        }

        for(int i = 0; i < to_insert.size(); i++){
            improvement_list.insert(improvement_list.begin() + get_index_to_insert_descending(improvement_list, to_insert[i]), to_insert[i]);
        }
        
    }

    return solution;
}


std::vector <int> local_delta(std::vector <int> solution, std::vector <std::vector <int>> & dataset, std::vector <std::vector <int>> & distance_matrix){
    std::vector <int> available_nodes;
    std::vector <int> unavailable_nodes = solution;
    std::sort(unavailable_nodes.begin(), unavailable_nodes.end());

    std::vector <int> node_position;
    std::vector <bool> node_in_cycle;

    //std::cout << "\n\nDELTA \n";
    int index = 0;
    for(int i = 0; i < dataset.size(); i++){
        if(index < unavailable_nodes.size() && unavailable_nodes[index] == i){
            index++;
        }
        else{
            available_nodes.push_back(i);
        }
        node_in_cycle.push_back(true);
        node_position.push_back(-1);
    }

    for(int i = 0; i < solution.size(); i++){
        node_in_cycle[solution[i]] = true;
        node_position[solution[i]] = i;
    }
    for(int i = 0; i < available_nodes.size(); i++){
        node_in_cycle[available_nodes[i]] = false;
        node_position[available_nodes[i]] = i;
    }

    // WRITE ONCE AGAIN, BUT CLEAN
    
    std::priority_queue<std::vector<int>, std::vector<std::vector<int>>, std::less<>> improvement_list;
    //std::vector <std::vector <int>> improvement_list; // <improvement, (0 - node replacement / 1 - edge exchange), (cycle node, available node) / (edge1 start, edge1 end, edge2 start, edge2 end)>
    // <improvement, 0 (node replacement), node, node before, node after, replacement node>
    // <improvement, 1 (edge exchange), edge1 start node, edge1 end node, edge2 start node, edge2 end node>

    // Initiate improvement list
    
    // Edges
    for(int i = 0; i < solution.size(); i++){
        for(int j = 0; j < i - 1; j++){
            if(abs(i - j) < 2 || abs(i - j) == solution.size() - 1) continue; // catch intersections
            int old_cost = 0;
            int cost_improvement = 0;

            old_cost += get_cost(solution[i], solution[(i+1) % solution.size()], dataset, distance_matrix);
            old_cost += get_cost(solution[j], solution[(j+1) % solution.size()], dataset, distance_matrix);
            old_cost += dataset[solution[i]][2];

            cost_improvement -= get_cost(solution[(j+1) % solution.size()], solution[(i+1) % solution.size()], dataset, distance_matrix);
            cost_improvement -= get_cost(solution[j], solution[i], dataset, distance_matrix);
            cost_improvement -= dataset[solution[(j+1) % solution.size()]][2];

            cost_improvement += old_cost;

            if(cost_improvement > 0){
                std::vector <int> move{cost_improvement, 1, solution[i], solution[(i+1) % solution.size()], solution[j], solution[(j+1) % solution.size()]};
                improvement_list.push(move);
            }
        }
    }

    // Reversed Edges
    for(int i = 0; i < solution.size(); i++){
        for(int j = 0; j < i - 1; j++){
            if(abs(i - j) < 2 || abs(i - j) == solution.size() - 1) continue; // catch intersections
            int old_cost = 0;
            int cost_improvement = 0;

            old_cost += get_cost(solution[i], solution[(i+1) % solution.size()], dataset, distance_matrix);
            old_cost += get_cost(solution[(j+1) % solution.size()], solution[j], dataset, distance_matrix);
            old_cost += dataset[solution[i]][2];

            cost_improvement -= get_cost(solution[j], solution[(i+1) % solution.size()], dataset, distance_matrix);
            cost_improvement -= get_cost(solution[(j+1) % solution.size()], solution[i], dataset, distance_matrix);
            cost_improvement -= dataset[solution[j]][2];

            cost_improvement += old_cost;

            if(cost_improvement > 0){
                std::vector <int> move{cost_improvement, 1, solution[i], solution[(i+1) % solution.size()], solution[(j+1) % solution.size()], solution[j]};
                improvement_list.push(move);
            }
        }
    }

    // Nodes
    for(int i = 0; i < solution.size(); i++){
        int old_cost = 0;
        old_cost += get_cost(solution[(i-1 + solution.size()) % solution.size()], solution[i], dataset, distance_matrix);
        old_cost += get_cost(solution[i], solution[(i+1) % solution.size()], dataset, distance_matrix);
        for(int j = 0; j < available_nodes.size(); j++){

            int cost_improvement = 0;
            cost_improvement -= get_cost(solution[(i-1 + solution.size()) % solution.size()], available_nodes[j], dataset, distance_matrix);
            cost_improvement -= get_cost(available_nodes[j], solution[(i+1) % solution.size()], dataset, distance_matrix);
            cost_improvement += old_cost;
            if(cost_improvement > 0){
                std::vector <int> move{cost_improvement, 0, solution[i], solution[(i-1 + solution.size()) % solution.size()], solution[(i+1) % solution.size()], available_nodes[j]};
                improvement_list.push(move);
            }
        }
    }

    std::vector <std::vector <int>> saved_moves;
    int old_node, prev_node, next_node, new_node;
    int old_node_index, prev_node_index, next_node_index, new_node_index;

    int edge1_start, edge1_end, edge2_start, edge2_end;
    int edge1_start_index, edge1_end_index, edge2_start_index, edge2_end_index;
    bool edge1_switched, edge2_switched;

    std::vector <int> current_move;
    while(!improvement_list.empty()){
        current_move = improvement_list.top();
        improvement_list.pop();

        bool applicable = false;
        int type = -1;

        if(current_move[1] == 0){


            type = 0;
            old_node = current_move[2];
            prev_node = current_move[3];
            next_node = current_move[4];
            new_node = current_move[5];


            //auto it = std::find(solution.begin(), solution.end(), old_node); 
            //if(it == solution.end()){ 
            //    continue;
            //}
            
            //old_node_index = it - solution.begin();
            
            if(!node_in_cycle[old_node] || node_in_cycle[new_node]){
                continue;
            }

            old_node_index = node_position[old_node];

            prev_node_index = (old_node_index - 1 + solution.size()) % solution.size();
            next_node_index = (old_node_index + 1) % solution.size();

            if(!(solution[prev_node_index] == next_node && solution[next_node_index] == prev_node) &&
                !(solution[prev_node_index] == prev_node && solution[next_node_index] == next_node)){ // If neither normal nor flipped
                continue;
            }
            
            //it = std::find(available_nodes.begin(), available_nodes.end(), new_node); 
            //if(it == available_nodes.end()){ 
            //    continue;
            //}
            
            //new_node_index = it - available_nodes.begin();
            new_node_index = node_position[new_node];

            applicable = true;

            // If gets to here, then applicable
        }
        else{
            type = 1;
            edge1_start = current_move[2];
            edge1_end = current_move[3];
            edge2_start = current_move[4];
            edge2_end = current_move[5];


            //auto it = std::find(solution.begin(), solution.end(), edge1_start); 
            //if(it == solution.end()){ 
            //    continue;
            //}
            //edge1_start_index = it - solution.begin();

            if(!node_in_cycle[edge1_start] || !node_in_cycle[edge2_start]){
                continue;
            }
            edge1_start_index = node_position[edge1_start];

            edge1_end_index = (edge1_start_index + 1) % solution.size();
            if(solution[edge1_end_index] == edge1_end){
                edge1_switched = false;
            }
            else{
                edge1_end_index = (edge1_start_index - 1 + solution.size()) % solution.size();
                if(solution[edge1_end_index] == edge1_end){
                    edge1_switched = true;
                }
                else{
                    continue;
                }
            }

            //it = std::find(solution.begin(), solution.end(), edge2_start); 
            //if(it == solution.end()){ 
            //    continue;
            //}
            //edge2_start_index = it - solution.begin();
            edge2_start_index = node_position[edge2_start];

            edge2_end_index = (edge2_start_index + 1) % solution.size();
            if(solution[edge2_end_index] == edge2_end){
                edge2_switched = false;
            }
            else{
                edge2_end_index = (edge2_start_index - 1 + solution.size()) % solution.size();
                if(solution[edge2_end_index] == edge2_end){
                    edge2_switched = true;
                }
                else{
                    continue;
                }
            }

            if((!edge1_switched && edge2_switched) || (edge1_switched && !edge2_switched)){ // save for later
                saved_moves.push_back(current_move);
                continue;
            }

            applicable = true;
            // If gets here, then applicable
        }

        if(applicable){
            for(int i = 0; i < saved_moves.size(); i++){
                improvement_list.push(saved_moves[i]);
            }
            saved_moves.clear();

            if(type == 0){

                solution[old_node_index] = new_node;
                available_nodes[new_node_index] = old_node;

                node_in_cycle[old_node] = false;
                node_in_cycle[new_node] = true;
                node_position[old_node] = new_node_index;
                node_position[new_node] = old_node_index;

                // Nodes
                std::vector <int> considered_node_indexes{prev_node_index, old_node_index, next_node_index};
                for(int k = 0; k < considered_node_indexes.size(); k++){
                    int i = considered_node_indexes[k];
                    int old_cost = 0;
                    old_cost += get_cost(solution[(i-1 + solution.size()) % solution.size()], solution[i], dataset, distance_matrix);
                    old_cost += get_cost(solution[i], solution[(i+1) % solution.size()], dataset, distance_matrix);
                    for(int j = 0; j < available_nodes.size(); j++){
                        int cost_improvement = 0;
                        cost_improvement -= get_cost(solution[(i-1 + solution.size()) % solution.size()], available_nodes[j], dataset, distance_matrix);
                        cost_improvement -= get_cost(available_nodes[j], solution[(i+1) % solution.size()], dataset, distance_matrix);
                        cost_improvement += old_cost;
                        if(cost_improvement > 0){
                            std::vector <int> move{cost_improvement, 0, solution[i], solution[(i-1 + solution.size()) % solution.size()], solution[(i+1) % solution.size()], available_nodes[j]};
                            improvement_list.push(move);
                        }
                    }
                }

                // For new available node
                for(int i = 0; i < solution.size(); i++){
                    int old_cost = 0;
                    old_cost += get_cost(solution[(i-1 + solution.size()) % solution.size()], solution[i], dataset, distance_matrix);
                    old_cost += get_cost(solution[i], solution[(i+1) % solution.size()], dataset, distance_matrix);
                    int j = new_node_index;
                    int cost_improvement = 0;
                    cost_improvement -= get_cost(solution[(i-1 + solution.size()) % solution.size()], available_nodes[j], dataset, distance_matrix);
                    cost_improvement -= get_cost(available_nodes[j], solution[(i+1) % solution.size()], dataset, distance_matrix);
                    cost_improvement += old_cost;
                    if(cost_improvement > 0){
                        std::vector <int> move{cost_improvement, 0, solution[i], solution[(i-1 + solution.size()) % solution.size()], solution[(i+1) % solution.size()], available_nodes[j]};
                        improvement_list.push(move);
                    }
                }

                // Edges
                std::vector <int> considered_edge_indexes{prev_node_index, old_node_index};
                for(int i = 0; i < solution.size(); i++){
                    for(int k = 0; k < considered_edge_indexes.size(); k++){
                        int j = considered_edge_indexes[k];
                        if(abs(i - j) < 2 || abs(i - j) == solution.size() - 1) continue; // catch intersections
                        int old_cost = 0;
                        int cost_improvement = 0;

                        old_cost += get_cost(solution[i], solution[(i+1) % solution.size()], dataset, distance_matrix);
                        old_cost += get_cost(solution[j], solution[(j+1) % solution.size()], dataset, distance_matrix);
                        old_cost += dataset[solution[i]][2];

                        cost_improvement -= get_cost(solution[(j+1) % solution.size()], solution[(i+1) % solution.size()], dataset, distance_matrix);
                        cost_improvement -= get_cost(solution[j], solution[i], dataset, distance_matrix);
                        cost_improvement -= dataset[solution[(j+1) % solution.size()]][2];

                        cost_improvement += old_cost;

                        if(cost_improvement > 0){
                            std::vector <int> move{cost_improvement, 1, solution[i], solution[(i+1) % solution.size()], solution[j], solution[(j+1) % solution.size()]};
                            improvement_list.push(move);
                        }
                    }
                }

                // Reversed edges
                for(int i = 0; i < solution.size(); i++){
                    for(int k = 0; k < considered_edge_indexes.size(); k++){
                        int j = considered_edge_indexes[k];
                        if(abs(i - j) < 2 || abs(i - j) == solution.size() - 1) continue; // catch intersections
                        int old_cost = 0;
                        int cost_improvement = 0;

                        old_cost += get_cost(solution[i], solution[(i+1) % solution.size()], dataset, distance_matrix);
                        old_cost += get_cost(solution[(j+1) % solution.size()], solution[j], dataset, distance_matrix);
                        old_cost += dataset[solution[i]][2];

                        cost_improvement -= get_cost(solution[j], solution[(i+1) % solution.size()], dataset, distance_matrix);
                        cost_improvement -= get_cost(solution[(j+1) % solution.size()], solution[i], dataset, distance_matrix);
                        cost_improvement -= dataset[solution[j]][2];

                        cost_improvement += old_cost;

                        if(cost_improvement > 0){
                            std::vector <int> move{cost_improvement, 1, solution[i], solution[(i+1) % solution.size()], solution[(j+1) % solution.size()], solution[j]};
                            improvement_list.push(move);
                        }
                    }
                }
            }
            else{
                if(!edge1_switched && !edge2_switched){ // normal placement
                    int rev1 = std::min(edge1_end_index, edge2_end_index), rev2 = std::max(edge1_end_index, edge2_end_index);
                    for(int k = rev1; k < rev2; k++){
                        node_position[solution[k]] = rev1 + (rev2 - 1 - node_position[solution[k]]);
                    }
                    std::reverse(solution.begin() + rev1, solution.begin() + rev2);

                }
                else if(edge1_switched && edge2_switched){ // flipped placement
                    // Behave as if both edges started 1 earlier
                    int temp;
                    temp = edge1_start_index;
                    edge1_start_index = edge1_end_index;
                    edge1_end_index = temp;
                    temp = edge2_start_index;
                    edge2_start_index = edge2_end_index;
                    edge2_end_index = temp;
                    int rev1 = std::min(edge1_end_index, edge2_end_index), rev2 = std::max(edge1_end_index, edge2_end_index);
                    for(int k = rev1; k < rev2; k++){
                        node_position[solution[k]] = rev1 + (rev2 - 1 - node_position[solution[k]]);
                    }
                    std::reverse(solution.begin() + rev1, solution.begin() + rev2);
                }
                else {
                    continue;
                }

                // Nodes
                std::vector <int> considered_node_indexes{edge1_start_index, edge2_start_index, edge1_end_index, edge2_end_index};
                for(int k = 0; k < considered_node_indexes.size(); k++){
                    int i = considered_node_indexes[k];
                    int old_cost = 0;
                    old_cost += get_cost(solution[(i-1 + solution.size()) % solution.size()], solution[i], dataset, distance_matrix);
                    old_cost += get_cost(solution[i], solution[(i+1) % solution.size()], dataset, distance_matrix);
                    for(int j = 0; j < available_nodes.size(); j++){
                        int cost_improvement = 0;
                        cost_improvement -= get_cost(solution[(i-1 + solution.size()) % solution.size()], available_nodes[j], dataset, distance_matrix);
                        cost_improvement -= get_cost(available_nodes[j], solution[(i+1) % solution.size()], dataset, distance_matrix);
                        cost_improvement += old_cost;
                        if(cost_improvement > 0){
                            std::vector <int> move{cost_improvement, 0, solution[i], solution[(i-1 + solution.size()) % solution.size()], solution[(i+1) % solution.size()], available_nodes[j]};
                            improvement_list.push(move);
                        }
                    }
                }

                // Edges
                std::vector <int> considered_edge_indexes{edge1_start_index, edge2_start_index};
                for(int i = 0; i < solution.size(); i++){
                    for(int k = 0; k < considered_edge_indexes.size(); k++){
                        int j = considered_edge_indexes[k];
                        if(abs(i - j) < 2 || abs(i - j) == solution.size() - 1) continue; // catch intersections
                        int old_cost = 0;
                        int cost_improvement = 0;

                        old_cost += get_cost(solution[i], solution[(i+1) % solution.size()], dataset, distance_matrix);
                        old_cost += get_cost(solution[j], solution[(j+1) % solution.size()], dataset, distance_matrix);
                        old_cost += dataset[solution[i]][2];

                        cost_improvement -= get_cost(solution[(j+1) % solution.size()], solution[(i+1) % solution.size()], dataset, distance_matrix);
                        cost_improvement -= get_cost(solution[j], solution[i], dataset, distance_matrix);
                        cost_improvement -= dataset[solution[(j+1) % solution.size()]][2];

                        cost_improvement += old_cost;

                        if(cost_improvement > 0){
                            std::vector <int> move{cost_improvement, 1, solution[i], solution[(i+1) % solution.size()], solution[j], solution[(j+1) % solution.size()]};
                            improvement_list.push(move);
                        }
                    }
                }

                // Reversed edges
                for(int i = 0; i < solution.size(); i++){
                    for(int k = 0; k < considered_edge_indexes.size(); k++){
                        int j = considered_edge_indexes[k];
                        if(abs(i - j) < 2 || abs(i - j) == solution.size() - 1) continue; // catch intersections
                        int old_cost = 0;
                        int cost_improvement = 0;

                        old_cost += get_cost(solution[i], solution[(i+1) % solution.size()], dataset, distance_matrix);
                        old_cost += get_cost(solution[(j+1) % solution.size()], solution[j], dataset, distance_matrix);
                        old_cost += dataset[solution[i]][2];

                        cost_improvement -= get_cost(solution[j], solution[(i+1) % solution.size()], dataset, distance_matrix);
                        cost_improvement -= get_cost(solution[(j+1) % solution.size()], solution[i], dataset, distance_matrix);
                        cost_improvement -= dataset[solution[j]][2];

                        cost_improvement += old_cost;

                        if(cost_improvement > 0){
                            std::vector <int> move{cost_improvement, 1, solution[i], solution[(i+1) % solution.size()], solution[(j+1) % solution.size()], solution[j]};
                            improvement_list.push(move);
                        }
                    }
                }
            }
        }
    }

    return solution;
}

std::vector <int> multiple_start_LS(std::vector <std::vector <int>> & dataset, std::vector <std::vector <int>> & distance_matrix){
    int iterations = 200; // 200

    std::vector <int> saved_solution = local_search(true, true, random_solution(dataset, distance_matrix), dataset, distance_matrix);
    int saved_cost = get_path_cost(saved_solution, dataset, distance_matrix);

    for(int i = 1; i < iterations; i++){
        std::vector <int> solution = random_solution(dataset, distance_matrix);
        solution = local_search(true, true, solution, dataset, distance_matrix);
        int new_cost = get_path_cost(solution, dataset, distance_matrix);
        if(new_cost < saved_cost){
            saved_cost = new_cost;
            saved_solution = solution;
        }
                
    }

    return saved_solution;
}

std::vector <int> iterated_LS(int stopping_time, int max_changes, int & number_of_runs, std::vector <std::vector <int>> & dataset, std::vector <std::vector <int>> & distance_matrix){
    std::vector <int> best_solution = random_solution(dataset, distance_matrix);
    int best_cost = get_path_cost(best_solution, dataset, distance_matrix);

    std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
    std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
    int elapsed_time = std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count();
    while(elapsed_time < stopping_time){
        number_of_runs++;
        
        // Perturbate solution 
        std::vector <int> solution = best_solution;


        // Create available nodes
        std::vector <int> available_nodes;
        std::vector <int> unavailable_nodes = solution;
        std::sort(unavailable_nodes.begin(), unavailable_nodes.end());
        int index = 0;
        for(int i = 0; i < dataset.size(); i++){
            if(index < unavailable_nodes.size() && unavailable_nodes[index] == i){
                index++;
            }
            else{
                available_nodes.push_back(i);
            }
        }


        // Do changes
        float ratio = 1 - (elapsed_time / stopping_time);
        int change_count = std::round(ratio * max_changes); // Ceil?
        //int change_count = max_changes;
        for(int i = 0; i < change_count; i++){
            if(i % 2){
                // Node exchange
                int old_node = rand() % solution.size();
                int new_node = rand() % available_nodes.size();

                int temp = solution[old_node];
                solution[old_node] = available_nodes[new_node];
                available_nodes[new_node] = temp;
            }
            else{
                // Edge exchange
                int node1 = rand() % solution.size();
                int node2 = rand() % solution.size();
                while(node2 == node1){
                    node2 = rand() % solution.size();
                }

                if(((node1 + 1) % solution.size()) > ((node2 + 1) % solution.size())){
                    int temp = node1;
                    node1 = node2;
                    node2 = temp;
                }
                // Edge exchange through flipping a subpath
                std::reverse(solution.begin() + ((node1 + 1) % solution.size()), solution.begin() + ((node2 + 1) % solution.size()));
                
            }
        }


        // Perform local search
        solution = local_search(true, true, solution, dataset, distance_matrix);

        // Compare
        int cost = get_path_cost(solution, dataset, distance_matrix);
        if(cost < best_cost){
            best_cost = cost;
            best_solution = solution;
        }

        // Update time
        end = std::chrono::steady_clock::now();
        elapsed_time = std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count();
    }

    return best_solution;
}

void compare_MSLS_ILS(std::vector <std::vector <int>> & dataset, std::vector <std::vector <int>> & distance_matrix, std::string filename = "", std::string dataset_name = "example"){
    int iterations = 20; // 20

    std::vector <std::vector <int>> solutions_MSLS;
    std::vector <int> costs_MSLS;

    std::chrono::steady_clock::time_point MSLS_begin = std::chrono::steady_clock::now();

    for(int i = 0; i < iterations; i++){
        //std::cout << i << "\n";
        solutions_MSLS.push_back(multiple_start_LS(dataset, distance_matrix)); 
    }

    std::chrono::steady_clock::time_point MSLS_end = std::chrono::steady_clock::now();

    int time_MSLS = std::chrono::duration_cast<std::chrono::milliseconds>(MSLS_end - MSLS_begin).count();

    for(int i = 0; i < solutions_MSLS.size(); i++){
        costs_MSLS.push_back(get_path_cost(solutions_MSLS[i], dataset, distance_matrix));
    }

    auto min_MSLS = *(std::min_element(costs_MSLS.begin(), costs_MSLS.end()));
    auto max_MSLS = *(std::max_element(costs_MSLS.begin(), costs_MSLS.end()));
    auto const count = static_cast<float>(costs_MSLS.size());
    auto avg_MSLS =  std::reduce(costs_MSLS.begin(), costs_MSLS.end()) / count;

    int best_path_index_MSLS = std::find(costs_MSLS.begin(), costs_MSLS.end(), min_MSLS) - costs_MSLS.begin();
    
    if(filename != ""){
        std::ofstream ofs;
        ofs.open(filename, std::ios_base::app);
        ofs << "\nDATASET " << dataset_name << "\n";
        ofs << "\n" << "MSLS" << "\n";
        ofs << "Best, Average, Worst scores:\n" << min_MSLS << " " << avg_MSLS << " " << max_MSLS << "\n";
        for(int k = 0; k < solutions_MSLS[best_path_index_MSLS].size() - 1; k++){
            ofs << solutions_MSLS[best_path_index_MSLS][k] << ", ";
        }
        ofs << solutions_MSLS[best_path_index_MSLS][solutions_MSLS[best_path_index_MSLS].size() - 1] << "\n";
        ofs << "Time: " << time_MSLS << "ms\n";
    }

    // Experimental ILS

    int stopping_time = time_MSLS / iterations;
    for(int j = 1; j < 2; j++){ // 15 Flat is best
        int max_changes = 10; //j * 5;
        std::vector <std::vector <int>> solutions_ILS;
        std::vector <int> costs_ILS;

        std::chrono::steady_clock::time_point ILS_begin = std::chrono::steady_clock::now();

        int number_of_runs = 0;
        for(int i = 0; i < iterations; i++){
            //std::cout << j << " " << i << "\n";
            solutions_ILS.push_back(iterated_LS(stopping_time, max_changes, number_of_runs, dataset, distance_matrix)); 
        }

        std::chrono::steady_clock::time_point ILS_end = std::chrono::steady_clock::now();
        int time_ILS = std::chrono::duration_cast<std::chrono::milliseconds>(ILS_end - ILS_begin).count();

        for(int i = 0; i < solutions_ILS.size(); i++){
            costs_ILS.push_back(get_path_cost(solutions_ILS[i], dataset, distance_matrix));
        }

        auto min_ILS = *(std::min_element(costs_ILS.begin(), costs_ILS.end()));
        auto max_ILS = *(std::max_element(costs_ILS.begin(), costs_ILS.end()));
        auto avg_ILS =  std::reduce(costs_ILS.begin(), costs_ILS.end()) / count;

        int best_path_index_ILS = std::find(costs_ILS.begin(), costs_ILS.end(), min_ILS) - costs_ILS.begin();

        if(filename != ""){
            std::ofstream ofs;
            ofs.open(filename, std::ios_base::app);

            ofs << "\nILS " << max_changes << "\n";
            ofs << "Best, Average, Worst scores:\n" << min_ILS << " " << avg_ILS << " " << max_ILS << "\n";
            for(int k = 0; k < solutions_ILS[best_path_index_ILS].size() - 1; k++){
                ofs << solutions_ILS[best_path_index_ILS][k] << ", ";
            }
            ofs << solutions_ILS[best_path_index_ILS][solutions_ILS[best_path_index_ILS].size() - 1] << "\n";
            ofs << "Time: " << time_ILS << "ms\n";
            ofs << "Number of runs per search: " << number_of_runs / iterations << "\n";
            
            ofs.close();
        }
    }
}

std::vector <int> greedy_weighted_partial(std::vector <int> starting_solution, std::vector <std::vector <int>> & dataset, std::vector <std::vector <int>> & distance_matrix){
    std::vector <int> available_nodes;
    std::vector <int> unavailable_nodes = starting_solution;
    std::sort(unavailable_nodes.begin(), unavailable_nodes.end());

    int index = 0;
    for(int i = 0; i < dataset.size(); i++){
        if(index < unavailable_nodes.size() && unavailable_nodes[index] == i){
            index++;
        }
        else{
            available_nodes.push_back(i);
        }
    }

    int path_target_length = std::floor(dataset.size() / 2);

    if(starting_solution.size() < 3){
        std::cout << "ERROR: Invalid starting path length\n";
        return starting_solution;
    }

    std::vector <int> solution = starting_solution;

    for(int i = solution.size(); i < path_target_length; i++){
        int highest_regret = -1000000;
        int chosen_available_index = -1;
        int chosen_solution_index = -1;

        for(int j = 0; j < available_nodes.size(); j++){
            int best_location_for_this_node = -1;
            std::vector <int> top_2_scores;
            top_2_scores.push_back(1000000);
            top_2_scores.push_back(1000000);

            for(int k = 0; k < solution.size(); k++){
                int extra_distance = get_cost(solution[k], available_nodes[j], dataset, distance_matrix) + \
                                    get_cost(available_nodes[j], solution[(k+1) % solution.size()], dataset, distance_matrix) - \
                                    get_cost(solution[k], solution[(k+1) % solution.size()], dataset, distance_matrix);

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
                chosen_available_index = j;
                chosen_solution_index = best_location_for_this_node;
            }
        }

        solution.insert(solution.begin() + chosen_solution_index, available_nodes[chosen_available_index]);
        available_nodes.erase(available_nodes.begin() + chosen_available_index);
    }

    return solution;
}

std::vector <int> large_scale_neighborhood_search(int stopping_time, bool random_starting, bool perform_local_search, int & number_of_runs, std::vector <std::vector <int>> & dataset, std::vector <std::vector <int>> & distance_matrix){
    std::vector <int> best_solution;
    if(random_starting){
        best_solution = random_solution(dataset, distance_matrix);
        best_solution = local_search(true, true, best_solution, dataset, distance_matrix); // Randomis always local searched
    }
    else{
        std::vector <int> solution = random_solution(dataset, distance_matrix);
        best_solution.push_back(solution[0]);
        best_solution.push_back(solution[1]);
        best_solution.push_back(solution[2]);
        best_solution = greedy_weighted_partial(best_solution, dataset, distance_matrix);
    }

    // Perform local search (optional)
    if(perform_local_search){
        best_solution = local_search(true, true, best_solution, dataset, distance_matrix);
    }

    int best_cost = get_path_cost(best_solution, dataset, distance_matrix);


    int how_many_to_remove = std::floor(dataset.size() * 0.3);
    // Optional local search

    std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
    std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
    int elapsed_time = std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count();
    while(elapsed_time < stopping_time){
        number_of_runs++;
        
        // Destroy solution 
        std::vector <int> solution = best_solution;

        // Remove nodes with probability <- node cost
        // For each node create pair start-end with (end = start + node cost) and (start = end of previos node + 1)
        // Generate number from 0 to last node end, and remove the node in which this random number resides
        std::vector <int> indexes_to_remove;

        std::vector <int> probability_ranges;
        for(int i = 0; i < solution.size(); i++){
            std::vector <int> probability_end;
            if(probability_ranges.empty()){
                probability_ranges.push_back(0 + dataset[solution[i]][2]);
            }
            else{
                probability_ranges.push_back(probability_ranges[probability_ranges.size()-1] + dataset[solution[i]][2]); // Take node cost as probability
            }
        }

        
        std::random_device dev;
        std::mt19937 rng(dev());
        std::uniform_int_distribution<std::mt19937::result_type> dist(0,probability_ranges[probability_ranges.size()-1]); // distribution in range [1, 6]

        while(indexes_to_remove.size() < how_many_to_remove){
            int random_value = dist(rng);

            // Find which node index that is
            int random_index = std::upper_bound(probability_ranges.begin(), probability_ranges.end(), random_value) - probability_ranges.begin();

            if(indexes_to_remove.empty()){
                indexes_to_remove.push_back(random_index);
                continue;
            }

            auto it = std::find(indexes_to_remove.begin(), indexes_to_remove.end(), random_index);
            if(it != indexes_to_remove.end()){ 
                // If it is already on the list to remove, reroll
                continue;
            }
            // Else add to the list to remove
            indexes_to_remove.push_back(random_index);
        }

        // Now remove these random nodes
        while(!indexes_to_remove.empty()){
            int index = indexes_to_remove[indexes_to_remove.size() - 1];
            indexes_to_remove.pop_back();
            solution.erase(solution.begin() + index); // From last to first so correct indexes should be removed
        }

        // Repair solution
        // Greedy heuristic
        // Modified greedy wieghted which takes partial solution instead of only a starting node
        solution = greedy_weighted_partial(solution, dataset, distance_matrix);

        // Perform local search (optional)
        
        if(perform_local_search){
            solution = local_search(true, true, solution, dataset, distance_matrix);
        }

        // Compare
        int cost = get_path_cost(solution, dataset, distance_matrix);
        if(cost < best_cost){
            best_cost = cost;
            best_solution = solution;
        }

        // Update time
        end = std::chrono::steady_clock::now();
        elapsed_time = std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count();
    }

    return best_solution;
}

void compare_MSLS_LSN(std::vector <std::vector <int>> & dataset, std::vector <std::vector <int>> & distance_matrix, std::string filename = "", std::string dataset_name = "example"){
    int iterations = 20; // 20

    std::vector <std::vector <int>> solutions_MSLS;
    std::vector <int> costs_MSLS;

    std::chrono::steady_clock::time_point MSLS_begin = std::chrono::steady_clock::now();

    for(int i = 0; i < iterations; i++){
        //std::cout << i << "\n";
        solutions_MSLS.push_back(multiple_start_LS(dataset, distance_matrix)); 
    }

    std::chrono::steady_clock::time_point MSLS_end = std::chrono::steady_clock::now();

    int time_MSLS = std::chrono::duration_cast<std::chrono::milliseconds>(MSLS_end - MSLS_begin).count();

    for(int i = 0; i < solutions_MSLS.size(); i++){
        costs_MSLS.push_back(get_path_cost(solutions_MSLS[i], dataset, distance_matrix));
    }

    auto min_MSLS = *(std::min_element(costs_MSLS.begin(), costs_MSLS.end()));
    auto max_MSLS = *(std::max_element(costs_MSLS.begin(), costs_MSLS.end()));
    auto const count = static_cast<float>(costs_MSLS.size());
    auto avg_MSLS =  std::reduce(costs_MSLS.begin(), costs_MSLS.end()) / count;

    int best_path_index_MSLS = std::find(costs_MSLS.begin(), costs_MSLS.end(), min_MSLS) - costs_MSLS.begin();
    
    if(filename != ""){
        std::ofstream ofs;
        ofs.open(filename, std::ios_base::app);
        ofs << "\nDATASET " << dataset_name << "\n";
        ofs << "\n" << "MSLS" << "\n";
        ofs << "Best, Average, Worst scores:\n" << min_MSLS << " " << avg_MSLS << " " << max_MSLS << "\n";
        for(int k = 0; k < solutions_MSLS[best_path_index_MSLS].size() - 1; k++){
            ofs << solutions_MSLS[best_path_index_MSLS][k] << ", ";
        }
        ofs << solutions_MSLS[best_path_index_MSLS][solutions_MSLS[best_path_index_MSLS].size() - 1] << "\n";
        ofs << "Time: " << time_MSLS << "ms\n";
    }

    // Experimental LSN with Random starting and no LS

    int stopping_time = time_MSLS / iterations;
    std::vector <std::vector <int>> solutions_LSN;
    std::vector <int> costs_LSN;

    std::chrono::steady_clock::time_point LSN_begin = std::chrono::steady_clock::now();

    int number_of_runs = 0;
    for(int i = 0; i < iterations; i++){
        solutions_LSN.push_back(large_scale_neighborhood_search(stopping_time, true, false, number_of_runs, dataset, distance_matrix)); 
    }

    std::chrono::steady_clock::time_point LSN_end = std::chrono::steady_clock::now();
    int time_LSN = std::chrono::duration_cast<std::chrono::milliseconds>(LSN_end - LSN_begin).count();

    for(int i = 0; i < solutions_LSN.size(); i++){
        costs_LSN.push_back(get_path_cost(solutions_LSN[i], dataset, distance_matrix));
    }

    auto min_LSN = *(std::min_element(costs_LSN.begin(), costs_LSN.end()));
    auto max_LSN = *(std::max_element(costs_LSN.begin(), costs_LSN.end()));
    auto avg_LSN =  std::reduce(costs_LSN.begin(), costs_LSN.end()) / count;

    int best_path_index_LSN = std::find(costs_LSN.begin(), costs_LSN.end(), min_LSN) - costs_LSN.begin();

    if(filename != ""){
        std::ofstream ofs;
        ofs.open(filename, std::ios_base::app);

        ofs << "\nLSN Random no LS" << "\n";
        ofs << "Best, Average, Worst scores:\n" << min_LSN << " " << avg_LSN << " " << max_LSN << "\n";
        for(int k = 0; k < solutions_LSN[best_path_index_LSN].size() - 1; k++){
            ofs << solutions_LSN[best_path_index_LSN][k] << ", ";
        }
        ofs << solutions_LSN[best_path_index_LSN][solutions_LSN[best_path_index_LSN].size() - 1] << "\n";
        ofs << "Time: " << time_LSN << "ms\n";
        ofs << "Number of runs per search: " << number_of_runs / iterations << "\n";
        
        ofs.close();
    }
    
    // Experimental LSN with Random starting and LS

    solutions_LSN.clear();
    costs_LSN.clear();

    LSN_begin = std::chrono::steady_clock::now();

    number_of_runs = 0;
    for(int i = 0; i < iterations; i++){
        solutions_LSN.push_back(large_scale_neighborhood_search(stopping_time, true, true, number_of_runs, dataset, distance_matrix)); 
    }

    LSN_end = std::chrono::steady_clock::now();
    time_LSN = std::chrono::duration_cast<std::chrono::milliseconds>(LSN_end - LSN_begin).count();

    for(int i = 0; i < solutions_LSN.size(); i++){
        costs_LSN.push_back(get_path_cost(solutions_LSN[i], dataset, distance_matrix));
    }

    min_LSN = *(std::min_element(costs_LSN.begin(), costs_LSN.end()));
    max_LSN = *(std::max_element(costs_LSN.begin(), costs_LSN.end()));
    avg_LSN =  std::reduce(costs_LSN.begin(), costs_LSN.end()) / count;

    best_path_index_LSN = std::find(costs_LSN.begin(), costs_LSN.end(), min_LSN) - costs_LSN.begin();

    if(filename != ""){
        std::ofstream ofs;
        ofs.open(filename, std::ios_base::app);

        ofs << "\nLSN Random LS" << "\n";
        ofs << "Best, Average, Worst scores:\n" << min_LSN << " " << avg_LSN << " " << max_LSN << "\n";
        for(int k = 0; k < solutions_LSN[best_path_index_LSN].size() - 1; k++){
            ofs << solutions_LSN[best_path_index_LSN][k] << ", ";
        }
        ofs << solutions_LSN[best_path_index_LSN][solutions_LSN[best_path_index_LSN].size() - 1] << "\n";
        ofs << "Time: " << time_LSN << "ms\n";
        ofs << "Number of runs per search: " << number_of_runs / iterations << "\n";
        
        ofs.close();
    }

    // Experimental LSN with Greedy starting, without Local Search

    solutions_LSN.clear();
    costs_LSN.clear();

    LSN_begin = std::chrono::steady_clock::now();

    number_of_runs = 0;
    for(int i = 0; i < iterations; i++){
        solutions_LSN.push_back(large_scale_neighborhood_search(stopping_time, false, false, number_of_runs, dataset, distance_matrix)); 
    }

    LSN_end = std::chrono::steady_clock::now();
    time_LSN = std::chrono::duration_cast<std::chrono::milliseconds>(LSN_end - LSN_begin).count();

    for(int i = 0; i < solutions_LSN.size(); i++){
        costs_LSN.push_back(get_path_cost(solutions_LSN[i], dataset, distance_matrix));
    }

    min_LSN = *(std::min_element(costs_LSN.begin(), costs_LSN.end()));
    max_LSN = *(std::max_element(costs_LSN.begin(), costs_LSN.end()));
    avg_LSN =  std::reduce(costs_LSN.begin(), costs_LSN.end()) / count;

    best_path_index_LSN = std::find(costs_LSN.begin(), costs_LSN.end(), min_LSN) - costs_LSN.begin();

    if(filename != ""){
        std::ofstream ofs;
        ofs.open(filename, std::ios_base::app);

        ofs << "\nLSN Greedy no LS" << "\n";
        ofs << "Best, Average, Worst scores:\n" << min_LSN << " " << avg_LSN << " " << max_LSN << "\n";
        for(int k = 0; k < solutions_LSN[best_path_index_LSN].size() - 1; k++){
            ofs << solutions_LSN[best_path_index_LSN][k] << ", ";
        }
        ofs << solutions_LSN[best_path_index_LSN][solutions_LSN[best_path_index_LSN].size() - 1] << "\n";
        ofs << "Time: " << time_LSN << "ms\n";
        ofs << "Number of runs per search: " << number_of_runs / iterations << "\n";
        
        ofs.close();
    }
    

    // Experimental LSN with Greedy starting and Local Search

    solutions_LSN.clear();
    costs_LSN.clear();

    LSN_begin = std::chrono::steady_clock::now();

    number_of_runs = 0;
    for(int i = 0; i < iterations; i++){
        solutions_LSN.push_back(large_scale_neighborhood_search(stopping_time, false, true, number_of_runs, dataset, distance_matrix)); 
    }

    LSN_end = std::chrono::steady_clock::now();
    time_LSN = std::chrono::duration_cast<std::chrono::milliseconds>(LSN_end - LSN_begin).count();

    for(int i = 0; i < solutions_LSN.size(); i++){
        costs_LSN.push_back(get_path_cost(solutions_LSN[i], dataset, distance_matrix));
    }

    min_LSN = *(std::min_element(costs_LSN.begin(), costs_LSN.end()));
    max_LSN = *(std::max_element(costs_LSN.begin(), costs_LSN.end()));
    avg_LSN =  std::reduce(costs_LSN.begin(), costs_LSN.end()) / count;

    best_path_index_LSN = std::find(costs_LSN.begin(), costs_LSN.end(), min_LSN) - costs_LSN.begin();

    if(filename != ""){
        std::ofstream ofs;
        ofs.open(filename, std::ios_base::app);

        ofs << "\nLSN Greedy LS" << "\n";
        ofs << "Best, Average, Worst scores:\n" << min_LSN << " " << avg_LSN << " " << max_LSN << "\n";
        for(int k = 0; k < solutions_LSN[best_path_index_LSN].size() - 1; k++){
            ofs << solutions_LSN[best_path_index_LSN][k] << ", ";
        }
        ofs << solutions_LSN[best_path_index_LSN][solutions_LSN[best_path_index_LSN].size() - 1] << "\n";
        ofs << "Time: " << time_LSN << "ms\n";
        ofs << "Number of runs per search: " << number_of_runs / iterations << "\n";
        
        ofs.close();
    }
}


void calculate_best_paths(std::vector <std::vector <int>> & dataset, std::vector <std::vector <int>> & distance_matrix, std::string filename = "", std::string dataset_name = "example"){
    int iterations = 20;
    std::vector <std::vector <int>> best_paths;
    std::vector <int> best_scores;
    std::vector <int> worst_scores;
    std::vector <float> average_scores;
    std::vector <int> times;
    
    std::vector <std::string> algorithm_names;


    
    // Create candidate nodes for each node
    // <index, distance>
    std::vector <std::vector <std::tuple <int, int>>> closest_nodes;
    for(int i = 0; i < dataset.size(); i++){
        std::vector <std::tuple <int, int>> neighbours_and_distances;
        for(int j = 0; j < dataset.size(); j++){  
            if(i == j){
                continue;
            }
            std::tuple <int, int> distance_and_index = std::make_tuple(get_cost(i, j, dataset, distance_matrix), j);
            neighbours_and_distances.push_back(distance_and_index);

        }

        
        std::sort(neighbours_and_distances.begin(), neighbours_and_distances.end()); 

        std::vector <std::tuple <int, int>> ten_best;

        for(int j = 0; j < 10; j++){
            ten_best.push_back(std::make_tuple(std::get<1>(neighbours_and_distances[j]), std::get<0>(neighbours_and_distances[j])));
        }

        closest_nodes.push_back(ten_best);
    }

    std::vector <std::vector <int>> base_solutions;
    for(int i = 0; i < iterations; i++){
        base_solutions.push_back(random_solution(dataset, distance_matrix));
    }


    /*
    algorithm_names.push_back("Nearest");
    algorithm_names.push_back("Greedy Cycle");
    algorithm_names.push_back("Greedy Regret");
    algorithm_names.push_back("Greedy Weighted");
    /// Then Greedy Weighted starting
    algorithm_names.push_back("Local Search101");
    algorithm_names.push_back("Local Search111");
    algorithm_names.push_back("Local Search001");
    algorithm_names.push_back("Local Search011");
    /// Then Random starting
    algorithm_names.push_back("Local Search100");
    //algorithm_names.push_back("Local Search110");
    //algorithm_names.push_back("Local Search000");
    //algorithm_names.push_back("Local Search010");
    */

    algorithm_names.push_back("Random Baseline");
    algorithm_names.push_back("Local Search110");
    algorithm_names.push_back("Local Delta");

    //algorithm_names.push_back("Cand Comp");
    //algorithm_names.push_back("Local Candidate");


    /*
    best_paths.push_back(nearest_solution(0, dataset, distance_matrix));
    best_paths.push_back(greedy_cycle_solution(0, dataset, distance_matrix));
    best_paths.push_back(greedy_2_regret_solution(0, dataset, distance_matrix));
    best_paths.push_back(greedy_weighted_solution(0, dataset, distance_matrix));

    best_paths.push_back(local_search(true, false, 0, dataset, distance_matrix));
    best_paths.push_back(local_search(true, true, 0, dataset, distance_matrix));
    best_paths.push_back(local_search(false, false, 0, dataset, distance_matrix));
    best_paths.push_back(local_search(false, true, 0, dataset, distance_matrix));

    best_paths.push_back(local_search(true, false, -1, dataset, distance_matrix));
    best_paths.push_back(local_search(true, true, -1, dataset, distance_matrix));
    //best_paths.push_back(local_search(false, false, -1, dataset, distance_matrix));
    //best_paths.push_back(local_search(false, true, -1, dataset, distance_matrix));
    */

    best_paths.push_back(base_solutions[0]);
    best_paths.push_back(local_search(true, true, base_solutions[0], dataset, distance_matrix));
    best_paths.push_back(local_delta(base_solutions[0], dataset, distance_matrix));

    //best_paths.push_back(local_candidate_moves(base_solutions[0], closest_nodes, dataset, distance_matrix));
    

    for(int i = 0; i < best_paths.size(); i++){
        best_scores.push_back(get_path_cost(best_paths[i], dataset, distance_matrix));
        worst_scores.push_back(get_path_cost(best_paths[i], dataset, distance_matrix));
        average_scores.push_back(get_path_cost(best_paths[i], dataset, distance_matrix));
    }

    for(int j = 0; j < algorithm_names.size(); j++){
        std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
        
        for(int i = 1; i < iterations; i++){
            //std::cout << i << "\n";
            
            std::vector <int> solution;
            int cost;
            if(algorithm_names[j] == "Random Baseline"){
                solution = base_solutions[i];
            }
            if(algorithm_names[j] == "Nearest"){
                solution = nearest_solution(i, dataset, distance_matrix);
            }
            if(algorithm_names[j] == "Greedy Cycle"){
                solution = greedy_cycle_solution(i, dataset, distance_matrix);
            }
            if(algorithm_names[j] == "Greedy Regret"){
                solution = greedy_2_regret_solution(i, dataset, distance_matrix);
            }
            if(algorithm_names[j] == "Greedy Weighted"){
                solution = greedy_weighted_solution(i, dataset, distance_matrix);
            }
            if(algorithm_names[j] == "Local Candidate"){
                solution = local_candidate_moves(base_solutions[i], closest_nodes, dataset, distance_matrix);
            }
            if(algorithm_names[j] == "Local Delta"){
                solution = local_delta(base_solutions[i], dataset, distance_matrix);
            }
            if(algorithm_names[j].substr(0, 12) == "Local Search"){
                bool steepest;
                bool edges;
                if(int(algorithm_names[j][12]) - 48){
                    steepest = true;
                }
                else{
                    steepest = false;
                }
                if(int(algorithm_names[j][13]) - 48){
                    edges = true;
                }
                else{
                    edges = false;
                }
                if(int(algorithm_names[j][14]) - 48){
                    solution = local_search(steepest, edges, greedy_weighted_solution(i, dataset, distance_matrix), dataset, distance_matrix);
                }
                else{
                    solution = local_search(steepest, edges, base_solutions[i], dataset, distance_matrix);
                }
                
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
            //ofs << "Best score:\n" << best_scores[j] << "\n";
            ofs << "Best, Average, Worst scores:\n" << best_scores[j] << " " << average_scores[j] << " " << worst_scores[j] << "\n";
            for(int k = 0; k < best_paths[j].size() - 1; k++){
                ofs << best_paths[j][k] << ", ";
            }
            ofs << best_paths[j][best_paths[j].size() - 1] << "\n";
            ofs << "Time: " << times[j] << "ms\n";
            //ofs << "Worst score:\n" << worst_scores[j] << "\n";
            //ofs << "Average score:\n" << average_scores[j] << "\n";
        }
        ofs.close();
    }
}

void global_convexity(std::vector <std::vector <int>> & dataset, std::vector <std::vector <int>> & distance_matrix, std::string filename = "", std::string dataset_name = "example"){
    int iterations = 1000; // 1000

    std::vector <std::vector <int>> solutions;

    std::vector <int> best_solution;
    int best_length = 10000000;

    std::vector <std::vector <int>> edge_and_node_count; // At [i][i] there is a node count of i; At [i][j] there is edge count i->j
    std::vector <int> empty_vector;
    for(int i = 0; i < dataset.size(); i++){
        empty_vector.push_back(0);
    }
    for(int i = 0; i < dataset.size(); i++){
        edge_and_node_count.push_back(empty_vector);
    }

    for(int i = 0; i < iterations; i++){
        std::vector <int> solution = random_solution(dataset, distance_matrix);
        solution = local_search(false, true, solution, dataset, distance_matrix);
        solutions.push_back(solution);

        int cost = get_path_cost(solution, dataset, distance_matrix);
        if(cost < best_length){
            best_length = cost;
            best_solution = solution;
        }

        for(int j = 0; j < solution.size(); j++){
            int node = solution[j];
            edge_and_node_count[node][node]++;
            int next_node = solution[(j+1)%solution.size()];
            edge_and_node_count[node][next_node]++;
        }
    }

    std::vector <std::vector <int>> edge_and_node_count_in_best;
    for(int i = 0; i < dataset.size(); i++){
        edge_and_node_count_in_best.push_back(empty_vector);
    }
    for(int j = 0; j < best_solution.size(); j++){
        int node = best_solution[j];
        edge_and_node_count_in_best[node][node]++;
        int next_node = best_solution[(j+1)%best_solution.size()];
        edge_and_node_count_in_best[node][next_node]++;
    }

    std::vector <float> node_best_sim;
    std::vector <float> edge_best_sim;
    std::vector <float> node_avg_sim;
    std::vector <float> edge_avg_sim;
    for(int i = 0; i < iterations; i++){
        std::vector <int> solution = solutions[i];
        float node_best = 0;
        float edge_best = 0;
        float node_avg = 0;
        float edge_avg = 0;

        for(int j = 0; j < solution.size(); j++){
            int node = solution[j];
            int next_node = solution[(j+1)%solution.size()];

            node_best += edge_and_node_count_in_best[node][node]; // 0 or 1
            edge_best += edge_and_node_count_in_best[node][next_node]; // 0 or 1

            node_avg += (edge_and_node_count[node][node] - 1); // Occurences except its own
            edge_avg += (edge_and_node_count[node][next_node] - 1); // Occurences except its own
        }

        node_best = node_best / solution.size();
        edge_best = edge_best / solution.size();

        node_avg = node_avg / (iterations-1);
        edge_avg = edge_avg / (iterations-1);
        node_avg = node_avg / solution.size();
        edge_avg = edge_avg / solution.size();

        node_best_sim.push_back(node_best);
        edge_best_sim.push_back(edge_best);
        node_avg_sim.push_back(node_avg);
        edge_avg_sim.push_back(edge_avg);
    }

    if(filename != ""){
        std::ofstream ofs;
        ofs.open(filename, std::ios_base::app);
        ofs << "\nDATASET " << dataset_name << "\n";
        ofs << "node_best edge_best node_avg edge_avg\n";
        for(int i = 0; i < iterations; i++){
            for(int j = 0; j < solutions[i].size(); j++){
                ofs << solutions[i][j] << " ";
            }
            ofs << "\n" << get_path_cost(solutions[i], dataset, distance_matrix) << "\n";
            ofs << node_best_sim[i] << " " << edge_best_sim[i] << " " << node_avg_sim[i] << " " << edge_avg_sim[i] << "\n";
        }
        ofs.close();
    }
}

void run_hybrid(int elite, int max_time, bool ls, std::vector <std::vector <int>> & dataset, std::vector <std::vector <int>> & distance_matrix, std::string filename = "", std::string dataset_name = "example"){

    //std::vector <std::vector <int>> solutions;
    //std::vector <int> costs;

    std::vector <std::pair <int, std::vector <int>>> solutions;

    int i = 0;
    while(i < elite){
        std::vector <int> solution = random_solution(dataset, distance_matrix);
        int cost = get_path_cost(solution, dataset, distance_matrix);
        std::pair <int, std::vector <int>> entry;
        entry.first = cost;
        entry.second = solution;
        // Find if already in solutions (by repeated cost):
        int j = 0;
        for(; j < solutions.size(); j++){
            if(solutions[j].first == entry.first) break;
        }
        if(j < solutions.size()) continue;

        solutions.push_back(entry);
        i++;
    }
    //std::cout << solutions.size() << "\n";
    //std::cout << solutions[0].first << "\n";
    //std::cout << solutions[solutions.size()-1].first << "\n";

    // SORT solutions
    
    std::sort(solutions.begin(), solutions.end());
    //std::cout << solutions.size() << "\n";
    //std::cout << solutions[0].first << "\n";
    //std::cout << solutions[solutions.size()-1].first << "\n";


    std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
    std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
    std::random_device dev;
    std::mt19937 rng(dev());
    std::uniform_int_distribution<std::mt19937::result_type> elite_dist(0, elite-1);
    std::uniform_int_distribution<std::mt19937::result_type> index_dist(0, solutions[0].second.size()-1);


    while(std::chrono::duration_cast<std::chrono::seconds>(end - begin).count() < max_time){
        //std::cout << std::chrono::duration_cast<std::chrono::seconds>(end - begin).count() << "\n";
        

        int parent_1_index = elite_dist(rng);
        int parent_2_index = elite_dist(rng);
        while(parent_2_index == parent_1_index) parent_2_index = elite_dist(rng);

        

        std::vector <int> parent_1 = solutions[parent_1_index].second;
        std::vector <int> parent_2 = solutions[parent_2_index].second;
        std::vector <int> child;

        // perform crossover

        // For each node in parent1:
        // If not in parent2, skip both
        // If in parent2, check if previos is the same
        // If yes, add to child in this order
        // If not, leave to add to child later

        // For now, just add nodes existing in both in order of parent1

        for(int k = 0; k < parent_1.size(); k++){
            int j = 0;
            for(; j < parent_2.size(); j++){
                if(parent_2[j] == parent_1[k]){
                    child.push_back(parent_1[k]);
                    break;
                }
            }
        }

        // preform LS
        
        while(child.size() < parent_1.size()){
            int node = index_dist(rng);
            int j = 0;
            for(; j < child.size(); j++){
                if(child[j] == node){
                    break;
                }
            }
            if(j < child.size()) continue;
            child.push_back(node);
        }

        // Random offset rotation???
        if(ls){
            child = local_search(true, true, child, dataset, distance_matrix);
        }
        
       
        // calculate cost
        int child_cost = get_path_cost(child, dataset, distance_matrix);

        // Do not put the child identical to an existing solution
        int j = 0;
        for(; j < solutions.size(); j++){
            if(child_cost == solutions[j].first){
                break;
            }
        }
        if(j < solutions.size()) continue;

        // replace the worst if better than it
        if(solutions[solutions.size()-1].first > child_cost){
            solutions[solutions.size()-1].first = child_cost;
            solutions[solutions.size()-1].second = child;

            // SORT solutions
            
            std::sort(solutions.begin(), solutions.end());
        }


        end = std::chrono::steady_clock::now();
    }

    //std::cout << "Best solution: " << solutions[0].first << "\n";
    //for(int j = 0; j < solutions[0].second.size(); j++){
    //   std::cout << solutions[0].second[j] << " ";
    //}
    //std::cout << "\n";


    if(filename != ""){
        std::ofstream ofs;
        ofs.open(filename, std::ios_base::app);
        //ofs << "\nDATASET " << dataset_name << "\n";
        ofs << "Elite: " << elite << "  Max time: " << max_time << "\n";
        std::cout << "Best solution: " << solutions[0].first << "\n";
        ofs << "Best solution: " << solutions[0].first << "\n";
        for(int j = 0; j < solutions[0].second.size(); j++){
            ofs << solutions[0].second[j] << " ";
        }
        ofs << "\n";
        ofs.close();
    }
}


int main(){
    std::srand(148253);

    std::string filename = "hybrid.txt";

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


    for(int i = 0; i < dataset_paths.size(); i++){
        std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();

        std::vector <std::vector <int>> dataset = read_dataset(dataset_paths[i]);
        std::vector <std::vector <int>> distance_matrix = create_distance_matrix(dataset);

        std::string dataset_name;
        dataset_name += (char(65 + i));
        //calculate_best_paths(dataset, distance_matrix, filename, dataset_name);
        //compare_MSLS_LSN(dataset, distance_matrix, filename, dataset_name);
        //global_convexity(dataset, distance_matrix, filename, dataset_name);
        std::ofstream ofs;
        ofs.open(filename, std::ios_base::app);
        ofs << "\nDATASET " << dataset_name << "\n";
        ofs.close();
        run_hybrid(20, 30, true, dataset, distance_matrix, filename, dataset_name);
        run_hybrid(60, 30, true, dataset, distance_matrix, filename, dataset_name);
        run_hybrid(100, 30, true, dataset, distance_matrix, filename, dataset_name);
        run_hybrid(20, 90, true, dataset, distance_matrix, filename, dataset_name);
        run_hybrid(60, 90, true, dataset, distance_matrix, filename, dataset_name);
        run_hybrid(100, 90, true, dataset, distance_matrix, filename, dataset_name);

        std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();

    }

    return 0;
}


