#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <string>
#include <cstdlib>
#include <algorithm>
#include <chrono>
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

std::vector <int> local_search(bool steepest_neighborhood, bool edges_exchange, int starting_node, std::vector <std::vector <int>> & dataset, std::vector <std::vector <int>> & distance_matrix){
    // if not steepest_neighborhood, then greedy
    // if not edges_exchange, then nodes exchange // But intra- and inter- should both be used
    // if starting_node == -1 then random starting solution

    std::vector <int> solution;
    if(starting_node != -1){
        solution = greedy_weighted_solution(starting_node, dataset, distance_matrix);
        //std::cout << get_path_cost(solution, dataset, distance_matrix) << "\n";
    }
    else{
        solution = random_solution(dataset, distance_matrix);
    }

    // sort the starting solution, then for each of 200 check if in sorted, if yes advance sorted, if no add to available
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

    // 1000 iterations instead of while(true)
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

                        cost_improvement -= get_cost(solution[i], solution[(j+1) % solution.size()], dataset, distance_matrix);
                        cost_improvement -= get_cost(solution[j], solution[(i+1) % solution.size()], dataset, distance_matrix);

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
                //std::cout << "improvement by: " << best_improvement;
                //std::cout << "\nold: " << get_path_cost(solution, dataset, distance_matrix);
                if(new_node_index != -1){ // Then new node is best
                
                    int temp = solution[replacing_node_index];
                    solution[replacing_node_index] = available_nodes[new_node_index];
                    available_nodes[new_node_index] = temp;
                }
                else{
                    if(edges_exchange){ // Then edge exchange worked
                        // Edge exchange through flipping a subpath
                        std::reverse(solution.begin() + ((first_edge_start_index + 1) % solution.size()), solution.begin() + ((second_edge_start_index + 1) % solution.size()));
                    }
                    else{ // Then node exchange worked
                    
                        //std::cout << ", by replacing nodes";
                        int temp = solution[first_to_switch_index];
                        solution[first_to_switch_index] = solution[second_to_switch_index];
                        solution[second_to_switch_index] = temp;
                    }
                }
                //std::cout << "\nnew: " << get_path_cost(solution, dataset, distance_matrix) << "\n";
            }
            else{ // No improvement, break the local search loop, local optimum found!
                return solution;
            }
        }
        else{
            int random_offset = 0; //rand() % solution.size(); TODO FIX
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
                int choice = 0; //rand() % 2; // TODO FIX

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
                        if((j_intra+1) % solution.size() >= i_intra){
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

                        cost_improvement -= get_cost(solution[i_intra], solution[(j_intra+1) % solution.size()], dataset, distance_matrix);
                        cost_improvement -= get_cost(solution[j_intra], solution[(i_intra+1) % solution.size()], dataset, distance_matrix);

                        cost_improvement += old_cost;
                        if(cost_improvement > 0){
                            std::reverse(solution.begin() + ((j_intra + 1) % solution.size()), solution.begin() + ((i_intra + 1) % solution.size()));
                            break;
                        }
                        
                    }

                    // Node rearrangement advance
                    else{
                        //std::cout << i_intra_iter << " " << j_intra << " " << i_intra << "\n";
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
                            std::cout << "Found improvement" << "\n";
                            break;
                        }
                    }
                }
                else{
                    
            
                    //std::cout << "hereere" << "\n";
                    return solution;
                    // No improvement can be found by greedy, so end
                }
            }
        }
    }

    return solution;
}

void calculate_best_paths(std::vector <std::vector <int>> & dataset, std::vector <std::vector <int>> & distance_matrix, std::string filename = "", std::string dataset_name = "example"){
    int iterations = 5;//200;
    std::vector <std::vector <int>> best_paths;
    std::vector <int> best_scores;
    std::vector <int> worst_scores;
    std::vector <float> average_scores;
    std::vector <int> times;
    
    std::vector <std::string> algorithm_names;
    //algorithm_names.push_back("Nearest");
    //algorithm_names.push_back("Greedy Cycle");
    //algorithm_names.push_back("Greedy Regret");
    algorithm_names.push_back("Greedy Weighted");
    algorithm_names.push_back("Local Search00");
    //algorithm_names.push_back("Local Search10");
    //algorithm_names.push_back("Local Search01");
    //algorithm_names.push_back("Local Search11");

    //best_paths.push_back(nearest_solution(0, dataset, distance_matrix));
    //best_paths.push_back(greedy_cycle_solution(0, dataset, distance_matrix));
    //best_paths.push_back(greedy_2_regret_solution(0, dataset, distance_matrix));
    best_paths.push_back(greedy_weighted_solution(0, dataset, distance_matrix));
    best_paths.push_back(local_search(false, false, 0, dataset, distance_matrix));
    //best_paths.push_back(local_search(true, false, 0, dataset, distance_matrix));
    //best_paths.push_back(local_search(false, true, 0, dataset, distance_matrix));
    //best_paths.push_back(local_search(true, true, 0, dataset, distance_matrix));

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
                
                solution = local_search(steepest, edges, i, dataset, distance_matrix);
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
            ofs << best_paths[j][best_paths[j].size() - 1] << "\n";
            ofs << "Worst score:\n" << worst_scores[j] << "\n";
            ofs << "Average score:\n" << average_scores[j] << "\n";
            ofs << "Time:\n" << times[j] << "ms\n";
        }
        ofs.close();
    }
}


int main(){
    std::srand(148253);

    std::string filename = "cpp_local_search_results.txt";


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
        std::vector <std::vector <int>> dataset = read_dataset(dataset_paths[i]);
        std::vector <std::vector <int>> distance_matrix = create_distance_matrix(dataset);

        std::string dataset_name;
        dataset_name += (char(65 + i));
        calculate_best_paths(dataset, distance_matrix, filename, dataset_name);

        
    }
    return 0;
}


