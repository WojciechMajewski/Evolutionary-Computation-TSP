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


int main(){

    std::srand(148253);

    std::vector <std::vector <int>> dataset_A = read_dataset("Data/TSPA.csv");

    std::vector <std::vector <int>> distance_matrix_A = create_distance_matrix(dataset_A);

    std::vector <int> solution_random = random_solution(dataset_A, distance_matrix_A);

    std::cout << get_path_cost(solution_random, dataset_A, distance_matrix_A);

    //std::cout << get_cost(0, 1, dataset_A, distance_matrix_A);


    //for(int i = 0; i < distance_matrix_A.size(); i++){
    //    for(int j = 0; j < distance_matrix_A[i].size(); j++){
    //        std::cout << distance_matrix_A[i][j] << "\n";
    //    }
    //}

    // Create extra cost list

    // Greedy heuristic algorithms


    return 0;
}


