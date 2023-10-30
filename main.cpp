#include <iostream>
#include <fstream>
#include <vector>
#include <string>

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


int main(){

    std::vector <std::vector <int>> dataset_A = read_dataset("Data/TSPA.csv");

    //for(int i = 0; i < dataset_A.size(); i++){
    //    for(int j = 0; j < dataset_A[i].size(); j++){
    //        std::cout << dataset_A[i][j];
    //    }
    //}

    // Create distance matrix

    // Create extra cost list

    // Greedy heuristic algorithms


    return 0;
}


