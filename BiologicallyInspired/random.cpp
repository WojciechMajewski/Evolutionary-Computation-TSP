
#include <fstream>
#include <random>
#include <stdlib.h>
#include <iostream>


int main(){
    std::srand(148253);

    std::string filename = "C:/Users/wojci/Documents/GitHub/Evolutionary-Computation-TSP/BiologicallyInspired/ATSP/br17.atsp";

    std::ifstream file(filename);

    std::vector <std::vector <int>> rows;
    std::vector <int> row;
    std::string row_string;

    while (getline(file, row_string)) {
        if(row_string.find(":", 0) != std::string::npos){
            if(row_string.find("DIMENSION", 0) != std::string::npos){
                std::cout << row_string;
            }
        }


    }
    file.close();


    return 0;
}