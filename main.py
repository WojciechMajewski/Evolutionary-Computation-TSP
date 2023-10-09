### TJWP - Travelling Jehovah's Witness Problem

### The dataset contains 3 columns - x, y coordinates and extra weight of the node, with each row representing a node


import csv

dataset_path = ""

with open('TSPA.csv', newline='') as f:
    reader = csv.reader(f, delimiter=";")
    
    data = [[int(row[0]), int(row[1]), int(row[2])] for row in reader if row]

print(data)

