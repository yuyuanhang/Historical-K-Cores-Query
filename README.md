# PHC
## Introduction
This project includes:

(1) two query algorithms for solving historical k-core problem in temporal graphs.

- a core decomposition algorithm ```Online-Query```.
- a PHC index-based query algorithm ```PHC-Query```.

(2) two PHC index construction algorithms.

- an algorithm ```PHC-Construct``` based on core neighbors.
- an optimized algorithm ```PHC-Construct^*``` based on core time neighbors.

(3) two PHC index update algorithms for edge insertion.

- an algorithm ```PHC-Update``` that directly updates earliest invalid time for each vertex.
- an optimized algorithm ```PHC-Update^*``` that updates earliest invalid time for candidates in a sub-core.

(4) an PHC index maintenance algorithm or edge expiration.

## Usage
### Compile program
```zsh
mkdir build
cd build
cmake ..
make
```

### Construction
To invoke ```PHC-Construct```, the following command line is required:
```zsh
./build/span_core -idx-bl/-idx-bl-4col [dataset] [index] [num_edge]
```
- The first parameter is set according to the format of dataset.
Specifically, when the format of dataset is ```u, v, t```, use ```-idx-bl```; When the format of dataset is ```u, v, w, t```, use ```-idx-bl-4col```.

- The second parameter is the path of dataset, only ```*.txt``` file is allowed.

- The third parameter is the path of index, which is stored in binary format.

- The fourth parameter is optional, which represents the number of edges to be loaded.
Specifically, if this parameter is not specified, all edges in the graph will be loaded.
If the parameter is specified, edges in the graph are sorted in chronological order, and only the first num_edge edges will be loaded.

To invoke ```PHC-Construct^*```, the following command line is required:
```zsh
./build/span_core -idx/-idx-4col [dataset] [index] [num_edge]
```
Here, parameters are the same as described above, the only difference is that the first parameter is ```-idx/-idx-4col```.

### Search
To invoke query algorithms, the following command line is required:
```zsh
./build/span_core -q-sd/-q-sd-4col [dataset] [index] [log] [t_range] [k_range] [num_query]
```
- Similarly, the first parameter is set according to the format of dataset.

- The second parameter is the path of dataset.

- The third parameter is the path of index, which is loaded into main memory.

- The fourth parameter is the path of log file, which records experimental results such as average running time and standard deviation.

- The fifth parameter is the length of generated time windows, i.e., |t_e - t_s|. ```t_range``` lies in [1, 100], which sets the length of generated time windows to ```t_range/100*t_max```.

- The sixth parameter specifies the value of ```k```. ```k_range``` lies in [1, 10], which sets ```k``` to ```k_range/10*k_max```.

- The seventh parameter is the number of generated time windows.

## Datasets
The information of used real-world datasets is provided in our paper. Currently, we support ```.fvecs``` file for base datasets and query datasets 
and ```.ivecs``` for groundtruth file.
