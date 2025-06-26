# PHC
## Update
There is a minor mistake in the PHC-Update* algorithm that can lead to incorrect results, which we have now corrected in our implementation.
Specifically, in Lines 1–3 of Algorithm 14 (“UpdateEIT”), the original pseudocode initializes CN only for the vertices in the candidate set C;
The revised version initializes CN for every vertex in the neighborhood of each vertex in C whose core value in the [t_{s}, t_{max}] window is at least k.

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

After the above command line is completed, the results of ```Online-Query``` and ```PHC-Query``` are both reported.
Specifically, the results of ```PHC-Query``` are first reported, and the results of ```Online-Query``` are then reported.

### Update
To invoke ```PHC-Update```, the following command line is required:
```zsh
./build/span_core -uidx-bl/-uidx-bl-4col [dataset] [index] [log] [ratio]
```
Here, parameters are the same as described above, the only difference is that the fifth parameter specifies the ratio between arrival edges and total edges.
Specifically, ```ratio``` lies in (0, 1]. Edges over time window [0, ratio * t_max) are used to construct PHC index, and edges over time window [ratio * t_max, t_max) are used as arrival edges to update PHC index.

To invoke ```PHC-Update^*```, the following command line is required:
```zsh
./build/span_core -uidx/-uidx-4col [dataset] [index] [log] [ratio]
```
Here, parameters are the same as described above, the only difference is that the first parameter is ```-uidx/-uidx-4col```.

### Handle expired edges
To update PHC index for expired edges, the following command line is required:
```zsh
./build/span_core -rm/-rm-4col [dataset] [index] [log] [ratio]
```
Specifically, ```ratio``` lies in (0, 1]. Edges over time window [0, ratio * t_max) are regarded as expired edges to update PHC index.

## Datasets
The information of used real-world datasets is provided in our paper.
