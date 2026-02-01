# hw2 Multi-way Min-Cut Partitioning

## Overview
This project implements a multi-way min-cut partitioning algorithm based on the **Fiduccia–Mattheyses (FM) framework**, extended from the original 2-way partitioning to support **n-way (e.g., 4-way)** partitioning.  
The implementation focuses on efficient gain updates, scalable data structures, and runtime optimization through parallel execution.

## Key Features
* Extension of FM algorithm from 2-way to n-way partitioning
* Generalized bucket list using a 3D data structure
* Support for multiple initial solutions with parallel execution
* Efficient gain update and restoration mechanism
* Recursive and iterative refinement for 4-way partitioning



## Compile and Run
  In ```src/```, enter the following command:
  
  ```
  make
  ```
  
  An executable file ```hw2``` will be generated in ```bin/```.
  
  If you want to remove it, please enter the following command:
  
  ``` 
  make clean
  ```

  Usage:
  
  ``` 
  ./hw2 <input file> <output file> <number of partitions>
```
