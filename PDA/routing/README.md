# hw4 Global Routing

## Overview
This project implements a global router using a negotiation-based rip-up and re-route framework, 
optimized with A* search and a dynamic cost function to minimize congestion and wirelength.

## Implementation Details

### Algorithm Overview
The router consists of three main stages:

1. **Initial Routing**  
   - Nets are sorted by HPWL and routed using A\* search.  
   - Preliminary usage information is collected and all nets are immediately re-routed for a balanced starting solution.

2. **Negotiation Loop**  
   - Congestion is resolved iteratively by penalizing overused edges using dynamic **history costs**.  
   - Nets with excessive detours or high congestion are prioritized for rip-up and re-route.  
   - Bounding box constraints gradually expand during iterations for better path exploration.

3. **Greedy Refinement**  
   - Nets exceeding their HPWL are further optimized.  
   - Full edges are treated as obstacles, historical bias is removed, ensuring wirelength reduction without introducing new overflows.

## Compile and Run
  In ```src/```, enter the following command:
  
  ```
  make
  ```
  
  An executable file ```hw4``` will be generated in ```bin/```.
  
  If you want to remove it, please enter the following command:
  
  ``` 
  make clean
  ```

  Usage:
  
  ``` 
  ./hw4 <input file> <output file>
```
