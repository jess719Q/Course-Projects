# hw3 Detailed Placement

## Overview
This project implements a detailed placement algorithm based on FastDP, 
aiming to reduce total wirelength (HPWL) while maintaining legal placement. 
The program parses LEF/DEF files, performs iterative optimization, and reports 
final wirelength and runtime for each testcase.

## Implementation Details

The FastDP algorithm is adapted with the following iterative flow:

```
1. Single-Segment Clustering
2. Repeat until iteration limit or timeout:
    - Shuffle component order
    - Global Swap
    - Vertical Swap
    - Local Reordering
3. Final Single-Segment Clustering
```

### Key design choices:

- Randomized traversal order is applied each iteration to improve solution diversity
- Best HPWL and placement are recorded and restored at the end
- The redundant final clustering step in original FastDP is removed, as it showed no further improvement
- Abacus legalization is applied after swap stages to ensure row legality

### Swap strategies

- Global Swap: Searches within optimal regions and nearby whitespace, using a decaying negative threshold
- Vertical Swap: Restricts candidates to nearby rows
- Local Reordering: Optimizes ordering within a sliding window of three components

## Compile and Run
  In ```src/```, enter the following command:
  
  ```
  make
  ```
  
  An executable file ```hw3``` will be generated in ```bin/```.
  
  If you want to remove it, please enter the following command:
  
  ``` 
  make clean
  ```

  Usage:
  
  ``` 
  ./hw3 <input LEF> <input DEF> <output DEF>
```
