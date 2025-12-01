# Round-Robin Scheduler

A C++ program that finds an optimized round-robin schedule for a free-for-all game with 4 players per game, using simulated annealing to minimize player pair repeats.

## Problem Description

- **N players** (must be a multiple of 4)
- **K rounds** of games
- Each game has **4 players**
- Goal: Minimize the number of times any two players play together more than once

## Default Configuration

- N = 32 players
- K = 10 rounds
- 8 games per round

## Compilation

### Using Make (Linux/Mac/Windows with MinGW/MSYS):
```bash
make
```

To clean build artifacts:
```bash
make clean
```

### Using g++ directly:
```bash
g++ -std=c++11 -O2 -o round_robin round_robin.cpp
```

### Using Visual Studio (Windows):
```bash
cl /EHsc /O2 round_robin.cpp
```

## Usage

```bash
./round_robin
```

The program will:
1. Run in an infinite loop, attempting to find better schedules
2. For each attempt:
   - Generate a random initial schedule
   - Use simulated annealing to optimize the schedule
   - Display statistics about pair repeats
   - Print the final optimized schedule
   - Save results to `results/` directory if the score is better than existing results
3. Stop automatically if an ideal schedule (cost = 0) is found
4. Press Ctrl+C to stop manually

Results are saved in the `results/` directory (created automatically if it doesn't exist) with the format: `N_K_score_time.txt`
- Example: `32_10_19_1764621860.txt` (32 players, 10 rounds, score 19, Unix timestamp)
- Only results with better (lower) scores than existing files are saved

## Algorithm

The program uses **Simulated Annealing** to optimize the schedule:

1. **Initial State**: Random schedule where all players are randomly assigned to games
2. **Neighbor Operation**: Swap two random players from different games within the same round (maintains constraint that each player appears exactly once per round)
3. **Cost Function**: Sum of (count - 1) for all pairs that appear more than once
4. **Acceptance**: Accept better solutions always, accept worse solutions with probability based on temperature
5. **Cooling**: Gradually decrease temperature to converge to a good solution

The optimization runs continuously, saving only results that improve upon the best score found so far.

## Output

The program outputs to console:
- Optimization progress (iterations, temperature, cost, acceptance rate)
- Statistics (total pairs, repeating pairs, maximum repeats, repeat distribution)
- N×N pair count matrix (showing how many times each pair of players plays together)
- Complete schedule showing which players are in each game for each round (1-indexed)

Results are also saved to files in the `results/` directory, containing:
- Header information (players, rounds, attempt number, costs)
- Statistics
- Pair count matrix
- Complete schedule

## Customization

You can modify the parameters in `main()`:
- `N`: Number of players (must be multiple of 4)
- `K`: Number of rounds

You can also adjust optimization parameters in the `optimize()` call in `main()`:
- `initial_temp`: Starting temperature (currently: 1000.0)
- `cooling_rate`: Temperature reduction factor (currently: 0.9999)
- `iterations_per_temp`: Iterations at each temperature (currently: 10000)
- `min_temp`: Minimum temperature to stop (currently: 0.01)

The default parameters in the `optimize()` function signature are:
- `initial_temp`: 1000.0
- `cooling_rate`: 0.995
- `iterations_per_temp`: 1000
- `min_temp`: 0.1

