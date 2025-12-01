#include <iostream>
#include <vector>
#include <random>
#include <algorithm>
#include <cmath>
#include <iomanip>
#include <map>
#include <chrono>
#include <fstream>
#include <sstream>
#include <string>
#ifdef _WIN32
#include <direct.h>
#include <windows.h>
#else
#include <sys/stat.h>
#include <sys/types.h>
#include <dirent.h>
#endif

class RoundRobinScheduler {
private:
    int N;  // Number of players (multiple of 4)
    int K;  // Number of rounds
    int games_per_round;
    
    // Schedule: rounds[round][game][player_index]
    std::vector<std::vector<std::vector<int>>> schedule;
    
    // Track pairs: pair_count[p1][p2] = number of times p1 and p2 play together
    std::vector<std::vector<int>> pair_count;
    
    std::mt19937 rng;
    
    // Calculate cost: sum of (count - 1) for all pairs that appear more than once
    int calculateCost() const {
        int cost = 0;
        for (int i = 0; i < N; i++) {
            for (int j = i + 1; j < N; j++) {
                if (pair_count[i][j] > 1) {
                    cost += pair_count[i][j] - 1;
                }
            }
        }
        return cost;
    }
    
    // Update pair counts based on current schedule
    void updatePairCounts() {
        // Reset pair counts
        for (int i = 0; i < N; i++) {
            std::fill(pair_count[i].begin(), pair_count[i].end(), 0);
        }
        
        // Count pairs in current schedule
        for (int round = 0; round < K; round++) {
            for (int game = 0; game < games_per_round; game++) {
                const auto& players = schedule[round][game];
                for (int i = 0; i < 4; i++) {
                    for (int j = i + 1; j < 4; j++) {
                        int p1 = players[i];
                        int p2 = players[j];
                        pair_count[p1][p2]++;
                        pair_count[p2][p1]++;
                    }
                }
            }
        }
    }
    
    // Generate random initial schedule
    void generateRandomSchedule() {
        std::vector<int> all_players(N);
        for (int i = 0; i < N; i++) {
            all_players[i] = i;
        }
        
        for (int round = 0; round < K; round++) {
            std::shuffle(all_players.begin(), all_players.end(), rng);
            
            for (int game = 0; game < games_per_round; game++) {
                for (int pos = 0; pos < 4; pos++) {
                    schedule[round][game][pos] = all_players[game * 4 + pos];
                }
            }
        }
        
        updatePairCounts();
    }
    
    // Validate that each round has exactly N unique players (each player appears once)
    bool validateSchedule() const {
        for (int round = 0; round < K; round++) {
            std::vector<bool> player_used(N, false);
            for (int game = 0; game < games_per_round; game++) {
                for (int pos = 0; pos < 4; pos++) {
                    int player = schedule[round][game][pos];
                    if (player < 0 || player >= N) {
                        return false;  // Invalid player index
                    }
                    if (player_used[player]) {
                        return false;  // Player appears twice in same round
                    }
                    player_used[player] = true;
                }
            }
            // Check that all players are used
            for (int i = 0; i < N; i++) {
                if (!player_used[i]) {
                    return false;  // Player missing from round
                }
            }
        }
        return true;
    }
    
    // Swap two players in the schedule (neighbor operation)
    void swapPlayers(int round1, int game1, int pos1, int round2, int game2, int pos2) {
        std::swap(schedule[round1][game1][pos1], schedule[round2][game2][pos2]);
    }
    
    // Get a random neighbor by swapping two random players within the same round
    // This ensures that each player appears exactly once per round
    // Swaps are only between different games (swapping within same game doesn't change pairs)
    void getRandomNeighbor(int& round1, int& game1, int& pos1, int& round2, int& game2, int& pos2) {
        round1 = rng() % K;
        round2 = round1;  // Same round to maintain constraint
        
        // Get two different games in the same round
        do {
            game1 = rng() % games_per_round;
            game2 = rng() % games_per_round;
        } while (game1 == game2);  // Ensure different games
        
        // Random positions within each game
        pos1 = rng() % 4;
        pos2 = rng() % 4;
    }
    
public:
    RoundRobinScheduler(int n, int k, unsigned int seed = 0) 
        : N(n), K(k), games_per_round(n / 4), rng(seed) {
        
        if (N % 4 != 0) {
            throw std::invalid_argument("N must be a multiple of 4");
        }
        
        // Initialize schedule
        schedule.resize(K);
        for (int round = 0; round < K; round++) {
            schedule[round].resize(games_per_round);
            for (int game = 0; game < games_per_round; game++) {
                schedule[round][game].resize(4);
            }
        }
        
        // Initialize pair count matrix
        pair_count.resize(N);
        for (int i = 0; i < N; i++) {
            pair_count[i].resize(N, 0);
        }
        
        generateRandomSchedule();
        
        // Validate initial schedule
        if (!validateSchedule()) {
            throw std::runtime_error("Initial schedule validation failed!");
        }
    }
    
    // Simulated annealing optimization
    void optimize(double initial_temp = 1000.0, double cooling_rate = 0.995, 
                  int iterations_per_temp = 1000, double min_temp = 0.1) {
        
        int current_cost = calculateCost();
        int best_cost = current_cost;
        auto best_schedule = schedule;
        
        double temperature = initial_temp;
        long long iteration = 0;
        
        std::cout << "Starting optimization...\n";
        std::cout << "Initial cost: " << current_cost << "\n\n";
        
        while (temperature > min_temp) {
            int accepted = 0;
            
            for (int i = 0; i < iterations_per_temp; i++) {
                // Get random neighbor
                int round1, game1, pos1, round2, game2, pos2;
                getRandomNeighbor(round1, game1, pos1, round2, game2, pos2);
                
                // Perform swap
                int p1 = schedule[round1][game1][pos1];
                int p2 = schedule[round2][game2][pos2];
                
                // Calculate cost change
                int cost_change = 0;
                
                // Remove old pairs
                for (int pos = 0; pos < 4; pos++) {
                    if (pos != pos1) {
                        int other = schedule[round1][game1][pos];
                        if (pair_count[p1][other] > 1) cost_change -= 1;
                        pair_count[p1][other]--;
                        pair_count[other][p1]--;
                    }
                }
                for (int pos = 0; pos < 4; pos++) {
                    if (pos != pos2) {
                        int other = schedule[round2][game2][pos];
                        if (pair_count[p2][other] > 1) cost_change -= 1;
                        pair_count[p2][other]--;
                        pair_count[other][p2]--;
                    }
                }
                
                // Perform swap
                swapPlayers(round1, game1, pos1, round2, game2, pos2);
                
                // Add new pairs
                for (int pos = 0; pos < 4; pos++) {
                    if (pos != pos1) {
                        int other = schedule[round1][game1][pos];
                        pair_count[p2][other]++;
                        pair_count[other][p2]++;
                        if (pair_count[p2][other] > 1) cost_change += 1;
                    }
                }
                for (int pos = 0; pos < 4; pos++) {
                    if (pos != pos2) {
                        int other = schedule[round2][game2][pos];
                        pair_count[p1][other]++;
                        pair_count[other][p1]++;
                        if (pair_count[p1][other] > 1) cost_change += 1;
                    }
                }
                
                int new_cost = current_cost + cost_change;
                
                // Accept or reject
                bool accept = false;
                if (cost_change <= 0) {
                    accept = true;
                } else {
                    double prob = std::exp(-cost_change / temperature);
                    if (std::uniform_real_distribution<double>(0.0, 1.0)(rng) < prob) {
                        accept = true;
                    }
                }
                
                if (accept) {
                    current_cost = new_cost;
                    accepted++;
                    
                    if (current_cost < best_cost) {
                        best_cost = current_cost;
                        best_schedule = schedule;
                    }
                } else {
                    // Revert swap
                    swapPlayers(round1, game1, pos1, round2, game2, pos2);
                    
                    // Revert pair counts
                    for (int pos = 0; pos < 4; pos++) {
                        if (pos != pos1) {
                            int other = schedule[round1][game1][pos];
                            pair_count[p2][other]--;
                            pair_count[other][p2]--;
                        }
                    }
                    for (int pos = 0; pos < 4; pos++) {
                        if (pos != pos2) {
                            int other = schedule[round2][game2][pos];
                            pair_count[p1][other]--;
                            pair_count[other][p1]--;
                        }
                    }
                    
                    // Restore old pairs
                    for (int pos = 0; pos < 4; pos++) {
                        if (pos != pos1) {
                            int other = schedule[round1][game1][pos];
                            pair_count[p1][other]++;
                            pair_count[other][p1]++;
                        }
                    }
                    for (int pos = 0; pos < 4; pos++) {
                        if (pos != pos2) {
                            int other = schedule[round2][game2][pos];
                            pair_count[p2][other]++;
                            pair_count[other][p2]++;
                        }
                    }
                }
                
                iteration++;
            }
            
            // Print progress - more frequent at low temperature, less frequent at high temperature
            long long log_interval;
            if (temperature > initial_temp * 0.5) {
                // High temperature: log every 5000 iterations per temp
                log_interval = static_cast<long long>(iterations_per_temp) * 5000;
            } else if (temperature > initial_temp * 0.1) {
                // Medium temperature: log every 2000 iterations per temp
                log_interval = static_cast<long long>(iterations_per_temp) * 2000;
            } else if (temperature > initial_temp * 0.01) {
                // Low temperature: log every 1000 iterations per temp
                log_interval = static_cast<long long>(iterations_per_temp) * 1000;
            } else {
                // Very low temperature: log every 500 iterations per temp
                log_interval = static_cast<long long>(iterations_per_temp) * 500;
            }
            
            if (iteration % log_interval == 0 || temperature <= min_temp * 1.01) {
                std::cout << "Iteration: " << iteration 
                          << ", Temperature: " << std::fixed << std::setprecision(3) << temperature
                          << ", Current cost: " << current_cost
                          << ", Best cost: " << best_cost
                          << ", Acceptance rate: " << std::setprecision(2) 
                          << (100.0 * accepted / iterations_per_temp) << "%\n";
            }
            
            temperature *= cooling_rate;
        }
        
        // Restore best schedule
        schedule = best_schedule;
        updatePairCounts();
        
        // Validate final schedule
        if (!validateSchedule()) {
            std::cerr << "WARNING: Final schedule validation failed!\n";
        }
        
        std::cout << "\nOptimization complete!\n";
        std::cout << "Final cost: " << calculateCost() << "\n";
    }
    
    // Print statistics
    void printStatistics() const {
        std::cout << "\n=== STATISTICS ===\n\n";
        
        int total_pairs = 0;
        int repeat_pairs = 0;
        int max_repeats = 0;
        std::map<int, int> repeat_distribution;
        
        for (int i = 0; i < N; i++) {
            for (int j = i + 1; j < N; j++) {
                int count = pair_count[i][j];
                total_pairs++;
                if (count > 1) {
                    repeat_pairs++;
                    max_repeats = std::max(max_repeats, count);
                    repeat_distribution[count]++;
                }
            }
        }
        
        std::cout << "Total unique pairs: " << total_pairs << "\n";
        std::cout << "Pairs that repeat: " << repeat_pairs << "\n";
        std::cout << "Maximum repeats: " << max_repeats << "\n";
        std::cout << "\nRepeat distribution:\n";
        for (const auto& entry : repeat_distribution) {
            std::cout << "  " << entry.first << " times: " << entry.second << " pairs\n";
        }
        
        std::cout << "\nCost (sum of excess repeats): " << calculateCost() << "\n";
    }
    
    // Get current cost
    int getCost() const {
        return calculateCost();
    }
    
    // Print NxN matrix of pair counts
    void printPairMatrix() const {
        std::cout << "\n=== PAIR COUNT MATRIX ===\n\n";
        for (int i = 0; i < N; i++) {
            for (int j = 0; j < N; j++) {
                std::cout << pair_count[i][j];
            }
            std::cout << "\n";
        }
    }
    
    // Print schedule to stream
    void printScheduleToStream(std::ostream& os) const {
        os << "\n=== SCHEDULE ===\n\n";
        for (int round = 0; round < K; round++) {
            os << "Round " << (round + 1) << ":\n";
            for (int game = 0; game < games_per_round; game++) {
                os << "  Game " << (game + 1) << ": ";
                for (int pos = 0; pos < 4; pos++) {
                    os << std::setw(3) << (schedule[round][game][pos] + 1);  // 1-indexed
                    if (pos < 3) os << ", ";
                }
                os << "\n";
            }
            os << "\n";
        }
    }
    
    // Print statistics to stream
    void printStatisticsToStream(std::ostream& os) const {
        os << "\n=== STATISTICS ===\n\n";
        
        int total_pairs = 0;
        int repeat_pairs = 0;
        int max_repeats = 0;
        std::map<int, int> repeat_distribution;
        
        for (int i = 0; i < N; i++) {
            for (int j = i + 1; j < N; j++) {
                int count = pair_count[i][j];
                total_pairs++;
                if (count > 1) {
                    repeat_pairs++;
                    max_repeats = std::max(max_repeats, count);
                    repeat_distribution[count]++;
                }
            }
        }
        
        os << "Total unique pairs: " << total_pairs << "\n";
        os << "Pairs that repeat: " << repeat_pairs << "\n";
        os << "Maximum repeats: " << max_repeats << "\n";
        os << "\nRepeat distribution:\n";
        for (const auto& entry : repeat_distribution) {
            os << "  " << entry.first << " times: " << entry.second << " pairs\n";
        }
        
        os << "\nCost (sum of excess repeats): " << calculateCost() << "\n";
    }
    
    // Print pair matrix to stream
    void printPairMatrixToStream(std::ostream& os) const {
        os << "\n=== PAIR COUNT MATRIX ===\n\n";
        for (int i = 0; i < N; i++) {
            for (int j = 0; j < N; j++) {
                os << pair_count[i][j];
            }
            os << "\n";
        }
    }
};

// Create directory if it doesn't exist
void createDirectory(const std::string& dir) {
#ifdef _WIN32
    _mkdir(dir.c_str());
#else
    mkdir(dir.c_str(), 0755);
#endif
}

// Find the best (lowest) score in the results directory
// Returns the best score, or -1 if no results found
// Filename format: N_K_score_time.txt
int findBestScore(const std::string& dir) {
    int best_score = -1;
    
#ifdef _WIN32
    WIN32_FIND_DATAA findData;
    std::string searchPath = dir + "\\*.txt";
    HANDLE hFind = FindFirstFileA(searchPath.c_str(), &findData);
    
    if (hFind != INVALID_HANDLE_VALUE) {
        do {
            std::string filename = findData.cFileName;
            // Parse score from filename: N_K_score_time.txt
            // Find the third underscore (after N, K, score)
            size_t pos1 = filename.find('_');
            if (pos1 != std::string::npos) {
                size_t pos2 = filename.find('_', pos1 + 1);
                if (pos2 != std::string::npos) {
                    size_t pos3 = filename.find('_', pos2 + 1);
                    if (pos3 != std::string::npos) {
                        std::string score_str = filename.substr(pos2 + 1, pos3 - pos2 - 1);
                        try {
                            int score = std::stoi(score_str);
                            if (best_score == -1 || score < best_score) {
                                best_score = score;
                            }
                        } catch (...) {
                            // Ignore files that don't match the pattern
                        }
                    }
                }
            }
        } while (FindNextFileA(hFind, &findData));
        FindClose(hFind);
    }
#else
    DIR* dirp = opendir(dir.c_str());
    if (dirp != nullptr) {
        struct dirent* entry;
        while ((entry = readdir(dirp)) != nullptr) {
            std::string filename = entry->d_name;
            // Check if it's a .txt file
            if (filename.length() > 4 && filename.substr(filename.length() - 4) == ".txt") {
                // Parse score from filename: N_K_score_time.txt
                // Find the third underscore (after N, K, score)
                size_t pos1 = filename.find('_');
                if (pos1 != std::string::npos) {
                    size_t pos2 = filename.find('_', pos1 + 1);
                    if (pos2 != std::string::npos) {
                        size_t pos3 = filename.find('_', pos2 + 1);
                        if (pos3 != std::string::npos) {
                            std::string score_str = filename.substr(pos2 + 1, pos3 - pos2 - 1);
                            try {
                                int score = std::stoi(score_str);
                                if (best_score == -1 || score < best_score) {
                                    best_score = score;
                                }
                            } catch (...) {
                                // Ignore files that don't match the pattern
                            }
                        }
                    }
                }
            }
        }
        closedir(dirp);
    }
#endif
    
    return best_score;
}

void printUsage(const char* program_name) {
    std::cerr << "Usage: " << program_name << " <N> <K>\n";
    std::cerr << "  N: Number of players (must be a multiple of 4)\n";
    std::cerr << "  K: Number of rounds\n";
    std::cerr << "\nExample: " << program_name << " 32 10\n";
}

int main(int argc, char* argv[]) {
    // Parse command line arguments
    if (argc != 3) {
        printUsage(argv[0]);
        return 1;
    }
    
    int N, K;
    try {
        N = std::stoi(argv[1]);
        K = std::stoi(argv[2]);
    } catch (const std::exception& e) {
        std::cerr << "Error: Invalid arguments. N and K must be integers.\n\n";
        printUsage(argv[0]);
        return 1;
    }
    
    // Validate N (must be multiple of 4)
    if (N <= 0 || N % 4 != 0) {
        std::cerr << "Error: N must be a positive multiple of 4.\n";
        std::cerr << "  Provided: N = " << N << "\n\n";
        printUsage(argv[0]);
        return 1;
    }
    
    // Validate K (must be positive)
    if (K <= 0) {
        std::cerr << "Error: K must be a positive integer.\n";
        std::cerr << "  Provided: K = " << K << "\n\n";
        printUsage(argv[0]);
        return 1;
    }
    
    std::cout << "Round-Robin Scheduler\n";
    std::cout << "====================\n";
    std::cout << "Players: " << N << "\n";
    std::cout << "Rounds: " << K << "\n";
    std::cout << "Games per round: " << (N / 4) << "\n";
    std::cout << "Running extended optimization in infinite loop... (Press Ctrl+C to stop)\n\n";
    
    // Create results directory
    createDirectory("results");
    
    int attempt = 0;
    
    while (true) {
        attempt++;
        std::cout << "\n=== Attempt " << attempt << " ===\n";
        
        // Use current time as seed (different for each attempt)
        unsigned int seed = std::chrono::system_clock::now().time_since_epoch().count() + attempt;
        RoundRobinScheduler scheduler(N, K, seed);
        
        int initial_cost = scheduler.getCost();
        std::cout << "Initial cost: " << initial_cost << "\n";
        
        // Check best score before optimization
        int best_score_before = findBestScore("results");
        if (best_score_before == -1) {
            std::cout << "No previous results found.\n";
        } else {
            std::cout << "Best score in directory: " << best_score_before << "\n";
        }
        std::cout << "\n";
        
        // Start timing
        auto start_time = std::chrono::high_resolution_clock::now();
        
        // Run extended optimization with faster parameters
        scheduler.optimize(
            1000.0,     // initial temperature (higher for more exploration)
            0.9999,     // cooling rate (faster cooling)
            10000,      // iterations per temperature (fewer iterations at each temp)
            0.01        // minimum temperature (higher minimum for faster completion)
        );
        
        // End timing
        auto end_time = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time);
        
        int final_cost = scheduler.getCost();
        std::cout << "\nFinal cost: " << final_cost << "\n";
        std::cout << "Optimization time: " << duration.count() << " ms ("
                  << std::fixed << std::setprecision(2) << (duration.count() / 1000.0) << " seconds)\n";
        
        // Check if this result is better than the best in the directory
        int best_score = findBestScore("results");
        if (best_score == -1) {
            std::cout << "No previous results found.\n";
        } else {
            std::cout << "Best score in directory: " << best_score << "\n";
        }
        
        // Only save if this is better (lower cost) than the best, or if no results exist
        bool should_save = (best_score == -1 || final_cost < best_score);
        std::string saved_filename = "";
        
        if (should_save) {
            // Generate filename: N_K_score_time.txt (time is Unix timestamp in seconds)
            auto now = std::chrono::system_clock::now();
            auto time_t = std::chrono::system_clock::to_time_t(now);
            std::stringstream filename;
            filename << "results/" << N << "_" << K << "_" << final_cost << "_" << time_t << ".txt";
            saved_filename = filename.str();
            
            // Save results to file
            std::ofstream outfile(saved_filename);
            if (outfile.is_open()) {
                outfile << "Round-Robin Scheduler Results\n";
                outfile << "============================\n";
                outfile << "Players: " << N << "\n";
                outfile << "Rounds: " << K << "\n";
                outfile << "Games per round: " << (N / 4) << "\n";
                outfile << "Attempt: " << attempt << "\n";
                outfile << "Initial cost: " << initial_cost << "\n";
                outfile << "Final cost: " << final_cost << "\n";
                
                scheduler.printStatisticsToStream(outfile);
                scheduler.printPairMatrixToStream(outfile);
                scheduler.printScheduleToStream(outfile);
                
                outfile.close();
                std::cout << "*** NEW BEST RESULT! ***\n";
                std::cout << "Results saved to: " << saved_filename << "\n";
            } else {
                std::cerr << "Error: Could not open file " << saved_filename << " for writing\n";
            }
        } else {
            std::cout << "Result not saved (cost " << final_cost << " >= best " << best_score << ")\n";
        }
        
        // Print to console
        scheduler.printStatistics();
        scheduler.printPairMatrix();
        
        // Check if ideal schedule found (cost == 0)
        if (final_cost == 0) {
            std::cout << "\n*** IDEAL SCHEDULE FOUND! No repeats! ***\n";
            if (!saved_filename.empty()) {
                std::cout << "Results saved to: " << saved_filename << "\n";
            }
            break;
        }
        
        std::cout << "\nContinuing to next attempt...\n";
    }
    
    return 0;
}

