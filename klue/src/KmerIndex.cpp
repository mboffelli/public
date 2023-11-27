#include <iostream>
#include <algorithm>
#include <random>
#include <string>
#include <sstream>
#include <iomanip>
#include <ctype.h>
#include <unordered_set>
#include <unordered_map>
#include <queue>
#include <functional>
#include "common.h"
#include "KmerIndex.h"
#include "ColoredCDBG.hpp"
#include <ExpressionParser.h>

// other helper functions
// pre: u is sorted
// checks if all elements in a sorted vector are unique
bool isUnique(const std::vector<int>& u) {} 

// removes duplicates from vector
std::vector<int> unique(const std::vector<int>& u) {} 

// returns the corresponding character ('A', 'C', 'G', or 'T') from a static array
const char Dna(int i) {} 

// calculates the Hamming distance (the number of positions at which the corresponding symbols are different) between a and b
int hamming(const char* a, const char* b) {} 

// generates a random temporary filename
std::string generate_tmp_file(std::string seed) {}

// Updated forward DFS (string-type)
// Perform DFS traversal starting from unitig tail
std::string extendUnitig(const UnitigMap<DataAccessor<void>, DataStorage<void>, false>& current, const int& current_color, std::string current_str, const std::unordered_set<int>& superset_colors, const int& k, std::unordered_set<std::string>& visited) {
    if (visited.find(current.getUnitigHead().rep().toString()) != visited.end()) { return "";  }
    visited.insert(current.getUnitigHead().rep().toString());

    std::string result = "";
    bool hasValidSuccessor = false;
    for (const auto& next : current.getSuccessors()) {
        UnitigColors::const_iterator it_next = next.getData()->getUnitigColors(next)->begin(next);
        UnitigColors::const_iterator it_next_end = next.getData()->getUnitigColors(next)->end();
        std::unordered_set<int> colors;
        for (; it_next != it_next_end; ++it_next) { colors.insert(it_next.getColorID()); }
        for (const auto& cc : colors) {
            if (std::find(superset_colors.begin(), superset_colors.end(), cc) != superset_colors.end() && cc == current_color && next.strand == current.strand) {
                hasValidSuccessor = true;
                std::string output;
                if (next.getUnitigKmer(next.dist).rep() == next.getUnitigKmer(next.dist).twin()) { output = next.getUnitigKmer(next.dist).rep().toString(); }
                else { output = next.getUnitigKmer(next.dist).twin().toString(); }
                if (next.getUnitigHead().rep().toString().substr(0, k - next.len) != current.getUnitigTail().toString().substr(current.getUnitigTail().toString().size() - k)) { output = next.getUnitigKmer(next.dist).toString(); }
                std::string extended_str = current_str + output.substr(output.length() - next.len);
                result += extendUnitig(next, current_color, extended_str, superset_colors, k, visited);
            }
        }
    }
    if (!hasValidSuccessor) { result += current_str + " "; } // at leaf node
    return result;
}

// Begin set operations
// Split user-inputted set operation commands by whitespace " "
std::vector<std::string> split(const std::string& s, char delimiter) {
    std::vector<std::string> tokens;
    std::string token;
    std::istringstream tokenStream(s);
    while (std::getline(tokenStream, token, delimiter)) {
        tokens.push_back(token);
    }
    return tokens;
}

// Recursively compute set operations (union and intersection are nodes) that correspond to specific sets (leaves)
std::set<int> computeSetOperation(const Node* root, const std::map<int, std::set<int>>& k_map) {
    if (!root) { return {}; }
    if (root->value != 'U' && root->value != 'I' && root->value != '\\') {
        auto it = k_map.find(root->value - 'A'); // 'A' maps to 0, 'B' to 1, etc.
        return (it != k_map.end()) ? it->second : std::set<int>{};
    }
    auto leftSet = computeSetOperation(root->left, k_map);
    auto rightSet = computeSetOperation(root->right, k_map);
    std::set<int> result;
    if (root->value == 'U') {
        std::set_union(leftSet.begin(), leftSet.end(), rightSet.begin(), rightSet.end(), std::inserter(result, result.end()));
    }
    else if (root->value == 'I') {
        std::set_intersection(leftSet.begin(), leftSet.end(), rightSet.begin(), rightSet.end(), std::inserter(result, result.end()));
    }
    else if (root->value == '\\') {
        std::set_difference(leftSet.begin(), leftSet.end(), rightSet.begin(), rightSet.end(), std::inserter(result, result.end()));
    }
    return result;
}
// End set operations

std::string stringOpt(size_t n, std::string option) {}

void KmerIndex::BuildReconstructionGraph(const ProgramOptions& opt) {}

void KmerIndex::BuildDistinguishingGraph(const ProgramOptions& opt, const std::vector<std::string>& transfasta, bool reconstruct) {
    // Initialize variables and set up temporary file handling based on reconstruction status
    // Involves setting k-mer length, creating temporary files, and initializing variables for sequence processing
    // initialization and setup ...
    
    // Process FASTA files based on the reconstruction flag and k-mer multiplicity
    // Involves reading sequences, converting non-standard nucleotides, trimming sequences, and writing to temporary files
    // read fasta and process sequences ...
    
    // Finalize sequence processing and clean up resources
    // Involves closing file streams, outputting warnings for non-standard characters, and memory management for file pointers
    // set CCDBG_Build options ...

    // Initialize thread pool and start processing k-mers from the graph
    // Involves setting up worker threads, distributing unitig processing tasks, and handling multithreading
    // initialize and start thread pool ...
    
    int range_discard = 0;
    int num_written = 0;
    std::string input_str = opt.input_set_operations;
    bool perform_set_operations = false;
    auto expressions = split(input_str, ' ');
    std::map<std::string, int> expr_to_int;
    if (!input_str.empty()) {
        perform_set_operations = true;
        ExpressionParser parser(input_str);  // create parser instance
        auto tokens = parser.tokenize(input_str);
        for (int i = 0; i < expressions.size(); ++i) { expr_to_int[expressions[i]] = i; }
        int maxWidth = 0;
        for (const auto& pair : expr_to_int) { maxWidth = std::max(maxWidth, static_cast<int>(pair.first.length())); }
        for (const auto& pair : expr_to_int) { std::cout << std::setw(maxWidth + 8) << pair.first << ": " << pair.second << "\n"; }
    }
    for (const auto& unitig : ccdbg) { // Iterate through all the unitigs in the de bruijn graph
        const UnitigColors* uc = unitig.getData()->getUnitigColors(unitig);
        const UnitigMap<DataAccessor<void>, DataStorage<void>, false> unitig_ = unitig;
        unitigs_v[n % unitigs_v.size()].push_back(std::make_pair(uc, unitig_)); // unitigs_v = vector of vectors of unitigs (b/c each thread contains a vector of unitigs)
        n++;
        if (unitigs_v[unitigs_v.size() - 1].size() >= thresh_size || n >= ccdbg.size()) { // If we're ready to start processing the unitigs using a series of workers
            for (size_t u_i = 0; u_i < unitigs_v.size(); u_i++) { // u_i = specific batch of unitigs that a worker will act on
                workers.emplace_back(
                    [&, u_i] {
                        std::ostringstream oss;
                        int _num_written = 0;
                        int _range_discard = 0;
                        for (auto unitig_x : unitigs_v[u_i]) { // Go through unitigs in the batch labeled u_i
                            auto uc = unitig_x.first;
                            auto& unitig = (unitig_x.second);
                            UnitigColors::const_iterator it_uc = uc->begin(unitig);
                            UnitigColors::const_iterator it_uc_end = uc->end();
                            std::map<int, std::set<int>> k_map; // key = color; value = list of positions (i.e. k-mers) along the current unitig (note: a k-mer is a position along a unitig)
                            std::unordered_set<int> superset_colors;
                            for (; it_uc != it_uc_end; ++it_uc) {
                                superset_colors.insert(it_uc.getColorID());
                                int color = color_map[it_uc.getColorID()];
                                k_map[color].insert(it_uc.getKmerPosition());
                                // DEBUG:
                                // std::cout << color << " " << unitig.getUnitigKmer(it_uc.getKmerPosition()).rep().toString() << " " << unitig.getUnitigKmer(it_uc.getKmerPosition()).toString() << " " << it_uc.getKmerPosition() << " " << unitig.strand << std::endl;
                            }
                            // begin extend 
                            if (opt.extend) {
                                std::unordered_set<std::string> visited;
                                k_map.clear();
                                for (const auto& c : superset_colors) {
                                    // extendUnitig(unitig, c, unitig.getUnitigTail().toString(), superset_colors, k, visited); // void
                                    std::string result = extendUnitig(unitig, c, unitig.getUnitigTail().toString(), superset_colors, k, visited);
                                    if (!result.empty()) {
                                        for (const auto& s : split(result, ' ')) {
                                            // [DEBUG] print extended contigs (all)
                                            // std::cout << ">" << c << "\n" << s << "\n";
                                            for (auto iter = s.begin(); iter != s.end(); ++iter) {
                                                int pos = std::distance(s.begin(), iter);
                                                k_map[c].insert(pos);
                                            }
                                        }
                                    }
                                }
                                // distinguish (default : find k-mer positions exlusive to each color)
                                std::string default_str = stringOpt(tmp_files.size(), "default"); // for n=8, generate "A\(AIBI...IH) B\(AIBI...IH) ... H\(AIBI...IH)"
                                auto all_expr = split(default_str, ' ');
                                for (const auto& expr : all_expr) {
                                    ExpressionParser expr_parser(expr);
                                    Node* root = expr_parser.parse();
                                    std::set<int> set_operation_result = computeSetOperation(root, k_map); // set of positions to keep
                                    // write out what remains among the contigs
                                    for (const auto& k_elem : k_map) {
                                        int curr_pos = -1;
                                        auto color = expr[0] - 'A'; // A=0, B=1, C=2, ..., H=7
                                        std::string colored_contig = "";
                                        for (const auto& pos : k_elem.second) {
                                            if (set_operation_result.count(pos)) {
                                                std::string km = unitig.getUnitigKmer(pos).toString();
                                                if (curr_pos == -1) { colored_contig = km; }
                                                else if (pos == curr_pos + 1) { colored_contig += km[km.length() - 1]; }
                                                else {
                                                    if (colored_contig.length() >= rb && colored_contig.length() <= re && k_elem.first == color) { oss << ">" << std::to_string(color) << "\n" << colored_contig << "\n"; _num_written++; }
                                                    else _range_discard++;
                                                    colored_contig = km;
                                                }
                                                curr_pos = pos;
                                            }
                                        }
                                        if (colored_contig != "") {
                                            if (colored_contig.length() >= rb && colored_contig.length() <= re && k_elem.first == color) { oss << ">" << std::to_string(color) << "\n" << colored_contig << "\n"; _num_written++; }
                                            else _range_discard++;
                                        }
                                    }
                                }
                            }
                            // END extend
                            std::set<int> positions_to_remove; // Positions (i.e. k-mers) along the current unitig that will be cut out
                            std::map<std::vector<int>, int> result_map; // Key = colors; Value = Position (i.e. k-mer)
                            std::stringstream ss; // For --combinations outputting aggregated colors
                            int int_to_print = 1;
                            if (perform_set_operations) {
                                for (const auto& expr : expressions) {
                                    ExpressionParser expr_parser(expr);
                                    Node* root = expr_parser.parse();
                                    std::set<int> set_operation_result = computeSetOperation(root, k_map); // set of positions to keep
                                    for (const auto& k_elem : k_map) {
                                        int curr_pos = -1;
                                        auto color = expr[0] - 'A'; // A=0, B=1, C=2, ..., H=7
                                        std::string colored_contig = "";
                                        for (const auto& pos : k_elem.second) {
                                            if (set_operation_result.count(pos)) {
                                                std::string km = unitig.getUnitigKmer(pos).toString();
                                                if (curr_pos == -1) { colored_contig = km; }
                                                else if (pos == curr_pos + 1) { colored_contig += km[km.length() - 1]; }
                                                else {
                                                    if (colored_contig.length() >= rb && colored_contig.length() <= re && k_elem.first == color) { oss << ">" << std::to_string(color) << "\n" << colored_contig << "\n"; _num_written++; }
                                                    else _range_discard++;
                                                    colored_contig = km;
                                                }
                                                curr_pos = pos;
                                            }
                                        }
                                        if (colored_contig != "") {
                                            if (colored_contig.length() >= rb && colored_contig.length() <= re && k_elem.first == color) { oss << ">" << std::to_string(color) << "\n" << colored_contig << "\n"; _num_written++; }
                                            else _range_discard++;
                                        }
                                    }
                                }
                                continue; // colored contigs already extracted, continue to next iteration
                            }
                            if (!opt.distinguish_union) { // If we don't specify --union (since if --union is specified, we don't actually need to do anything, remove any k-mers/positions, etc.)
                                if (!opt.distinguish_all_but_one_color && !opt.distinguish_combinations) { // Workflow: Find k-mers unique/exclusive to each color
                                    std::string default_str = stringOpt(tmp_files.size(), "default"); // for n=8, generate "A\(AIBI...IH) B\(AIBI...IH) ... H\(AIBI...IH)"
                                    auto all_expr = split(default_str, ' ');
                                    for (const auto& expr : all_expr) {
                                        ExpressionParser expr_parser(expr);
                                        Node* root = expr_parser.parse();
                                        std::set<int> set_operation_result = computeSetOperation(root, k_map); // set of positions to keep
                                        // write out what remains among the contigs
                                        for (const auto& k_elem : k_map) {
                                            int curr_pos = -1;
                                            auto color = expr[0] - 'A'; // A=0, B=1, C=2, ..., H=7
                                            std::string colored_contig = "";
                                            for (const auto& pos : k_elem.second) {
                                                if (set_operation_result.count(pos)) {
                                                    std::string km = unitig.getUnitigKmer(pos).toString();
                                                    if (curr_pos == -1) { colored_contig = km; }
                                                    else if (pos == curr_pos + 1) { colored_contig += km[km.length() - 1]; }
                                                    else {
                                                        if (colored_contig.length() >= rb && colored_contig.length() <= re && k_elem.first == color) { oss << ">" << std::to_string(color) << "\n" << colored_contig << "\n"; _num_written++; }
                                                        else _range_discard++;
                                                        colored_contig = km;
                                                    }
                                                    curr_pos = pos;
                                                }
                                            }
                                            if (colored_contig != "") {
                                                if (colored_contig.length() >= rb && colored_contig.length() <= re && k_elem.first == color) { oss << ">" << std::to_string(color) << "\n" << colored_contig << "\n"; _num_written++; }
                                                else _range_discard++;
                                            }
                                        }
                                    }          
                                    continue; // colored contigs already extracted, continue to next iteration
                                }
                                else if (!opt.distinguish_union && !opt.distinguish_combinations) { // workflow: opt.distinguish_all_but_one_color (e.g. if we have 8 colors, for each color, output all k-mers except those that are 8-colored)
                                    // method 1
                                    std::string all_but_one_str = stringOpt(tmp_files.size(), "all-but-one"); // generate "A\(AIBI...IH) B\(AIBI...IH) ... H\(AIBI...IH)"
                                    auto all_expr = split(all_but_one_str, ' ');
                                    for (const auto& expr : all_expr) {
                                        ExpressionParser expr_parser(expr);
                                        Node* root = expr_parser.parse();
                                        std::set<int> set_operation_result = computeSetOperation(root, k_map); // set of positions to keep
                                        // write out what remains among the contigs
                                        for (const auto& k_elem : k_map) {
                                            int curr_pos = -1;
                                            auto color = expr[0] - 'A'; // A=0, B=1, C=2, ..., H=7
                                            std::string colored_contig = "";
                                            for (const auto& pos : k_elem.second) {
                                                if (set_operation_result.count(pos)) {
                                                    std::string km = unitig.getUnitigKmer(pos).toString();
                                                    if (curr_pos == -1) { colored_contig = km; }
                                                    else if (pos == curr_pos + 1) { colored_contig += km[km.length() - 1]; }
                                                    else {
                                                        if (colored_contig.length() >= rb && colored_contig.length() <= re && k_elem.first == color) { oss << ">" << std::to_string(color) << "\n" << colored_contig << "\n"; _num_written++; }
                                                        else _range_discard++;
                                                        colored_contig = km;
                                                    }
                                                    curr_pos = pos;
                                                }
                                            }
                                            if (colored_contig != "") {
                                                if (colored_contig.length() >= rb && colored_contig.length() <= re && k_elem.first == color) { oss << ">" << std::to_string(color) << "\n" << colored_contig << "\n"; _num_written++; }
                                                else _range_discard++;
                                            }
                                        }
                                    }
                                    continue; // colored contigs already extracted, continue to next iteration
                                }
                                else if (!opt.distinguish_all_but_one_color) { // workflow: opt.distinguish_combinations (output every combination [i.e. each area in the Venn diagram], with colors separated by underscores)
                                    std::set<int> exclusive_positions;
                                    for (const auto& k_elem : k_map) { // Iterate over each color and get exclusive positions for that color
                                        bool is_exclusive = true;
                                        std::set<int> current_positions = k_elem.second;
                                        for (const auto& other_elem : k_map) { // Check if the positions of this color appear in other colors
                                            if (k_elem.first != other_elem.first) {
                                                std::set<int> intersection;
                                                std::set_intersection(current_positions.begin(), current_positions.end(),
                                                    other_elem.second.begin(), other_elem.second.end(),
                                                    std::inserter(intersection, intersection.begin()));
                                                if (!intersection.empty()) {
                                                    is_exclusive = false;
                                                    break;
                                                }
                                            }
                                        }
                                        if (is_exclusive) { exclusive_positions.insert(current_positions.begin(), current_positions.end()); }
                                    }
                                    positions_to_remove.insert(exclusive_positions.begin(), exclusive_positions.end()); // Add the exclusive positions to the positions_to_remove
                                    if (k_map.size() >= 2) {
                                        std::vector<int> consolidated_key;
                                        for (auto it = k_map.begin(); it != k_map.end(); ++it) { // Iterate over k_map and build the consolidated key
                                            consolidated_key.push_back(it->first);
                                        }
                                        int integer_value = *(k_map.begin()->second.begin());
                                        result_map[consolidated_key] = integer_value; // {0 1 2} : 0 (color : position)                                        
                                    }
                                    for (const auto& k_elem : result_map) {
                                        int curr_pos = -1;
                                        std::string colored_contig = "";
                                        auto color = k_elem.first;
                                        auto pos = k_elem.second;
                                        for (int i = 0; i < color.size(); ++i) {
                                            ss << color[i];
                                            if (i < color.size() - 1) {
                                                ss << "_";
                                            }
                                        }
                                    }
                                }
                            } // end if !opt.distinguish_union
                            
                            // Now, write out what remains among the contigs
                            for (const auto& k_elem : k_map) {
                                int curr_pos = -1;
                                std::string colored_contig = "";
                                auto color = k_elem.first;
                                std::string color_key = std::to_string(color);
                                //std::string contig_metadata = " :" + unitig.dist + "," + unitig.len + "," + unitig.size + "," + unitig.strand;
                                if (opt.distinguish_combinations) { color_key = ss.str(); int_to_print++; }
                                if (int_to_print == k_map.size()) {
                                    for (const auto& pos : k_elem.second) {
                                        if (!positions_to_remove.count(pos)) {
                                            std::string km = unitig.getUnitigKmer(pos).toString();
                                            if (curr_pos == -1) { // How to correspond color?
                                                colored_contig = km;
                                            }
                                            else if (pos == curr_pos + 1) {
                                                colored_contig += km[km.length() - 1];
                                            }
                                            else {
                                                if (colored_contig.length() >= rb && colored_contig.length() <= re) { oss << ">" << color_key << "\n" << colored_contig << "\n"; _num_written++; }
                                                else _range_discard++;
                                                colored_contig = km;
                                            }
                                            curr_pos = pos;
                                        }
                                    }
                                    if (colored_contig != "") {
                                        if (colored_contig.length() >= rb && colored_contig.length() <= re) { oss << ">" << color_key << "\n" /*<< sequence_to_prepend*/ << colored_contig /*<< sequence_to_append*/ << "\n"; _num_written++; }
                                        else _range_discard++;
                                    }
                                }
                            }
                        }
                        std::unique_lock<std::mutex> lock(mutex_unitigs);
                        o << oss.str();
                        num_written += _num_written;
                        range_discard += _range_discard;
                    }
                );
            }
            for (auto& t : workers) t.join();
            workers.clear();
            for (size_t u_i = 0; u_i < unitigs_v.size(); u_i++) {
                unitigs_v[u_i].clear();
            }
        }
    }
    o.flush();
    ccdbg.clear(); // Free memory associated with the colored compact dBG
    ncolors = tmp_files.size(); // Record the number of "colors"
    for (auto tmp_file : tmp_files) std::remove(tmp_file.c_str()); // Remove temp files needed to make colored graph
    if (!opt.stream_out) {
        of.close();
    }
    tmp_files.clear();
    if (opt.verbose) {
        if (range_discard > 0) std::cerr << "[build] Number of output sequences filtered out due to length: " << range_discard << std::endl;
        std::cerr << "[build] Number of output sequences written: " << num_written << std::endl;
    }
}
