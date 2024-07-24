#include <iostream>
#include <vector>
#include <numeric>
#include <algorithm>
#include <fstream>
#include <cstring>
#include <stdlib.h>
#include <cmath>

#include "Fasta.h"
#include "Utils.h"

static constexpr long long     MISMATCH = -1;
static constexpr long long        MATCH =  1;
static constexpr long long GAPEXTENSION = -2;
static constexpr long long      GAPOPEN =  0;

static std::vector<std::string> realigned_sequences;
utils::Fasta tmp_block;

static long long score(std::vector<std::string> sequences, unsigned j) {
    unsigned counts[128];
    memset(counts, 0, sizeof(counts));
    for (unsigned i = 0; i != sequences.size(); ++i) {
        char curr_base = sequences[i][j];
        if (curr_base == 'u' || curr_base == 'U') {
            curr_base = 't';
        }
        ++counts[curr_base];
    }

    unsigned const a = counts['a'] + counts['A'];
    unsigned const g = counts['g'] + counts['G'];
    unsigned const c = counts['c'] + counts['C'];
    unsigned const t = counts['t'] + counts['T'];
    unsigned const gap = counts['-'];
    unsigned const n = std::accumulate(counts, counts + 128, 0u) - a - g - c - t - gap;
    return ((a + c) * (g + t) + a * c + g * t) * MISMATCH + (a * (a - 1) + c * (c - 1) + g * (g - 1) + t * (t - 1)) / 2 * MATCH + (a + c + g + t + n) * gap * GAPEXTENSION;
}

static long long score(std::vector<std::string> sequences, unsigned l, unsigned r) {
    long long s = 0;
    for (unsigned i = l; i != r; ++i)
        s += score(sequences, i);
    return s;
}


bool executeCommand(const std::string& command) {
    FILE* pipe = popen(command.c_str(), "r");
    if (!pipe) {
        std::cerr << "Error: Unable to execute command." << std::endl;
        return false;
    }

    char buffer[128];
    std::string result = "";
    while (!feof(pipe)) {
        if (fgets(buffer, 128, pipe) != NULL)
            result += buffer;
    }
    pclose(pipe);

    // 检查命令是否执行成功
    if (!result.empty()) {
        std::cerr << "Error message: " << result << std::endl;
        return false;
    }

    return true;
}

static std::vector<unsigned int> find_all_matched_columns(utils::Fasta const &fasta) {
    std::vector<unsigned int> all_matched_columns;
    for (int j = 0; j < fasta.sequences[0].size(); j++) {
        char curr_base = fasta.sequences[0][j];
        bool all_matched = true;
        if (curr_base != '-') {
            for (int i = 0; i < fasta.sequences.size(); i++) {
                if (curr_base != fasta.sequences[i][j]) {
                    all_matched = false;
                    break;
                }
            }
            if (all_matched) {
                all_matched_columns.push_back(j);
            }
        }
    }
    return all_matched_columns;
}

static float a_log4_a(int a, float b) {
    return (a/b) * (log(a/b)/log(4));
}

static std::vector<unsigned int> find_clear_columns(utils::Fasta const &fasta) {
    std::vector<unsigned int> clear_columns;
    float row = fasta.sequences.size();
    unsigned counts[128];
    for (unsigned j = 0; j != fasta.sequences[0].size(); ++j) {
        memset(counts, 0, sizeof(counts));
        for (unsigned i = 0; i != fasta.sequences.size(); ++i) {
            char curr_base = fasta.sequences[i][j];
            if (curr_base == 'u' || curr_base == 'U') {
                curr_base = 't';
            }
            ++counts[curr_base];
        }
        unsigned const a = counts['a'] + counts['A'];
        unsigned const g = counts['g'] + counts['G'];
        unsigned const c = counts['c'] + counts['C'];
        unsigned const t = counts['t'] + counts['T'];
        unsigned const gap = counts['-'];
        if (gap / row <= 0.3) {
            float entropy = a_log4_a(a, row-gap) + a_log4_a(g, row-gap) + a_log4_a(c, row-gap) + a_log4_a(t, row-gap);
            if (entropy > -0.5) {
                clear_columns.push_back(j);
            }
        }
    }
    return clear_columns;
}

std::vector<unsigned int> find_cutting_sites_of_low_quality_blocks(utils::Fasta const &fasta,
                                                                   std::vector<unsigned int> matched_columns,
                                                                   int min_block_length) {
    std::vector<unsigned int> cutting_sites;
    for (int i = 1; i < matched_columns.size(); i++) {
        if ((matched_columns[i] - matched_columns[i-1]) >= min_block_length) {
            cutting_sites.push_back(matched_columns[i-1]);
            cutting_sites.push_back(matched_columns[i]);
        }
    }
    return cutting_sites;
}

std::vector<std::string> cut_blocks(utils::Fasta const &fasta, int start_site, int end_site) {
    std::vector<std::string> blocks;
    for (int i = 0; i < fasta.sequences.size(); i++) {
        blocks.push_back(fasta.sequences[i].substr(start_site + 1, (end_site - start_site - 1)));
    }
    return blocks;
}

std::vector<std::string> split_sequences(utils::Fasta const &fasta, int start_site, int end_site) {
    std::vector<std::string> blocks;
    for (int i = 0; i < fasta.sequences.size(); i++) {
        blocks.push_back(fasta.sequences[i].substr(start_site, (end_site - start_site + 1)));
    }
    return blocks;
}

void join_blocks(std::vector<std::string> curr_block) {
    if (realigned_sequences.size() == 0) {
        realigned_sequences.resize(curr_block.size());
    }
    for (int i = 0; i < curr_block.size(); i++) {
        realigned_sequences[i] += curr_block[i];
    }
}

std::vector<std::string> realign_block(utils::Fasta const &fasta, int start_site, int end_site) {
    std::vector<std::string> low_quality_block = cut_blocks(fasta, start_site, end_site);
    double before_realign_score = score(low_quality_block, 0, low_quality_block[0].size() - 1);

    std::ofstream ofs("tmp/tmp_block.fasta");
    if (!ofs) {
        std::cerr << "Error: cannot open file tmp/tmp_block.fasta" << std::endl;
        std::cerr << "Please make sure the file path is correct and has appropriate permissions." << std::endl;
        exit(1);
    }
    tmp_block.identifications = fasta.identifications;
    tmp_block.sequences = low_quality_block;
    tmp_block.write_to(ofs);
    ofs.close();

    // system("mafft --retree 1 tmp/tmp_block.fasta > tmp/aligned_block.fasta 2>&1");
    std::string command = "mafft --retree 1 tmp/tmp_block.fasta > tmp/aligned_block.fasta 2>&1";
    if (!executeCommand(command)) {
        std::cerr << "Please make sure mafft can be used normally." << std::endl; 
        exit(1);
    }
    
    utils::Fasta after_realigned_block = read_from("tmp/aligned_block.fasta");
    double after_realign_score = score(after_realigned_block.sequences, 0, after_realigned_block.sequences[0].size() - 1);
    if (after_realign_score > before_realign_score) {
        return after_realigned_block.sequences;
    } else {
        return low_quality_block;
    }
}

void realign_alignment(utils::Fasta const &fasta, std::vector<unsigned int> cutting_sites, int min_length) {

    realigned_sequences.resize(fasta.sequences.size());
    std::vector<std::string> realigned_block, high_quality_block;
    realigned_block.resize(fasta.sequences.size());
    high_quality_block.resize(fasta.sequences.size());

    int i = 0;
    if (cutting_sites.size() > 2) {
        while (i < cutting_sites.size() - 1) {
            if (i == 0) {
                if (cutting_sites[0] > min_length) {
                    realigned_sequences = realign_block(fasta, -1, cutting_sites[0] + 1);
                } else {
                    realigned_sequences = split_sequences(fasta, 0, cutting_sites[0]);
                }
                realigned_block = realign_block(fasta, cutting_sites[i], cutting_sites[i+1]);
                join_blocks(realigned_block);
            } else {
                high_quality_block = split_sequences(fasta, cutting_sites[i-1], cutting_sites[i]);
                join_blocks(high_quality_block);
                realigned_block = realign_block(fasta, cutting_sites[i], cutting_sites[i+1]);
                join_blocks(realigned_block);
                if (i == cutting_sites.size() - 2) {
                    if (fasta.sequences[0].size() - cutting_sites.back() - 1 > min_length) {
                        realigned_block = realign_block(fasta, cutting_sites.back() - 1, fasta.sequences[0].size());
                        join_blocks(realigned_block);
                    } else {
                        high_quality_block = split_sequences(fasta, cutting_sites.back(), fasta.sequences.size());
                        join_blocks(high_quality_block);
                    }
                }
            }
            i += 2;
        }
    } else if (cutting_sites.size() != 0) {
        if (cutting_sites[0] > min_length) {
            realigned_sequences = realign_block(fasta, -1, cutting_sites[0]+1);
        } else {
            realigned_sequences = split_sequences(fasta, 0, cutting_sites[0]);
        }
        realigned_block = realign_block(fasta, cutting_sites[i], cutting_sites[i+1]);
        join_blocks(realigned_block);
        if (fasta.sequences[0].size() - cutting_sites.back() - 1 > min_length) {
            realigned_block = realign_block(fasta, cutting_sites.back() - 1, fasta.sequences[0].size());
            join_blocks(realigned_block);
        } else {
            high_quality_block = split_sequences(fasta, cutting_sites.back(), fasta.sequences.size());
            join_blocks(high_quality_block);
        }
    } else {
        for (int i = 0; i < fasta.sequences.size(); i++) {
            realigned_sequences[i] = fasta.sequences[i];
        }
        // printf("No regions need optimization.\n");
    }

}

void local_realignment(const std::string initial_alignment, const int min_match_block_size, 
                       const int min_entropy_block_size, const std::string local_realigned_output){

    utils::Fasta alignment = read_from(initial_alignment);
    std::vector<std::string> matched_results;
    std::vector<std::string> entropy_results;
        
    printCurrentLocalTime();
    std::cout << "Matched local realigning..." << std::endl;
    // Find all matched columns in the initial alignment
    std::vector<unsigned int> matched_columns = find_all_matched_columns(alignment);

    // Filter out the start and end sites of low-quality blocks
    std::vector<unsigned int> blocks_cutting_sites = find_cutting_sites_of_low_quality_blocks(alignment, matched_columns,min_match_block_size);

    // Start to realign the initial alignment
    realign_alignment(alignment, blocks_cutting_sites, min_match_block_size);
    matched_results.swap(realigned_sequences);
    realigned_sequences.clear();
    printCurrentLocalTime();
    std::cout << "Done." << std::endl;

    printCurrentLocalTime();
    std::cout << "Entropy local realigning..." << std::endl;
    // Find all clear columns in the initial alignment
    std::vector<unsigned int> entropy_columns = find_clear_columns(alignment);

    // Filter out the start and end sites of low-quality blocks
    std::vector<unsigned int> blocks_cutting_sites_entropy = find_cutting_sites_of_low_quality_blocks(alignment, entropy_columns,min_entropy_block_size);

    // Start to realign the initial alignment
    realign_alignment(alignment, blocks_cutting_sites_entropy, min_entropy_block_size);
    entropy_results.swap(realigned_sequences);
    realigned_sequences.clear();
    printCurrentLocalTime();
    std::cout << "Done." << std::endl;

    // Save the realigned sequences
    std::ofstream ofs(local_realigned_output);
    if (!ofs) {
        std::cerr << "Error: cannot open file " << local_realigned_output << "." << std::endl;
        std::cerr << "Please make sure the file path is correct and has appropriate permissions." << std::endl;
        // return 1;
    }
    unsigned int matched_score = score(matched_results, 0, matched_results[0].size() - 1);
    unsigned int entropy_score = score(entropy_results, 0, entropy_results[0].size() - 1);
        
    printCurrentLocalTime();
    if(matched_score > entropy_score) {
        alignment.sequences = std::move(matched_results);
        alignment.write_to(ofs);
        ofs.close();
        std::cout << "Matched result is better.\n";
    } else {
        alignment.sequences = std::move(entropy_results);
        alignment.write_to(ofs);
        ofs.close();
        std::cout << "Entropy result is better.\n";
    }
    matched_results.clear();
    entropy_results.clear();
    printCurrentLocalTime();
    std::cout << "Done." << std::endl; 
}



void global_realignment(std::string unaligned_seq, std::string output, std::string initial_alignment) {
    std::string command = "global_realignment -i ";
    command += unaligned_seq;
    command += " -o ";
    command += output;
    command += " -a ";
    command += initial_alignment;
    command += " -t 1";
    // Print the command
    // system(command.c_str());
    if (!executeCommand(command)) {
        std::cerr << "Global realignment met an error." << std::endl; 
        exit(1);
    }
}

int main(int argc, char **argv) {
    // Check if the program was called with no arguments or with the "-h" option
    if (argc == 1 || (argc == 2 && std::string(argv[1]) == "-h")) {
        displayHelp();
        return 0;
    }

    std::string unaligned_seq, initial_alignment, min_match_block_size, min_entropy_block_size;
    std::string output = "realign_n_result.fas";
    std::string pattern = "1";
    bool has_min_match_block_size = false;
    bool has_min_entropy_block_size = false;


    for (int i = 1; i < argc; i += 2) {
        std::string option = argv[i];
        if (i + 1 < argc) {
            std::string value = argv[i + 1];
            if (option == "-r") {
                unaligned_seq = value;
            } else if (option == "-a") {
                initial_alignment = value;
            } else if (option == "-o") {
                output = value;
            } else if (option == "-p") {
                pattern = value;
            } else if (option == "-m") {
                min_match_block_size = value;
                has_min_match_block_size = true;
            } else if (option == "-e") {
                min_entropy_block_size = value;
                has_min_entropy_block_size = true;
            } else {
                std::cerr << "** Unknown option: " << option << std::endl;
                displayHelp();
                return 1;
            }
        } else {
            std::cerr << "** Missing value for option: " << option << std::endl;
            displayHelp();
            return 1;
        }
    }

    // Check if required -i and -a options are provided
    if (unaligned_seq.empty() || initial_alignment.empty()) {
        std::cerr << "** Error: Both -r and -a options are required." << std::endl;
        displayHelp();
        return 1;
    }

    // Check if required -i and -a options are provided
    if (!(pattern == "1" || pattern == "2")) {
        std::cerr << "** Error: -p can only be set to 1 or 2." << std::endl;
        displayHelp();
        return 1;
    }

    // Print the arguments
    std::cout << std::endl << "**" << std::endl;
    std::cout << "** Unaligned sequences: " << unaligned_seq << std::endl;
    std::cout << "** Initial alignment  : " << initial_alignment << std::endl;
    std::cout << "** Output file        : " << output << std::endl;
    std::cout << "** Realignment pattern: " << pattern << std::endl;
    std::cout << "**" << std::endl << std::endl;

    system("mkdir tmp");

    utils::Fasta raw_data = read_from(unaligned_seq);

    if (!has_min_entropy_block_size && !has_min_match_block_size) {
        if (pattern == "1") {
            // Local realignment begin
            double avg_similarity = compute_dataset_similarity(initial_alignment);
            int min_match_block_size_INT = compute_min_block_size(avg_similarity, raw_data.sequences.size(), "match");
            int min_entropy_block_size_INT = compute_min_block_size(avg_similarity, raw_data.sequences.size(), "entropy");
            printCurrentLocalTime();
            std::cout << "Average similarity: " << avg_similarity << std::endl;
            printCurrentLocalTime();
            std::cout << "Minimum match-based solit block size: " << min_match_block_size_INT << std::endl; 
            printCurrentLocalTime();
            std::cout << "Minimum entropy-based splt block size: " << min_entropy_block_size_INT << std::endl;  
            local_realignment(initial_alignment, min_match_block_size_INT, min_entropy_block_size_INT, "tmp/local_realigned.fas");
            // Local realignment end

            // Global realignment begin 
            printCurrentLocalTime();
            std::cout << "Global realigning..." << std::endl;
            global_realignment(unaligned_seq, output, "tmp/local_realigned.fas");
            printCurrentLocalTime();
            std::cout << "Done." << std::endl;
            // Global realignment end
        }
        if (pattern == "2") {
            // Global realignment begin
            printCurrentLocalTime();
            std::cout << "Global realigning..." << std::endl;
            global_realignment(unaligned_seq, "tmp/global_realigned.fas", initial_alignment);
            printCurrentLocalTime();
            std::cout << "Done." << std::endl;
            // Global realignment end

            // Local realignment begin
            double avg_similarity = compute_dataset_similarity("tmp/global_realigned.fas");
            int min_match_block_size_INT = compute_min_block_size(avg_similarity, raw_data.sequences.size(), "match");
            int min_entropy_block_size_INT = compute_min_block_size(avg_similarity, raw_data.sequences.size(), "entropy");
            printCurrentLocalTime();
            std::cout << "Average similarity: " << avg_similarity << std::endl;
            printCurrentLocalTime();
            std::cout << "Minimum match-based solit block size: " << min_match_block_size_INT << std::endl; 
            printCurrentLocalTime();
            std::cout << "Minimum entropy-based splt block size: " << min_entropy_block_size_INT << std::endl;  
            local_realignment("tmp/global_realigned.fas", min_match_block_size_INT, min_entropy_block_size_INT, output);
            // Local realignment end
        }
    } else if (has_min_entropy_block_size && !has_min_match_block_size) {
        if (pattern == "1") {
            // Local realignment begin
            double avg_similarity = compute_dataset_similarity(initial_alignment);
            int min_match_block_size_INT = compute_min_block_size(avg_similarity, raw_data.sequences.size(), "match");
            int min_entropy_block_size_INT = std::stoi(min_entropy_block_size);
            printCurrentLocalTime();
            std::cout << "Average similarity: " << avg_similarity << std::endl;
            printCurrentLocalTime();
            std::cout << "Minimum match-based solit block size: " << min_match_block_size_INT << std::endl; 
            printCurrentLocalTime();
            std::cout << "Minimum entropy-based splt block size: " << min_entropy_block_size_INT << std::endl;  
            local_realignment(initial_alignment, min_match_block_size_INT, min_entropy_block_size_INT, "tmp/local_realigned.fas");
            // Local realignment end

            // Global realignment begin 
            printCurrentLocalTime();
            std::cout << "Global realigning..." << std::endl;
            global_realignment(unaligned_seq, output, "tmp/local_realigned.fas");
            printCurrentLocalTime();
            std::cout << "Done." << std::endl;
            // Global realignment end
        }
        if (pattern == "2") {
            // Global realignment begin
            printCurrentLocalTime();
            std::cout << "Global realigning..." << std::endl;
            global_realignment(unaligned_seq, "tmp/global_realigned.fas", initial_alignment);
            printCurrentLocalTime();
            std::cout << "Done." << std::endl;
            // Global realignment end

            // Local realignment begin
            double avg_similarity = compute_dataset_similarity("tmp/global_realigned.fas");
            int min_match_block_size_INT = compute_min_block_size(avg_similarity, raw_data.sequences.size(), "match");
            int min_entropy_block_size_INT = std::stoi(min_entropy_block_size);
            printCurrentLocalTime();
            std::cout << "Average similarity: " << avg_similarity << std::endl;
            printCurrentLocalTime();
            std::cout << "Minimum match-based solit block size: " << min_match_block_size_INT << std::endl; 
            printCurrentLocalTime();
            std::cout << "Minimum entropy-based splt block size: " << min_entropy_block_size_INT << std::endl;  
            local_realignment("tmp/global_realigned.fas", min_match_block_size_INT, min_entropy_block_size_INT, output);
            // Local realignment end
        }
    } else if (!has_min_entropy_block_size && has_min_match_block_size) {
        if (pattern == "1") {
            // Local realignment begin
            double avg_similarity = compute_dataset_similarity(initial_alignment);
            int min_match_block_size_INT = std::stoi(min_match_block_size);
            int min_entropy_block_size_INT = compute_min_block_size(avg_similarity, raw_data.sequences.size(), "entropy");
            printCurrentLocalTime();
            std::cout << "Average similarity: " << avg_similarity << std::endl;
            printCurrentLocalTime();
            std::cout << "Minimum match-based solit block size: " << min_match_block_size_INT << std::endl; 
            printCurrentLocalTime();
            std::cout << "Minimum entropy-based splt block size: " << min_entropy_block_size_INT << std::endl;  
            local_realignment(initial_alignment, min_match_block_size_INT, min_entropy_block_size_INT, "tmp/local_realigned.fas");
            // Local realignment end

            // Global realignment begin 
            printCurrentLocalTime();
            std::cout << "Global realigning..." << std::endl;
            global_realignment(unaligned_seq, output, "tmp/local_realigned.fas");
            printCurrentLocalTime();
            std::cout << "Done." << std::endl;
            // Global realignment end
        }
        if (pattern == "2") {
            // Global realignment begin
            printCurrentLocalTime();
            std::cout << "Global realigning..." << std::endl;
            global_realignment(unaligned_seq, "tmp/global_realigned.fas", initial_alignment);
            printCurrentLocalTime();
            std::cout << "Done." << std::endl;
            // Global realignment end

            // Local realignment begin
            double avg_similarity = compute_dataset_similarity("tmp/global_realigned.fas");
            int min_match_block_size_INT = std::stoi(min_match_block_size);
            int min_entropy_block_size_INT = compute_min_block_size(avg_similarity, raw_data.sequences.size(), "entropy");
            printCurrentLocalTime();
            std::cout << "Average similarity: " << avg_similarity << std::endl;
            printCurrentLocalTime();
            std::cout << "Minimum match-based solit block size: " << min_match_block_size_INT << std::endl; 
            printCurrentLocalTime();
            std::cout << "Minimum entropy-based splt block size: " << min_entropy_block_size_INT << std::endl;  
            local_realignment("tmp/global_realigned.fas", min_match_block_size_INT, min_entropy_block_size_INT, output);
            // Local realignment end
        }
    } else if (has_min_entropy_block_size && has_min_match_block_size) {
        if (pattern == "1") {
            // Local realignment begin
            int min_match_block_size_INT = std::stoi(min_match_block_size);
            int min_entropy_block_size_INT = std::stoi(min_entropy_block_size);
            printCurrentLocalTime();
            std::cout << "Minimum match-based solit block size: " << min_match_block_size_INT << std::endl; 
            printCurrentLocalTime();
            std::cout << "Minimum entropy-based splt block size: " << min_entropy_block_size_INT << std::endl;  
            local_realignment(initial_alignment, min_match_block_size_INT, min_entropy_block_size_INT, "tmp/local_realigned.fas");
            // Local realignment end

            // Global realignment begin 
            printCurrentLocalTime();
            std::cout << "Global realigning..." << std::endl;
            global_realignment(unaligned_seq, output, "tmp/local_realigned.fas");
            printCurrentLocalTime();
            std::cout << "Done." << std::endl;
            // Global realignment end
        }
        if (pattern == "2") {
            // Global realignment begin
            printCurrentLocalTime();
            std::cout << "Global realigning..." << std::endl;
            global_realignment(unaligned_seq, "tmp/global_realigned.fas", initial_alignment);
            printCurrentLocalTime();
            std::cout << "Done." << std::endl;
            // Global realignment end

            // Local realignment begin
            int min_match_block_size_INT = std::stoi(min_match_block_size);
            int min_entropy_block_size_INT = std::stoi(min_entropy_block_size);
            printCurrentLocalTime();
            std::cout << "Minimum match-based solit block size: " << min_match_block_size_INT << std::endl; 
            printCurrentLocalTime();
            std::cout << "Minimum entropy-based splt block size: " << min_entropy_block_size_INT << std::endl;  
            local_realignment("tmp/global_realigned.fas", min_match_block_size_INT, min_entropy_block_size_INT, output);
            // Local realignment end
        }
    }
    
    system("rm -r tmp");
    return 0;
}
