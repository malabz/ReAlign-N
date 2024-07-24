#pragma once

#include <chrono>
#include <iostream>

#include "Fasta.h"

// print help message
void displayHelp() {
    std::cout << std::endl;
    std::cout << "Usage: /.realign_n [-r] path [-a] path [-o] path [-m] mode" << std::endl;
    std::cout << std::endl;
    std::cout << "  Necessary arguments:" << std::endl;
    std::cout << "    -r  Specify the path of raw data, a file in FASTA format." << std::endl;
    std::cout << "    -a  Specify the path of initial alignment, a file in FASTA format." << std::endl;
    std::cout << std::endl;
    std::cout << "  Optional arguments:" << std::endl; 
    std::cout << "    -o  Specify the output for ReAlign-N, a file in FASTA format." << std::endl;
    std::cout << "    -m  Specify the minium split distance of match (default based on the similarity)." << std::endl;
    std::cout << "    -e  Specify the minium split distance of entropy (default based on the similarity)." << std::endl;
    std::cout << "    -p  Specify the pattern of ReAlign-N (default pattern: 1)." << std::endl;
    std::cout << "          1 for local realignment followed by global realignment." << std::endl;
    std::cout << "          2 for global realignment followed by local realignment." << std::endl;
    std::cout << "    -h  Print the help message." << std::endl;
}

// print current local time
void printCurrentLocalTime() {
    // 获取当前系统时钟的时间点
    auto now = std::chrono::system_clock::now();

    // 将时间点转换为时间类型
    std::time_t now_time = std::chrono::system_clock::to_time_t(now);

    // 将时间类型转换为本地时间
    std::tm* local_time = std::localtime(&now_time);

    // 输出本地时间信息
    char buffer[80];
    std::strftime(buffer, 80, "[%Y-%m-%d %H:%M:%S] ", local_time);
    std::cout << buffer;
}

static utils::Fasta read_from(char *const file_path) {
    std::ifstream file(file_path);
    if (!file) {
        std::cerr << "Error: cannot open file " << file_path << std::endl;
        std::cerr << "Please check that the file path is correct and make sure the file exists." << std::endl;
        exit(1);
    }
    utils::Fasta fasta(file);
    file.close();
    return fasta;
}

static utils::Fasta read_from(std::string file_path) {
    std::ifstream file(file_path);
    if (!file) {
        std::cerr << "Error: cannot open file " << file_path << std::endl;
        std::cerr << "Please check that the file path is correct and make sure the file exists." << std::endl;
        exit(1);
    }
    utils::Fasta fasta(file);
    file.close();
    return fasta;
}

std::string toLower(const std::string& str) {
    std::string lowerStr = str;
    std::transform(lowerStr.begin(), lowerStr.end(), lowerStr.begin(), ::tolower);
    return lowerStr;
}

static double compare_seq_similarity(const std::string& s1, const std::string& s2) {
    if (s1.size() != s2.size()) {
        throw std::invalid_argument("Strings must be of the same length");
    }

    std::string s1Lower = toLower(s1);
    std::string s2Lower = toLower(s2);

    int n = s1Lower.size();
    int cnt = 0;

    for (int i = 0; i < n; ++i) {
        if (s1Lower[i] == s2Lower[i]) {
            cnt++;
        }
    }

    return static_cast<double>(cnt) / n;
}

static double compute_dataset_similarity(const std::string initial_alignment) {
    utils::Fasta alignment = read_from(initial_alignment);
    
    double total_similarity = 0.0;
    int total_pairs = 0;

    for (size_t i = 0; i < alignment.sequences.size(); ++i) {
        for (size_t j = i + 1; j < alignment.sequences.size(); ++j) {
            double curr_similarity = compare_seq_similarity(alignment.sequences[i], alignment.sequences[j]);
            total_similarity += curr_similarity;
            total_pairs += 1;
        }
    }
    
    double average_similarity = total_similarity / total_pairs;
    // 保留两位小数并四舍五入
    return std::round(average_similarity * 100) / 100.0; 
}

static unsigned int compute_min_block_size(const double similarity, const int sequences_size, const std::string mode) {
    if (sequences_size >= 200) {
        if (similarity >= 0.9) {
            return 10;
        } else {
            return 50;
        }
    } else {
        if (mode == "match") {
            if (similarity >= 0.99) {
                return 5;
            } else if (similarity >= 0.96) {
                return 10;
            } else if (similarity >= 0.93) {
                return 15;
            } else if (similarity >= 0.9) {
                return 20;
            } else if (similarity >= 0.8) {
                return 30;
            } else {
                return 50;
            }
        } else if (mode == "entropy"){
            if (similarity >= 0.97) {
                return 5;
            } else if (similarity >= 0.9) {
                return 50;
            } else {
                return 5;
            }
        }
    }
}