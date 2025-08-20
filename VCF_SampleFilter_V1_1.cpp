//BT June 12, 2025

//Four-Thread Architecture:

//Reader thread: Streams lines from input file
//Worker threads: Process lines (configurable number)
//Writer thread: Streams output to file
//Main thread: Coordinates everything


//Memory Optimizations:

//Uses std::move() to avoid unnecessary string copies
//Reserves vector space to prevent reallocations
//Processes and immediately releases each line


//Progress Monitoring: Shows progress every 10,000 lines so you can see it's working

//g++ -std=c++11 -O3 -o vcf_filter vcf_filter.cpp -lz -lpthread

//# Use fewer threads initially to test
//./vcf_filter -i input.vcf.gz -o output.vcf -s samples.txt -t 2


#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <unordered_set>
#include <thread>
#include <mutex>
#include <queue>
#include <condition_variable>
#include <algorithm>
#include <atomic>
#include <zlib.h>

class VCFSampleFilter {
private:
    std::unordered_set<std::string> target_samples;
    std::vector<int> sample_indices;
    std::string input_file;
    std::string output_file;
    std::string sample_file;
    bool compress_output;
    int num_threads;
    
    // Thread-safe queues with size limits to prevent memory overflow
    static const size_t MAX_QUEUE_SIZE = 1000;
    std::queue<std::string> input_queue;
    std::queue<std::string> output_queue;
    std::mutex input_mutex;
    std::mutex output_mutex;
    std::condition_variable input_cv;
    std::condition_variable output_cv;
    std::atomic<bool> finished_reading{false};
    std::atomic<bool> finished_processing{false};
    std::atomic<size_t> lines_processed{0};
    
    // Check if file is gzipped
    bool is_gzipped(const std::string& filename) {
        std::ifstream file(filename, std::ios::binary);
        if (!file) return false;
        
        unsigned char magic[2];
        file.read(reinterpret_cast<char*>(magic), 2);
        return (magic[0] == 0x1f && magic[1] == 0x8b);
    }
    
    // Load sample names from file
    void load_samples() {
        std::ifstream file(sample_file);
        if (!file) {
            throw std::runtime_error("Cannot open sample file: " + sample_file);
        }
        
        std::string sample;
        while (std::getline(file, sample)) {
            // Remove whitespace
            sample.erase(std::remove_if(sample.begin(), sample.end(), ::isspace), sample.end());
            if (!sample.empty()) {
                target_samples.insert(sample);
            }
        }
        
        if (target_samples.empty()) {
            throw std::runtime_error("No samples found in sample file");
        }
        
        std::cout << "Loaded " << target_samples.size() << " target samples" << std::endl;
    }
    
    // Parse header and find sample indices
    std::string process_header(const std::string& header_line) {
        std::istringstream iss(header_line);
        std::string token;
        std::vector<std::string> fields;
        
        while (std::getline(iss, token, '\t')) {
            fields.push_back(token);
        }
        
        // Find FORMAT column (should be at index 8)
        int format_idx = -1;
        for (int i = 0; i < fields.size(); i++) {
            if (fields[i] == "FORMAT") {
                format_idx = i;
                break;
            }
        }
        
        if (format_idx == -1) {
            throw std::runtime_error("FORMAT column not found in header");
        }
        
        // Sample columns start after FORMAT
        std::vector<std::string> output_fields(fields.begin(), fields.begin() + format_idx + 1);
        
        for (int i = format_idx + 1; i < fields.size(); i++) {
            if (target_samples.count(fields[i])) {
                sample_indices.push_back(i);
                output_fields.push_back(fields[i]);
            }
        }
        
        if (sample_indices.empty()) {
            throw std::runtime_error("No matching samples found in VCF header");
        }
        
        std::cout << "Found " << sample_indices.size() << " matching samples out of " 
                  << (fields.size() - format_idx - 1) << " total samples" << std::endl;
        
        // Reconstruct header line
        std::ostringstream oss;
        for (int i = 0; i < output_fields.size(); i++) {
            if (i > 0) oss << "\t";
            oss << output_fields[i];
        }
        
        return oss.str();
    }
    
    // Process a data line
    std::string process_data_line(const std::string& line) {
        std::istringstream iss(line);
        std::string token;
        std::vector<std::string> fields;
        fields.reserve(50); // Reserve space to avoid reallocations
        
        while (std::getline(iss, token, '\t')) {
            fields.push_back(std::move(token));
        }
        
        if (fields.size() < 9) {
            return line; // Invalid line, return as-is
        }
        
        // Build output line efficiently
        std::ostringstream oss;
        // First 9 columns (up to and including FORMAT)
        for (int i = 0; i < 9; i++) {
            if (i > 0) oss << "\t";
            oss << fields[i];
        }
        
        // Add selected sample columns
        for (int idx : sample_indices) {
            oss << "\t";
            if (idx < fields.size()) {
                oss << fields[idx];
            } else {
                oss << "."; // Missing data
            }
        }
        
        return oss.str();
    }
    
    // Reader thread - reads from file and feeds work queue
    void reader_thread() {
        try {
            if (is_gzipped(input_file)) {
                read_gz_stream();
            } else {
                read_regular_stream();
            }
        } catch (const std::exception& e) {
            std::cerr << "Reader error: " << e.what() << std::endl;
        }
        
        finished_reading = true;
        input_cv.notify_all();
    }
    
    void read_gz_stream() {
        gzFile file = gzopen(input_file.c_str(), "rb");
        if (!file) {
            throw std::runtime_error("Cannot open gzipped file: " + input_file);
        }
        
        char buffer[65536];
        while (gzgets(file, buffer, sizeof(buffer))) {
            std::string line(buffer);
            if (!line.empty() && line.back() == '\n') {
                line.pop_back();
            }
            
            // Wait if queue is full
            std::unique_lock<std::mutex> lock(input_mutex);
            input_cv.wait(lock, [this] { return input_queue.size() < MAX_QUEUE_SIZE; });
            
            input_queue.push(std::move(line));
            lock.unlock();
            input_cv.notify_one();
        }
        
        gzclose(file);
    }
    
    void read_regular_stream() {
        std::ifstream file(input_file);
        if (!file) {
            throw std::runtime_error("Cannot open file: " + input_file);
        }
        
        std::string line;
        while (std::getline(file, line)) {
            // Wait if queue is full
            std::unique_lock<std::mutex> lock(input_mutex);
            input_cv.wait(lock, [this] { return input_queue.size() < MAX_QUEUE_SIZE; });
            
            input_queue.push(std::move(line));
            lock.unlock();
            input_cv.notify_one();
        }
    }
    
    // Worker thread function
    void worker_thread() {
        while (true) {
            std::unique_lock<std::mutex> lock(input_mutex);
            input_cv.wait(lock, [this] { return !input_queue.empty() || finished_reading; });
            
            if (input_queue.empty() && finished_reading) {
                break;
            }
            
            if (input_queue.empty()) {
                continue;
            }
            
            std::string line = std::move(input_queue.front());
            input_queue.pop();
            lock.unlock();
            input_cv.notify_one();
            
            // Process the line
            std::string processed_line;
            if (line.empty() || line[0] == '#') {
                if (line.find("#CHROM") == 0) {
                    processed_line = process_header(line);
                } else {
                    processed_line = std::move(line);
                }
            } else {
                processed_line = process_data_line(line);
            }
            
            // Add to output queue (wait if full)
            std::unique_lock<std::mutex> output_lock(output_mutex);
            output_cv.wait(output_lock, [this] { return output_queue.size() < MAX_QUEUE_SIZE; });
            
            output_queue.push(std::move(processed_line));
            output_lock.unlock();
            output_cv.notify_one();
            
            lines_processed++;
            if (lines_processed % 10000 == 0) {
                std::cout << "Processed " << lines_processed << " lines\r" << std::flush;
            }
        }
    }
    
    // Writer thread - writes output as it becomes available
    void writer_thread() {
        try {
            if (compress_output) {
                write_gz_stream();
            } else {
                write_regular_stream();
            }
        } catch (const std::exception& e) {
            std::cerr << "Writer error: " << e.what() << std::endl;
        }
    }
    
    void write_gz_stream() {
        gzFile out_file = gzopen(output_file.c_str(), "wb");
        if (!out_file) {
            throw std::runtime_error("Cannot create output file: " + output_file);
        }
        
        while (true) {
            std::unique_lock<std::mutex> lock(output_mutex);
            output_cv.wait(lock, [this] { return !output_queue.empty() || finished_processing; });
            
            if (output_queue.empty() && finished_processing) {
                break;
            }
            
            if (output_queue.empty()) {
                continue;
            }
            
            std::string line = std::move(output_queue.front());
            output_queue.pop();
            lock.unlock();
            output_cv.notify_one();
            
            gzprintf(out_file, "%s\n", line.c_str());
        }
        
        gzclose(out_file);
    }
    
    void write_regular_stream() {
        std::ofstream out_file(output_file);
        if (!out_file) {
            throw std::runtime_error("Cannot create output file: " + output_file);
        }
        
        while (true) {
            std::unique_lock<std::mutex> lock(output_mutex);
            output_cv.wait(lock, [this] { return !output_queue.empty() || finished_processing; });
            
            if (output_queue.empty() && finished_processing) {
                break;
            }
            
            if (output_queue.empty()) {
                continue;
            }
            
            std::string line = std::move(output_queue.front());
            output_queue.pop();
            lock.unlock();
            output_cv.notify_one();
            
            out_file << line << "\n";
        }
    }
    
public:
    VCFSampleFilter(const std::string& input, const std::string& output, 
                   const std::string& samples, bool compress = false, int threads = 1)
        : input_file(input), output_file(output), sample_file(samples), 
          compress_output(compress), num_threads(threads) {}
    
    void filter() {
        std::cout << "Loading samples..." << std::endl;
        load_samples();
        
        std::cout << "Starting streaming filter with " << num_threads << " worker threads..." << std::endl;
        
        // Start reader thread
        std::thread reader(&VCFSampleFilter::reader_thread, this);
        
        // Start worker threads
        std::vector<std::thread> workers;
        for (int i = 0; i < num_threads; i++) {
            workers.emplace_back(&VCFSampleFilter::worker_thread, this);
        }
        
        // Start writer thread
        std::thread writer(&VCFSampleFilter::writer_thread, this);
        
        // Wait for reader to finish
        reader.join();
        
        // Wait for all workers to finish
        for (auto& worker : workers) {
            worker.join();
        }
        
        // Signal writer to finish
        finished_processing = true;
        output_cv.notify_all();
        
        // Wait for writer to finish
        writer.join();
        
        std::cout << "\nFiltering complete! Processed " << lines_processed << " lines" << std::endl;
    }
};

void print_usage(const char* program_name) {
    std::cout << "Usage: " << program_name << " [options]\n"
              << "Options:\n"
              << "  -i, --input FILE      Input VCF file (.vcf or .vcf.gz)\n"
              << "  -o, --output FILE     Output VCF file\n"
              << "  -s, --samples FILE    File containing sample names (one per line)\n"
              << "  -z, --compress        Compress output with gzip\n"
              << "  -t, --threads NUM     Number of threads (default: 1)\n"
              << "  -h, --help           Show this help message\n";
}

int main(int argc, char* argv[]) {
    std::string input_file, output_file, sample_file;
    bool compress_output = false;
    int num_threads = 1;
    
    // Parse command line arguments
    for (int i = 1; i < argc; i++) {
        std::string arg = argv[i];
        
        if (arg == "-h" || arg == "--help") {
            print_usage(argv[0]);
            return 0;
        } else if (arg == "-i" || arg == "--input") {
            if (i + 1 < argc) {
                input_file = argv[++i];
            } else {
                std::cerr << "Error: " << arg << " requires a filename" << std::endl;
                return 1;
            }
        } else if (arg == "-o" || arg == "--output") {
            if (i + 1 < argc) {
                output_file = argv[++i];
            } else {
                std::cerr << "Error: " << arg << " requires a filename" << std::endl;
                return 1;
            }
        } else if (arg == "-s" || arg == "--samples") {
            if (i + 1 < argc) {
                sample_file = argv[++i];
            } else {
                std::cerr << "Error: " << arg << " requires a filename" << std::endl;
                return 1;
            }
        } else if (arg == "-z" || arg == "--compress") {
            compress_output = true;
        } else if (arg == "-t" || arg == "--threads") {
            if (i + 1 < argc) {
                num_threads = std::stoi(argv[++i]);
                if (num_threads < 1) {
                    std::cerr << "Error: Number of threads must be positive" << std::endl;
                    return 1;
                }
            } else {
                std::cerr << "Error: " << arg << " requires a number" << std::endl;
                return 1;
            }
        } else {
            std::cerr << "Error: Unknown option " << arg << std::endl;
            print_usage(argv[0]);
            return 1;
        }
    }
    
    // Check required arguments
    if (input_file.empty() || output_file.empty() || sample_file.empty()) {
        std::cerr << "Error: Input file, output file, and sample file are required" << std::endl;
        print_usage(argv[0]);
        return 1;
    }
    
    try {
        VCFSampleFilter filter(input_file, output_file, sample_file, compress_output, num_threads);
        filter.filter();
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }
    
    return 0;
}
