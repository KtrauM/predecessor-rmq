#include <stdint.h>

#include <bit>
#include <fstream>
#include <iostream>
#include <sstream>
#include <unordered_map>
#include <vector>
#include <deque> 
#include <chrono>
#include <cmath>

#define st first
#define nd second

using namespace std;

class BitVector {
    public:
    BitVector(vector<bool> &input) {
        size_t N = input.size();
        size_t logN = bit_width(N) - 1;
        block_size = logN / 2;
        super_block_size = block_size * block_size;

        block.resize(N / block_size);
        block_number.resize(N / block_size);
        super_block.resize(N / super_block_size);
        
        vector<uint32_t> inblock;
        uint32_t inblock_zeroes = 0;
        size_t cur_block = 0;
        size_t cur_super_block = 0;
        for (int i = 0; i < N; i++) {
            block_number[cur_block] += input[i];
            block_number[cur_block] <<= 1;
            inblock_zeroes += input[i] ? 1 : 0;
            inblock.push_back(inblock_zeroes);

            if (input[i] == 0) {
                block[cur_block]++;
                super_block[cur_super_block]++;
            }
            
            if ((i + 1) % block_size == 0) {
                lookup[block_number[cur_block]] = inblock;
                inblock_zeroes = 0;
                inblock.clear();
                cur_block++;
                block[cur_block] = block[cur_block - 1];
            }
            
            if ((i + 1) % super_block_size == 0) {
                cur_super_block++;
                block[cur_block] = 0;
            }
        }
    }

    uint32_t rank(uint32_t x, uint32_t bit) {
        uint32_t super_idx = x / super_block_size;
        uint32_t block_idx = x / block_size;
        uint32_t block_offset = x % block_size;
        uint32_t zeroes = super_block[super_idx] + block[block_idx] + lookup[block_number[block_idx]][block_offset];
        return bit ? x - zeroes : zeroes;
    }

    int select(int x, int bit) {
        
    }

    private:
    size_t block_size, super_block_size;
    vector<uint32_t> block, block_number, super_block;
    unordered_map<uint32_t, vector<uint32_t>> lookup;


};

class NaiveBitVector {
    public:
    NaiveBitVector(vector<bool> &input) : v(input) {}
    
    uint64_t rank(uint64_t x, uint64_t bit) {
        uint64_t cnt = 0;
        for (uint64_t i = 0; i < x; i++) {
            if (!v[i]) {
                cnt++;
            }
        }
        return !bit ? cnt : x - cnt;
    }

    uint64_t select(uint64_t x, uint64_t bit) {
        uint64_t cnt = 0;
        if (x == 0) {
            return 0;
        }
        for (uint64_t i = 0; i < v.size(); i++) {
            if (v[i] == bit) {
                cnt++;
            }
            if (cnt == x) {
                return i;
            }
        }
        return v.size() - 1;
    }
    private:
    vector<bool> v;
};

class EliasFano {
    public:
    EliasFano(vector<uint64_t> &input) {
        size_t N = input.size();
        uint64_t U = input[N - 1];
        size_t upper_size = ceil(log2(N));
        lower_size = ceil(log2(U) - log2(N));

        lower_mask = (1 << lower_size) - 1;
        upper_mask = ((1 << (upper_size + lower_size)) - 1) ^ lower_mask;

        vector<bool> upper_bits;
        uint64_t current_bucket = 0;
        for (size_t i = 0; i < N; i++) {
            uint64_t upper_part = (input[i] & upper_mask) >> lower_size;
            while (current_bucket != upper_part) {
                upper_bits.push_back(false);
                current_bucket++;
            }
            upper_bits.push_back(true);
            lower.push_back(input[i] & lower_mask);
        }
        upper_bits.resize(2 * N);
        upper = new NaiveBitVector(upper_bits);
    }

    uint64_t access(int64_t index) {
        uint64_t upper_bits = upper->select(index + 1, 1) - index;
        return (upper_bits << lower_size) + lower[index];
    }

    uint64_t predecessor(uint64_t x) {
        uint64_t upper_x = (x & upper_mask) >> lower_size;
        uint64_t lower_x = x & lower_mask;
        int64_t start = upper_x ? upper->select(upper_x, 0) - upper_x + 1 : 0;
        int64_t end = upper->select(upper_x + 1, 0) - upper_x;
        for (int64_t i = end; i >= start; i--) {
            uint64_t candidate = access(i);
            if (candidate <= x) {
                return candidate;
            }
        }
        if (start > 0) {
            return access(start - 1);
        }
        return UINT64_MAX;
    }

    private:
    uint64_t lower_size, lower_mask, upper_mask;
    vector<uint64_t> lower;
    NaiveBitVector *upper;

};

class NaiveRMQ {
   public:
    NaiveRMQ(vector<uint64_t> &input) {
        size_t N = input.size();
        table.resize(N);
        for (size_t i = 0; i < N; i++) {
            table[i].resize(N);
            table[i][i] = i;
            for (size_t j = i + 1; j < N; j++) {
                if (input[table[i][j - 1]] > input[j]) {
                    table[i][j] = j;
                } else {
                    table[i][j] = table[i][j - 1];
                }
            }
        }
    }
    size_t query(size_t start, size_t end) { return table[start][end]; }

   private:
    vector<vector<size_t>> table;
};

class SparseTableRMQ {
   public:
    SparseTableRMQ(vector<uint64_t> &input) : values(input) {
        size_t N = input.size();
        size_t logN = bit_width(N);
        table.resize(N);
        for (size_t i = 0; i < N; i++) {
            table[i].resize(logN);
            table[i][0] = i;
        }
        for (size_t j = 1; j < logN; j++) {
            for (size_t i = 0; i < N; i++) {
                size_t offset = 1ull << (j - 1);
                size_t right_index = min(N - 1, i + offset);
                if (input[table[i][j - 1]] <= input[table[right_index][j - 1]]) {
                    table[i][j] = table[i][j - 1];
                } else {
                    table[i][j] = table[right_index][j - 1];
                }
            }
        }
    }

    size_t query(size_t start, size_t end) {
        if (start == end) {
            return start;
        }
        size_t level = bit_width(end - start + 1) - 1;
        size_t offset = level == 0 ? 0 : 1 << level;
        size_t left_rmq = table[start][level];
        size_t right_rmq = table[end - offset + 1][level];
        return values[left_rmq] <= values[right_rmq] ? left_rmq : right_rmq;
    }

   private:
    vector<uint64_t> values;
    vector<vector<size_t>> table;
};

class LinearSpaceRMQ {
   public:
    LinearSpaceRMQ(vector<uint64_t> input) : data(input) {
        block_size = bit_width(input.size()) / 4 + 1;
        uint64_t block_min = UINT64_MAX;
        size_t block_min_idx = 0;
        size_t block_idx = 0;
        size_t j = 0;
        vector<uint64_t> tree_stack;
        uint64_t tree_number = 1;
        while (block_idx * block_size + j < input.size()) {
            while (j < block_size && block_idx * block_size + j < input.size()) {
                size_t cur_idx = block_idx * block_size + j;
                while (!tree_stack.empty() && tree_stack.back() > input[cur_idx]) {
                    tree_stack.pop_back();
                    tree_number <<= 1;
                }
                tree_stack.push_back(input[cur_idx]);
                tree_number <<= 1;
                tree_number++;

                if (block_min > input[cur_idx]) {
                    block_min = input[cur_idx];
                    block_min_idx = cur_idx;
                }
                j++;
            }
            // padding for the last block's tree number
            while (j++ < block_size) {
                tree_number <<= 1;
            }

            if (cartesian_trees.find(tree_number) == cartesian_trees.end()) {
                cartesian_trees[tree_number] = vector<vector<size_t>>(block_size);
                for (size_t start = 0; start < block_size; start++) {
                    cartesian_trees[tree_number][start] = vector<size_t>(block_size);
                    cartesian_trees[tree_number][start][start] = start;
                    for (size_t end = start + 1; end < block_size; end++) {
                        size_t prev_min_idx = cartesian_trees[tree_number][start][end - 1];
                        size_t cur_min_idx = end;
                        size_t offset = block_size * block_idx;
                        if (input[offset + prev_min_idx] > input[offset + cur_min_idx]) {
                            cartesian_trees[tree_number][start][end] = cur_min_idx;
                        } else {
                            cartesian_trees[tree_number][start][end] = prev_min_idx;
                        }
                    }
                }
            }

            block_tree_number.push_back(tree_number);
            blocks.push_back(block_min_idx);
            tree_number = 1;
            tree_stack.clear();

            block_min = UINT64_MAX;
            block_min_idx = 0;

            j = 0;
            block_idx++;
        }
        vector<uint64_t> block_vals;
        for (size_t block : blocks) {
            block_vals.push_back(input[block]);
        }
        sparse_table = new SparseTableRMQ(block_vals);
    }

    size_t query(size_t start, size_t end) {
        int32_t first_block = start / block_size;
        int32_t first_block_start = start % block_size;
        bool isEndWithinTheSameBlockAsStart = (first_block + 1) * block_size > end;
        int32_t first_block_end = isEndWithinTheSameBlockAsStart ? end % block_size : block_size - 1;
        int32_t last_block = end / block_size;
        int32_t last_block_start = isEndWithinTheSameBlockAsStart ? end % block_size : 0;
        int32_t last_block_end = end % block_size;
        int32_t first_full_block =
            first_block_start == 0 && !isEndWithinTheSameBlockAsStart ? first_block : first_block + 1;
        int32_t last_full_block = last_block - 1;
        int32_t min_idx = first_full_block <= last_full_block
                             ? blocks[sparse_table->query(first_full_block, last_full_block)]
                             : -1;
        uint64_t min_val = min_idx != -1 ? data[min_idx] : UINT64_MAX;

        size_t left_partial_min_idx =
            first_block * block_size +
            cartesian_trees[block_tree_number[first_block]][first_block_start][first_block_end];

        size_t right_partial_min_idx =
            last_block * block_size +
            cartesian_trees[block_tree_number[last_block]][last_block_start][last_block_end];

        if (min_val > data[left_partial_min_idx] ||
            (min_val == data[left_partial_min_idx] && min_idx > left_partial_min_idx)) {
            min_idx = left_partial_min_idx;
            min_val = data[left_partial_min_idx];
        }

        if (min_val > data[right_partial_min_idx] ||
            (min_val == data[right_partial_min_idx] && min_idx > right_partial_min_idx)) {
            min_idx = right_partial_min_idx;
            min_val = data[right_partial_min_idx];
        }

        return min_idx;
    }

   private:
    size_t block_size;
    vector<size_t> blocks;
    unordered_map<uint64_t, vector<vector<size_t>>> cartesian_trees;
    vector<uint64_t> block_tree_number;
    vector<uint64_t> data;
    SparseTableRMQ *sparse_table;
};

int main(int argc, char *argv[]) {
    string query_type = argv[1];
    string input_path = argv[2];
    string output_path = argv[3];

    ifstream input_file(input_path);
    ofstream output_file(output_path);

    size_t n;
    input_file >> n;

    vector<uint64_t> input(n);
    for (size_t i = 0; i < n; i++) {
        input_file >> input[i];
    }

    if (query_type == "pd") {
        vector<uint64_t> query;
        vector<uint64_t> answers;
        uint64_t q;
        while(input_file >> q) {
            query.push_back(q);
        }

        EliasFano *coding = new EliasFano(input);
        for (uint64_t q: query) {
            if (input_path == "access.txt") {
                answers.push_back(coding->access(q));
            } else {
                answers.push_back(coding->predecessor(q));
            }
        }
        for (int i = 0; i < answers.size(); i++) {
            cout << query[i] << " " << answers[i] << endl;
        }

    } else if (query_type == "rmq") {
        string line;
        getline(input_file, line);  // Read the remaining newline character

        vector<pair<size_t, size_t>> query;
        vector<int32_t> answers; 
        while (getline(input_file, line)) {
            stringstream ss(line);
            string token;

            pair<size_t, size_t> q;

            getline(ss, token, ',');
            q.st = stoull(token);

            getline(ss, token, ',');
            q.nd = stoull(token);

            query.push_back(q);
        }

        auto start = std::chrono::high_resolution_clock::now();

        // NaiveRMQ *naive = new NaiveRMQ(input);
        // SparseTableRMQ *sparse = new SparseTableRMQ(input);
        LinearSpaceRMQ *linear = new LinearSpaceRMQ(input);

        for (pair<int64_t, int64_t> q : query) {
            // size_t result = output_path == "naive.txt" ? naive->query(q.st, q.nd)
            //                                            : (output_path == "sparse.txt" ? sparse->query(q.st, q.nd)
            //                                                                           : linear->query(q.st, q.nd));
            size_t result = linear->query(q.st, q.nd);
            answers.push_back(result);
        }

        auto end = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();

        printf("RESULT algo=%s name=murat_kurnaz time=%ld space=%ld\n", query_type.c_str(), duration, sizeof(linear));
        
        
        for (int32_t ans: answers) {
            output_file << ans << endl;
        }
    }

    input_file.close();
    output_file.close();
    return 0;
}