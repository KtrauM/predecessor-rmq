#include <stdint.h>

#include <fstream>
#include <iostream>
#include <sstream>
#include <vector>
#include <bit>
#include <unordered_map>

#define st first
#define nd second

using namespace std;

class NaiveRMQ {
   public:
    NaiveRMQ(vector<int64_t> &input) {
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
    int64_t query(size_t start, size_t end) {
        return table[start][end];
    }

   private:
    vector<vector<int64_t>> table;
};

class SparseTableRMQ {
    public:
    SparseTableRMQ(vector<int64_t> &input) : values(input) {
        size_t N = input.size();
        size_t logN = bit_width(N);
        table.resize(N);
        for (size_t i = 0; i < N; i++) {
            table[i].resize(logN);
            table[i][0] = i;
        }
        for (size_t i = 0; i < N; i++) {
            for (size_t j = 1; j < logN; j++) {
                size_t offset = 1 << (j - 1);
                size_t right_index = min(N - 1, i + offset);
                if (input[table[i][j - 1]] < input[table[right_index][j - 1]]) {
                    table[i][j] = table[i][j - 1];
                } else {
                    table[i][j] = table[right_index][j - 1];
                }
            }
        }    
    }

    int64_t query(size_t start, size_t end) {
        size_t level = bit_width(end - start) - 1;
        size_t offset = 1 << level;
        size_t left_rmq = table[start][level];
        size_t right_rmq = table[end - offset][level];
        return values[left_rmq] < values[right_rmq] ? left_rmq : right_rmq;
    }

   private:
    vector<int64_t> values;
    vector<vector<int64_t>> table;
};

struct Node{
    int64_t value;
    Node *left, *right;
    Node *parent;

    Node(int64_t x) {
        value = x;
        parent = NULL;
        left = NULL;
        right = NULL;
    }
};

class LinearSpaceRMQ {
    public:
    LinearSpaceRMQ(vector<int64_t> input) : data(input) {
        block_size = bit_width(input.size()) / 4 + 1;
        int64_t block_min = INT64_MAX;
        size_t block_idx = 0;
        size_t j = 0;
        vector<int64_t> tree_stack;
        uint64_t tree_number = 1;
        while (block_idx * block_size + j < input.size()) {
            printf("Input: %d, block_idx: %d, block_size: %d\n", input.size(), block_idx, block_size);
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
                    block_idx = cur_idx;
                }
                j++;
            }
            // padding for the last block's tree number
            while (j++ < block_size) {
                tree_number <<= 1;
            }

            if (cartesian_trees.find(tree_number) == cartesian_trees.end()) {
                cartesian_trees[tree_number] = vector<vector<size_t>>(block_size);
                for (int start = 0; start < block_size; start++) {
                    cartesian_trees[tree_number][start] = vector<size_t>(block_size);
                    cartesian_trees[tree_number][start][start] = block_size * block_idx + start;
                    for (int end = start + 1; end < block_size; end++) {
                        size_t prev_min_idx = cartesian_trees[tree_number][start][end - 1];
                        size_t cur_min_idx = block_size * block_idx + end;
                        if (input[prev_min_idx] > input[cur_min_idx]) {
                            cartesian_trees[tree_number][start][end] = cur_min_idx;
                        } else {
                            cartesian_trees[tree_number][start][end] = prev_min_idx;
                        }
                    }
                }
            }


            block_tree_number.push_back(tree_number);

            tree_number = 1;
            tree_stack.clear();
            
            block_min = INT64_MAX;

            j = 0;
            block_idx++;
        }
        vector<int64_t> block_vals;
        for (auto block: blocks) {
            block_vals.push_back(input[block]);
        }
        sparse_table = new SparseTableRMQ(block_vals);
    }

    size_t query(size_t start, size_t end) {
        size_t first_block = start / block_size;
        size_t first_block_start = first_block * block_size + start % block_size;
        size_t first_block_end = (first_block + 1) * block_size - 1;
        size_t last_block = end / block_size; 
        size_t last_block_start = last_block * block_size;
        size_t last_block_end = last_block * block_size + end % block_size;
        size_t first_full_block = first_block_start % block_size == 0 ? first_block : first_block + 1;
        size_t last_full_block = last_block_end % block_size == 0 ? last_block - 1 : last_block;

        size_t min_idx = sparse_table->query(first_full_block, last_full_block);
        int64_t min_val = data[min_idx];

        size_t left_partial_min_idx = cartesian_trees[block_tree_number[first_block]][first_block_start][first_block_end];
        size_t right_partial_min_idx = cartesian_trees[block_tree_number[last_block]][last_block_start][last_block_end];

        if (min_val > data[left_partial_min_idx]) {
            min_idx = left_partial_min_idx;
            min_val = data[left_partial_min_idx];
        }

        if (min_val > data[right_partial_min_idx]) {
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
    vector<int64_t> data;
    SparseTableRMQ *sparse_table;
};


int main(int argc, char *argv[]) {
    string query_type = argv[1];
    string input_path = argv[2];
    string output_path = argv[3];

    ifstream input_file(input_path);
    ofstream output_file(output_path);

    int64_t n;
    input_file >> n;

    vector<int64_t> input(n);

    for (size_t i = 0; i < n; i++) {
        input_file >> input[i];
    }

    if (query_type == "pd") {
        vector<pair<int64_t, int64_t>> query;

    } else if (query_type == "rmq") {
        string line;
        getline(input_file, line);  // Read the remaining newline character

        vector<pair<int, int>> query;

        while (getline(input_file, line)) {
            stringstream ss(line);
            string token;

            pair<int, int> q;

            getline(ss, token, ',');
            q.st = stoull(token);

            getline(ss, token, ',');
            q.nd = stoull(token);

            query.push_back(q);
        }

        NaiveRMQ *naive = new NaiveRMQ(input);
        SparseTableRMQ *sparse = new SparseTableRMQ(input);
        LinearSpaceRMQ *linear = new LinearSpaceRMQ(input);
        vector<int64_t> v = {8, 2, 5, 1, 9, 11, 10, 20, 22, 4};

        for (pair<int64_t, int64_t> q : query) {
            printf("%ld\n", linear->query(q.st, q.nd));
            output_file << linear->query(q.st, q.nd) << endl;
        }
    }

    input_file.close();
    output_file.close();
    // printf("RESULT algo=%s name%s time=%d space=%d", query_type);

    return 0;
}