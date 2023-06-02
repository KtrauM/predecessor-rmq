#include <stdint.h>

#include <fstream>
#include <iostream>
#include <sstream>
#include <vector>
#include <bit>

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

        for (pair<int64_t, int64_t> q : query) {
            printf("%ld\n", sparse->query(q.st, q.nd));
            output_file << sparse->query(q.st, q.nd) << endl;
        }
    }

    input_file.close();
    output_file.close();
    // printf("RESULT algo=%s name%s time=%d space=%d", query_type);

    return 0;
}