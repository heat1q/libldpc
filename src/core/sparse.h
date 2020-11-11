#pragma once

#include <vector>
#include <forward_list>
#include <iostream>
#include <sstream>
#include <map>

namespace ldpc
{
    template <typename T>
    struct edge
    {
        int rowIndex;
        int colIndex;
        T value;
    };

    struct node
    {
        int nodeIndex;
        int edgeIndex;
    };

    /**
     * @brief Sparse matrix with minimal functionality required for encoding/decoding
     * LDPC codes and other linear block codes.
     * 
     * @tparam T finite field
     */
    template <typename T>
    class sparse_csr
    {
    public:
        sparse_csr() = default;
        sparse_csr(const int m, const int n)
            : numCols(n),
              numRows(m),
              colN(numCols, std::vector<node>()),
              rowN(numRows, std::vector<node>())
        {
        }

        void read_from_file(const std::string &filename, int skipLines);

        void multiply_left(const std::vector<T> &left, std::vector<T> &result) const; // for encoding
        std::vector<T> multiply_left(const std::vector<T> &left) const;

        void multiply_right(const std::vector<T> &right, std::vector<T> &result) const; // for syndrome
        std::vector<T> multiply_right(const std::vector<T> &right) const;

        // Number of columns
        int num_cols() const { return numCols; }
        // Number of rows
        int num_rows() const { return numRows; }
        // Row index neighbors of columns.
        const std::vector<std::vector<node>> &col_neighbor() const { return colN; }
        // Column index neighbors of rows.
        const std::vector<std::vector<node>> &row_neighbor() const { return rowN; }
        // Non-zero entries, i.e. edges.
        const std::vector<edge<T>> &nz_entry() const { return nonZeroVals; }
        // Returns true if the matrix is empty.
        bool empty() const { return ((numCols == 0) && (numRows == 0)); }
        int rank() const;

        // operations for gaussian elimination on sparse matrices over gf(2)

        static void swap_rows(std::vector<std::forward_list<int>> &checkNodeN, std::vector<std::forward_list<int>> &varNodeN, const int i, const int j);
        static void swap_cols(std::vector<std::forward_list<int>> &checkNodeN, std::vector<std::forward_list<int>> &varNodeN, const int i, const int j);
        static void add_rows(std::vector<std::forward_list<int>> &checkNodeN, std::vector<std::forward_list<int>> &varNodeN, const int dest, const int src);
        static void add_cols(std::vector<std::forward_list<int>> &checkNodeN, std::vector<std::forward_list<int>> &varNodeN, const int dest, const int src);
        static void zero_row(std::vector<std::forward_list<int>> &checkNodeN, std::vector<std::forward_list<int>> &varNodeN, const int m);
        static void zero_col(std::vector<std::forward_list<int>> &checkNodeN, std::vector<std::forward_list<int>> &varNodeN, const int n);

    private:
        int numCols;                         // number of columns
        int numRows;                         // number of rows
        std::vector<std::vector<node>> colN; // column neighbors, with row and edge index
        std::vector<std::vector<node>> rowN; // row neigbors, with col and edge index
        std::vector<edge<T>> nonZeroVals;    // edges, i.e. non-zero entries with row and col index
    };

    /**
     * @brief Read a sparse CSR file.
     * 
     * @throw runtime_error
     * @tparam T finite field
     * @param filename Input file
     * @param skipLines Number of lines to skip from head
     */
    template <typename T>
    void sparse_csr<T>::read_from_file(const std::string &filename, int skipLines)
    {
        int index = 0;

        std::ifstream infile(filename);
        std::string line;

        if (!infile.good())
            throw std::runtime_error("can not open file for reading");

        std::map<int, std::vector<node>> cols;
        std::map<int, std::vector<node>> rows;

        numCols = 0;
        numRows = 0;
        nonZeroVals.clear();

        // skip lines at beginning of file
        while (skipLines-- > 0)
        {
            std::getline(infile, line);
        }

        while (getline(infile, line))
        {
            edge<T> entry;

            std::istringstream record(line);
            record >> entry.rowIndex;
            record >> entry.colIndex;
            record >> entry.value;

            // no value is given
            if (entry.value == 0)
            {
                entry.value = T(1);
            }

            nonZeroVals.push_back(entry);

            cols[entry.colIndex].push_back(node({entry.rowIndex, index}));
            rows[entry.rowIndex].push_back(node({entry.colIndex, index}));

            // find the number of columns and rows from indices
            numCols = std::max(numCols, entry.colIndex);
            numRows = std::max(numRows, entry.rowIndex);

            ++index;
        }

        ++numCols;
        ++numRows;

        colN = std::vector<std::vector<node>>(numCols, std::vector<node>());
        rowN = std::vector<std::vector<node>>(numRows, std::vector<node>());

        // copy assign the elements of the map to the vector
        for (auto e : cols)
            colN[e.first] = e.second;
        for (auto e : rows)
            rowN[e.first] = e.second;
    }

    /**
     * @brief Multiply vector from left handside with matrix over field T
     * 
     * @tparam T finite field
     * @param left Row vector
     * @param result Result row vector
     */
    template <typename T>
    void sparse_csr<T>::multiply_left(const std::vector<T> &left, std::vector<T> &result) const
    {
        for (int j = 0; j < numCols; ++j)
        {
            for (const auto &n : colN[j])
            {
                result[j] += left[n.nodeIndex] * nonZeroVals[n.edgeIndex].value;
            }
        }
    }

    /**
     * @brief Multiply vector from left handside with matrix over field T
     * 
     * @tparam T finite field
     * @param left Row vector
     * @return std::vector<T> Result row vector
     */
    template <typename T>
    std::vector<T> sparse_csr<T>::multiply_left(const std::vector<T> &left) const
    {
        std::vector<T> result(numCols);
        multiply_left(left, result);
        return result;
    }

    /**
     * @brief Multiply vector from right handside with matrix over field T
     * 
     * @tparam T finite field
     * @param right Column vector
     * @param result Result column vector
     */
    template <typename T>
    void sparse_csr<T>::multiply_right(const std::vector<T> &right, std::vector<T> &result) const
    {
        for (int i = 0; i < numRows; ++i)
        {
            for (const auto &n : rowN[i])
            {
                result[i] += right[n.nodeIndex] * nonZeroVals[n.edgeIndex].value;
            }
        }
    }

    /**
     * @brief Multiply vector from right handside with matrix over field T
     * 
     * @tparam T finite field
     * @param right Column vector
     * @return std::vector<T> Result column vector
     */
    template <typename T>
    std::vector<T> sparse_csr<T>::multiply_right(const std::vector<T> &right) const
    {
        std::vector<T> result(numRows);
        multiply_right(right, result);
        return result;
    }

    /**
     * @brief Calculate rank over gf(2)
     * 
     * @tparam T finite field gf(2)
     * @return int Rank
     */
    template<typename T>
    int sparse_csr<T>::rank() const
    {
        int rank = num_cols();

        std::vector<std::forward_list<int>> checkNodeN(num_rows()); // direct check node neighbourhood with node oriented index
        std::vector<std::forward_list<int>> varNodeN(num_cols());// direct variable node neighbourhood with node oriented index

        for (int i = 0; i < num_rows(); ++i)
        {
            for (auto r : row_neighbor()[i]) 
                checkNodeN[i].push_front(r.nodeIndex);
        }

        for (int j = 0; j < num_cols(); ++j)
        {
            for (auto c : col_neighbor()[j]) 
                varNodeN[j].push_front(c.nodeIndex);
        }
        
        for (int row = 0; row < rank; ++row)
        {
            // check what value h[row][row] has
            auto it = std::find(varNodeN[row].begin(), varNodeN[row].end(), row);
            if (it != varNodeN[row].end()) // values is non-zero
            {
                // now add current row to all rows where a non-zero entry is in the current col, to remove 1
                auto vn = varNodeN[row];
                for (auto ri : vn)
                {
                    if (ri > row)
                    {
                        sparse_csr<T>::add_rows(checkNodeN, varNodeN, ri, row);
                    }
                }
            }
            else // value is zero
            {
                // if there is a row below it with non-zero entry in same col, swap current rows
                bool isZero = true;
                // find first row with non-zero entry
                for (auto ri : varNodeN[row])
                {
                    if (ri > row)
                    {
                        sparse_csr<T>::swap_rows(checkNodeN, varNodeN, ri, row);
                        isZero = false;
                        break;
                    }
                }

                // if all elements in current col below h[row][row] are zero, swap col it with rank-1 col
                if (isZero)
                {
                    --rank;
                    // copy last col
                    sparse_csr<T>::zero_col(checkNodeN, varNodeN, row);
                    sparse_csr<T>::add_cols(checkNodeN, varNodeN, row, rank);
                }

                --row;
            }
        }
        
        return rank;
    }

    template<typename T>
    void sparse_csr<T>::swap_rows(std::vector<std::forward_list<int>> &checkNodeN, std::vector<std::forward_list<int>> &varNodeN, const int first, const int second)
    {
        sparse_csr<T>::add_rows(checkNodeN, varNodeN, first, second);
        sparse_csr<T>::add_rows(checkNodeN, varNodeN, second, first);
        sparse_csr<T>::add_rows(checkNodeN, varNodeN, first, second);
    }

    template<typename T>
    void sparse_csr<T>::swap_cols(std::vector<std::forward_list<int>> &checkNodeN, std::vector<std::forward_list<int>> &varNodeN, const int first, const int second)
    {
        sparse_csr<T>::add_cols(checkNodeN, varNodeN, first, second);
        sparse_csr<T>::add_cols(checkNodeN, varNodeN, second, first);
        sparse_csr<T>::add_cols(checkNodeN, varNodeN, first, second);
    }

    template<typename T>
    void sparse_csr<T>::add_rows(std::vector<std::forward_list<int>> &checkNodeN, std::vector<std::forward_list<int>> &varNodeN, const int dest, const int src)
    {
        auto &row = checkNodeN[dest];
        for (auto vn : checkNodeN[src]) // check vn index of src row
        {
            auto it = std::find(row.begin(), row.end(), vn); // dest already contains current vn index?
            if (it == row.end()) // gf2 addition
            {
                row.push_front(vn); // add var index to row
                varNodeN[vn].push_front(dest); // add row index to var
            }
            else // same index means removal
            {
                row.remove(vn); 
                varNodeN[vn].remove(dest);
            }
        }
    }

    template<typename T>
    void sparse_csr<T>::add_cols(std::vector<std::forward_list<int>> &checkNodeN, std::vector<std::forward_list<int>> &varNodeN, const int dest, const int src)
    {
        auto &col = varNodeN[dest];
        for (auto cn : varNodeN[src]) // check cn index of src col
        {
            auto it = std::find(col.begin(), col.end(), cn); // dest already contains current cn index?
            if (it == col.end()) // gf2 addition
            {
                col.push_front(cn); // add cn index to col
                checkNodeN[cn].push_front(dest); // add col index to cn
            }
            else // same index means removal
            {
                col.remove(cn); 
                checkNodeN[cn].remove(dest);
            }
        }
    }

    template<typename T>
    void sparse_csr<T>::zero_row(std::vector<std::forward_list<int>> &checkNodeN, std::vector<std::forward_list<int>> &varNodeN, const int m)
    {
        for (auto vn : checkNodeN[m]) // from selected row, for each vn index, remove m from vn
        {
            varNodeN[vn].remove(m);
        }
        checkNodeN[m].clear();
    }

    template<typename T>
    void sparse_csr<T>::zero_col(std::vector<std::forward_list<int>> &checkNodeN, std::vector<std::forward_list<int>> &varNodeN, const int n)
    {
        for (auto cn : varNodeN[n]) // from selected col, for each cn index, remove n from cn
        {
            checkNodeN[cn].remove(n);
        }
        varNodeN[n].clear();
    }
} // namespace ldpc
