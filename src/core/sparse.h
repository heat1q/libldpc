#pragma once

#include <vector>
#include <iostream>
#include <sstream>

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

        std::vector<T> multiply_left(const std::vector<T> &left) const;   // for encoding
        std::vector<T> multiply_right(const std::vector<T> &right) const; // for syndrome

        const int num_cols() const { return numCols; }
        const int num_rows() const { return numRows; }
        const std::vector<std::vector<node>> &col_neighbor() const { return colN; }
        const std::vector<std::vector<node>> &row_neighbor() const { return rowN; }
        const std::vector<edge<T>> &nz_entry() const { return nonZeroVals; }

    private:
        int numCols;                         // number of columns
        int numRows;                         // number of rows
        std::vector<std::vector<node>> colN; // column neighbors, with row and edge index
        std::vector<std::vector<node>> rowN; // row neigbors, with col and edge index
        std::vector<edge<T>> nonZeroVals;    // edges, i.e. non-zero entries with row and col index
    };

    template <typename T>
    void sparse_csr<T>::read_from_file(const std::string &filename, int skipLines)
    {
        int index = 0;

        std::ifstream infile(filename);
        std::string line;

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

            colN[entry.colIndex].push_back(node({entry.rowIndex, index}));
            rowN[entry.rowIndex].push_back(node({entry.colIndex, index}));

            ++index;
        }
    }

    /**
     * @brief Multiply vector from left handside with matrix over field T
     * 
     * @tparam T finite field
     * @param left Row vector
     * @return std::vector<T> 
     */
    template <typename T>
    std::vector<T> sparse_csr<T>::multiply_left(const std::vector<T> &left) const
    {
        std::vector<T> result(numCols);

        for (int j = 0; j < numCols; ++j)
        {
            for (const auto &n : colN[j])
            {
                result[j] += left[n.nodeIndex] * nonZeroVals[n.edgeIndex].value;
            }
        }

        return result;
    }

    /**
     * @brief Multiply vector from right handside with matrix over field T
     * 
     * @tparam T finite field
     * @param right Column vector
     * @return std::vector<T> 
     */
    template <typename T>
    std::vector<T> sparse_csr<T>::multiply_right(const std::vector<T> &right) const
    {
        std::vector<T> result(numRows);

        for (int i = 0; i < numRows; ++i)
        {
            for (const auto &n : rowN[i])
            {
                result[i] += right[n.nodeIndex] * nonZeroVals[n.edgeIndex].value;
            }
        }

        return result;
    }
} // namespace ldpc
