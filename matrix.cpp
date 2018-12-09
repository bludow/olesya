#include <cassert>
#include <cmath>
#include <iostream>
#include <vector>

/// A matrix row.
typedef std::vector<double> Row;
/// A whole matrix.
typedef std::vector<Row> Matrix;

namespace {
    /// Subtract matrix rows: a -= c * b
    /// Row sizes must be equal.
    void subtract(Row& a, const Row& b, double c) {
        assert(a.size() == b.size());
        for (size_t i = 0; i < a.size(); ++i) {
            a[i] -= b[i] * c;
        }
    }

    /// Multiply matrix row by a constant (in-place)
    void multiply(Row& r, double c) {
        for (size_t i = 0; i < r.size(); ++i) {
            r[i] *= c;
        }
    }
}


/// Compute an inverse matrix by using the Gauss method.
/// This algorithm tries to be as stable as possible by choosing
/// always the row with the largest element for subtraction.
class GaussSolver {

public:
    /// Construct the Gauss method solver for a given matrix.
    GaussSolver(Matrix&& m)
        : mat(m)
        , unit()
    {
        for (size_t i = 0; i < mat.size(); ++i) {
            // Check if the matrix is a square one
            assert(mat[i].size() == mat.size());
            // Fill unit matrix
            unit.push_back(Row(mat.size(), 0.0));
            unit[i][i] = 1.0;
        }
    }

    /// Compute an inverse matrix and return it.
    Matrix invert() {
        gauss_down();
        gauss_up();
        return unit;
    }

private:

    // Copying and assignment is prohibited
    GaussSolver(const GaussSolver&) = delete;
    GaussSolver& operator= (const GaussSolver&) = delete;

private:

    /// Swap matrix rows i and j. Do the same with unit matrix.
    void swap_rows(size_t i, size_t j) {
        mat[i].swap(mat[j]);
        unit[i].swap(unit[j]);
    }

    /// Find the row with maximal module element in rows [a, b)
    /// and swap its row with ith row.
    void find_max(size_t i, size_t a, size_t b) {
        size_t mj = i;
        double m = fabs(mat[i][i]);
        for (size_t j = a; j < b; ++j) {
            double x = fabs(mat[j][i]);
            if (x > m) {
                m = x;
                mj = j;
            }
        }
        if (mj != i) {
            swap_rows(mj, i);
        }
    }

    /// Subtract ith row from all rows in [a, b)
    /// to zero ith element of them. Copy same operations
    /// on the unit matrix.
    void subtract_all(size_t i, size_t a, size_t b) {
        double relem = 1.0 / mat[i][i];
        multiply(mat[i], relem);
        multiply(unit[i], relem);
        for (size_t j = a; j < b; ++j) {
            double coeff = mat[j][i];
            subtract(mat[j], mat[i], coeff);
            subtract(unit[j], unit[i], coeff);
        }

    }
    /// Gauss method top-to-bottom pass.
    /// Duplicate on unit matrix.
    void gauss_down() {
        for (size_t i = 0; i < mat.size(); ++i) {
            find_max(i, i + 1, mat.size());
            subtract_all(i, i + 1, mat.size());
        }
    }

    /// Gauss method bottom-to-top pass.
    /// Duplicate on unit matrix.
    void gauss_up() {
        for (size_t i = mat.size(); i > 0; --i) {
            find_max(i - 1, 0, i - 1);
            subtract_all(i - 1, 0, i - 1);
        }
    }

private:
    /// The matrix we work on.
    Matrix mat;
    /// The unit matrix we duplicate operations on.
    Matrix unit;
};


int main() {
    // A simple test for 3x3 matrix.
    Matrix m;
    for (size_t i = 0; i < 3; ++i) {
        m.push_back(Row(3, 0.0));
        m[i][i] = 2.0;
    }

    GaussSolver solver(std::move(m));
    Matrix r = solver.invert();
    for (size_t i = 0; i < r.size(); ++i) {
        for (size_t j = 0; j < r[i].size(); ++j) {
            std::cout << "\t" << r[i][j];
        }
        std::cout << std::endl;
    }
    std::cout << std::endl;

    return 0;
}
