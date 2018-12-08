#include <iostream>
#include <ctime>
#include <blocked_range.h>
#include <parallel_reduce.h>

using namespace std;
using namespace tbb;


namespace {
    const double x_min = 2;
    const double x_max = 6;
    const double y_min = 1;
    const double y_max = 4;

    /// Generate random numbers from `a` to `b`.
    inline double random(const double a, const double b) {
        return a + (b - a) * rand() / RAND_MAX;
    }

    /// Check if the point is inside of the area.
    inline bool is_in_area(const double x, const double y) {
        return ( (x <= 0)
              && (y >= (-1.5 * x) - 1)
              && (y <= -x)
               )
            || ( (x >= 0)
              && (y >= ( 1.5 * x) - 1)
              && (y <=  x)
               );
    }

    /// Perform the area computation within the given count range.
    inline unsigned area_hits(const unsigned begin, const unsigned end) {
        unsigned hits = 0;
        for (unsigned i = begin; i < end; ++i) {
            const double x = random(X_MIN, X_MAX);
            const double y = random(Y_MIN, Y_MAX);
            if (is_in_area(x, y)) {
                ++hits;
            }
        }
        return hits;
    }

    /// Compute the area for given hits.
    inline double area_of(const unsigned hits, const unsigned total) {
        return (x_max - x_min) * (y_max - y_min) * hits / total;
    }

} // local namespace


/// Compute of the area covered by function by Monte-Carlo method.
class AreaComputer {

public:

    /// Initialize the computation for N points.
    Computer(const unsigned n)
        : N(n)
        , sum(0)
    { }

    /// Initialize the computation for N points by splitting a bigger one.
    Computer(const Computer c, const unsigned split)
        : N(n)
        , sum(0)
    { }

    /// Perform the actual computation (the functor).
    void operator ()(const blocked_range<unsigned>& r) {
        sum = area_hits(r.begin(), r.end(), N);
    }

    /// Join results of two computations.
    void join(const AreaComputer& other) {
        sum += other.sum;
    }

    /// Get the total result.
    double result() const {
        return sum;
    }

private:

    /// Point count.
    const unsigned N;
    /// Area sum.
    double sum;

};


int main() {
    const unsigned N = 5000;
    srand(time(nullptr));

    // Method 1: direct computation.
    {
        int area = round(area_of(area_hits(0, N), N));
        cout << area << endl;
    }

    // Method 2: parallel computation using lambda functions.
    {
        int hits = parallel_reduce(
            blocked_range<unsigned>(0, N),
            0,

            [=](const blocked_range<unsigned>& r, unsigned hits) {
                return hits + area_hits(r.begin(), r.end());
            },

            [](unsigned n1, unsigned n2) {
                return n1 + n2;
            }
        );
        int area = round(area_of(hits, N));
        cout << area << endl;
    }

    // Method 3: use class.
    {
        AreaComputer ac(N);
        parallel_reduce(blocked_range<unsigned>(0, N), ac);
        int area = round(ac.result());
        cout << area << endl;
    }

    // That's all folks
    return 0;
}
