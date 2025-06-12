#include <gtest/gtest.h>
#include <beam/LinAlg/Vector.hpp>

using beam::Vector;

// Test default‐size construction, element assignment & access, and size()
TEST(VectorTest, SizeAndElementAccess) {
    Vector<int> v(3);
    EXPECT_EQ(v.size(), 3u);

    // assign and read back
    v(0) = 10;
    v(1) = 20;
    v(2) = 30;
    EXPECT_EQ(v(0), 10);
    EXPECT_EQ(v(1), 20);
    EXPECT_EQ(v(2), 30);
}

// Test initializer_list constructor and const access
TEST(VectorTest, InitializerListAndConstAccess) {
    Vector<double> w{1.1, 2.2, 3.3};
    EXPECT_EQ(w.size(), 3u);

    // non-const access
    EXPECT_DOUBLE_EQ(w(0), 1.1);
    EXPECT_DOUBLE_EQ(w(1), 2.2);
    EXPECT_DOUBLE_EQ(w(2), 3.3);

    // const reference access
    const Vector<double>& cw = w;
    EXPECT_DOUBLE_EQ(cw(0), 1.1);
    EXPECT_DOUBLE_EQ(cw(1), 2.2);
    EXPECT_DOUBLE_EQ(cw(2), 3.3);
}

// Test out‐of‐bounds indexing throws std::out_of_range
TEST(VectorTest, OutOfRangeAccess) {
    Vector<int> v(2);
    EXPECT_THROW(v(2), std::out_of_range);
    EXPECT_THROW(v( 100 ), std::out_of_range);

    const Vector<int> cv(1);
    EXPECT_THROW(cv(1), std::out_of_range);
}

// Test zero-length vector
TEST(VectorTest, ZeroLength) {
    Vector<int> empty(0);
    EXPECT_EQ(empty.size(), 0u);
    EXPECT_THROW(empty(0), std::out_of_range);
}

// main() is provided by GTest’s gtest_main library
