#include <gtest/gtest.h>

#include "general/array.hpp"

using namespace beam;

TEST(ArrayTest, DefaultConstructor)
{
    Array<int> arr;
    EXPECT_EQ(arr.Size(), 0);
}

TEST(ArrayTest, SizeConstructor)
{
    Array<int> arr(5);
    EXPECT_EQ(arr.Size(), 5);
    for (int i = 0; i < 5; ++i) {
        arr[i] = i;
        EXPECT_EQ(arr[i], i);
    }
}

TEST(ArrayTest, CopyConstructor)
{
    Array<int> arr1(3);
    arr1[0] = 1; arr1[1] = 2; arr1[2] = 3;
    Array<int> arr2 = arr1;
    EXPECT_EQ(arr2.Size(), 3);
    for (int i = 0; i < 3; ++i) {
        EXPECT_EQ(arr2[i], arr1[i]);
    }
}

TEST(ArrayTest, AssignmentOperator)
{
    Array<int> arr1(2);
    arr1[0] = 10; arr1[1] = 20;
    Array<int> arr2;
    arr2 = arr1;
    EXPECT_EQ(arr2.Size(), 2);
    EXPECT_EQ(arr2[0], 10);
    EXPECT_EQ(arr2[1], 20);
}

TEST(ArrayTest, ElementAccess)
{
    Array<int> arr(4);
    for (int i = 0; i < 4; ++i) arr[i] = i * 2;
    EXPECT_EQ(arr[0], 0);
    EXPECT_EQ(arr[1], 2);
    EXPECT_EQ(arr[2], 4);
    EXPECT_EQ(arr[3], 6);
}

TEST(ArrayTest, ConstElementAccess)
{
    Array<int> arr(2);
    arr[0] = 7; arr[1] = 8;
    const Array<int>& carr = arr;
    EXPECT_EQ(carr[0], 7);
    EXPECT_EQ(carr[1], 8);
}

// Add more tests if Array supports push_back, pop_back, clear, resize, etc.