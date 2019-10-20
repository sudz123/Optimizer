#include <gtest/gtest.h>
#include <Eigen/Dense>
#include "../../Optimizer/optimizer"

#include <string>

namespace util
{
double plus_ten_square(double x)
{
    return std::pow(x + 10, 2);
}

template <typename T>
auto transform_precision(T in, int precision) -> T
{
    std::stringstream ss;
    ss << std::setprecision(precision) << in;
    return std::stod(ss.str());
}
} // namespace util

TEST(SingleVariable, BoundingPhase)
{
    double initial_point = 5.4;
    auto result = Optimizer::BoundingPhase(util::plus_ten_square, initial_point);
    EXPECT_EQ(util::transform_precision(result(0), 5), util::transform_precision(-26.1, 5));
    EXPECT_EQ(util::transform_precision(result(1), 5), util::transform_precision( -2.1, 5));
    EXPECT_EQ(result.size(), 2);
}

TEST(SingleVariable, ExhaustiveSearch)
{
    double initial_point = 5.4;
    auto result = Optimizer::Exhaustive(util::plus_ten_square, initial_point);
    EXPECT_EQ(util::transform_precision(result(0), 5), util::transform_precision(-10.6, 5));
    EXPECT_EQ(util::transform_precision(result(1), 5), util::transform_precision( -9.6, 5));
    EXPECT_EQ(result.size(), 2);
}

TEST(SingleVariable, Derivatives)
{
    double initial_point = 5.4;
    auto result = Optimizer::Derivative(util::plus_ten_square, initial_point);
    EXPECT_EQ(util::transform_precision(result(0), 5), util::transform_precision(30.8,     5));
    EXPECT_EQ(util::transform_precision(result(1), 6), util::transform_precision( 1.99975, 6));
    EXPECT_EQ(result.size(), 2);
}

int main(int argc, char* argv[])
{
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
