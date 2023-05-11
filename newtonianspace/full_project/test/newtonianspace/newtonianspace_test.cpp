#include "newtonianspace/NewtonianSpace.hpp"

#include <gtest/gtest.h>

TEST(NewtonianSpaceTest, SetMass) {
  EXPECT_EQ(45==NewtonianSpace::getMass(), NewtonianSpace::setMass(45));
}

/*
TEST(NewtonianSpaceTest, MultiplyMultipliesTwoInts) {
  EXPECT_EQ(12, Calc::Multiply(3, 4));
}
*/
