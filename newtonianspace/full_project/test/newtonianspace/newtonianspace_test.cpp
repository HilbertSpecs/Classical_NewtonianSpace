#include "newtonianspace/NewtonianSpace.hpp"
#include <gtest/gtest.h>
#include <gmock/gmock.h>

class MockNewtonianSpace : public NewtonianSpace {
public:
    MOCK_METHOD(void, setMass, (double mass), (override));
    MOCK_METHOD(double, getMass, (), (const, override));
};

//using ::testing::AtLeast;
using ::testing::Return;
using ::testing::AtLeast;

TEST(NewtonianSpaceTest, SetMass) {
    MockNewtonianSpace newtonspace;
    EXPECT_CALL(newtonspace, setMass(45)).Times(AtLeast(1));
    EXPECT_CALL(newtonspace, getMass()).Times(3).WillRepeatedly(Return(45));
    //EXPECT_EQ(45==testMass.getMass(), testMass.setMass(45));
}

/*
TEST(NewtonianSpaceTest, MultiplyMultipliesTwoInts) {
  EXPECT_EQ(12, Calc::Multiply(3, 4));
}
*/
