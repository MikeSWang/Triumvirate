#include <vector>
#include <complex>

#include <gtest/gtest.h>

#include "fftlog.hpp"

// Test suite: HankelTransformTest

// Test fixture
class HankelTransformTest : public ::testing::Test {
 protected:
  void SetUp() override {
    // Create an instance with fixed parameters.
    const double MU = .5;
    const double Q = .1;
    const bool THREADED = true;

    this->ht = new trv::maths::HankelTransform(MU, Q, THREADED);

    // Initialise the instance.
    this->sample_pts = {.1, .2, .4, .8, 1.6, 3.2, 6.4, 12.8};
    this->kr_c = 1.;
    this->lowring = true;
    this->extrap = trva::ExtrapOption::NONE;
    this->extrap_exp = 2.;

    this->ht->initialise(
      this->sample_pts, this->kr_c, this->lowring,
      this->extrap, this->extrap_exp
    );
  }

  void TearDown() override {
    // Delete the instance.
    if (this->ht != nullptr) {
      delete this->ht; this->ht = nullptr;
    }
  }

  // Test data members
  trv::maths::HankelTransform* ht = nullptr;
  std::vector<double> sample_pts;
  double kr_c;
  bool lowring;
  trva::ExtrapOption extrap;
  double extrap_exp;
};

// Test method: test_initialise
TEST_F(HankelTransformTest, test_initialise) {
  // Call the initialise() function.
  ASSERT_NO_THROW(
    ht->initialise(sample_pts, kr_c, lowring, extrap, extrap_exp)
  );
}

// Test method: test_biased_transform
TEST_F(HankelTransformTest, test_biased_transform) {
  // Define input and output arrays.
  std::complex<double> a[8] = {1., 2., 3., 4., 5., 6., 7., 8.};
  std::complex<double> b[8];

  // Call the biased_transform() function.
  ASSERT_NO_THROW(
    ht->biased_transform(a, b);
  );

  // Check the output array values
  EXPECT_NEAR(b[0].real(), 2.492884552, 1.e-9);
  EXPECT_NEAR(b[7].real(), 2.080611583, 1.e-9);
}

// Test run: all
int main(int argc, char** argv) {
  testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
