#include "elff/general/random.hpp"
#include <gtest/gtest.h>
#include <limits>

using ELFF::CounterRandom;
namespace {
// Independent transcription of the supplied ParticleRandom transform.
std::uint64_t reference_mix(std::uint64_t x) {
  x += 0x9e3779b97f4a7c15ULL;
  x = (x ^ (x >> 30)) * 0xbf58476d1ce4e5b9ULL;
  x = (x ^ (x >> 27)) * 0x94d049bb133111ebULL;
  return x ^ (x >> 31);
}
std::uint64_t reference_fold(std::uint64_t a, std::uint64_t b) {
  return reference_mix(a ^ reference_mix(b));
}
double reference_uniform(std::uint64_t seed, std::uint64_t stream, std::uint64_t i) {
  const auto r=reference_fold(reference_fold(seed,stream),i);
  return double(r>>11)*0x1.0p-53;
}
double reference_normal(std::uint64_t seed, std::uint64_t stream, std::uint64_t i) {
  double u=reference_uniform(seed,stream,2*i);
  if(u<=0) u=0x1.0p-53;
  return std::sqrt(-2.*std::log(u))*std::cos(2.*3.14159265358979323846*reference_uniform(seed,stream,2*i+1));
}
TEST(CounterRandomTest, ExactIntegerVectorsAndZeroUniform) {
  const CounterRandom random(1,0x3333333333333333ULL);
  const std::uint64_t expected[]={0x4f20f7ec5cf695e7ULL,0x6da6be9294ad1f2fULL,
    0x7aa471b9e3125c14ULL,0x112bb94c252d3108ULL};
  for(unsigned i=0;i<4;++i) EXPECT_EQ(random.u64(i),expected[i]);
  EXPECT_EQ(CounterRandom(0,0).u64(0),0xa7f72697a2731486ULL);
  const auto max=std::numeric_limits<std::uint64_t>::max();
  EXPECT_EQ(CounterRandom(max,max).u64(max),0xcc56cdcdc256e448ULL);
  EXPECT_EQ(random.u64(0xbe7948e16e9eca05ULL),0u);
  EXPECT_EQ(random.uniform01(0xbe7948e16e9eca05ULL),0.);
  EXPECT_EQ(random.uniform01_open(0xbe7948e16e9eca05ULL),0x1.0p-53);
}
TEST(CounterRandomTest, PreservesOriginalSamplingAndCounterOverflow) {
  const std::uint64_t values[]={0,1,17,0xffffffffffffffffULL};
  for(auto seed:values) for(auto stream:values) for(auto index:values) {
    const CounterRandom random(seed,stream);
    EXPECT_EQ(random.uniform01(index),reference_uniform(seed,stream,index));
    EXPECT_DOUBLE_EQ(random.normal(index),reference_normal(seed,stream,index));
    for(unsigned dimensions=1;dimensions<=3;++dimensions) {
      const auto g=random.gaussian_vector(index,dimensions);
      const auto p=random.domain_position(index,{-2,3,5},{7,8,9},dimensions);
      for(unsigned d=0;d<3;++d) {
        EXPECT_DOUBLE_EQ(g[d],d<dimensions ? reference_normal(seed,stream,3*index+d) : 0.);
        const double origin[]={-2,3,5}, extent[]={7,8,9};
        EXPECT_DOUBLE_EQ(p[d],d<dimensions ? origin[d]+extent[d]*reference_uniform(seed,stream,3*index+d) : 0.);
      }
    }
  }
}
}
