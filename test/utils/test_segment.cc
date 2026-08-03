#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>
#include <cmath>

#include <utils/segment.hh>

using namespace fb::utils;
using Catch::Approx;

TEST_CASE("Segment default-constructs as UNKNOWN", "[utils][Segment]")
{
  Segment seg;
  CHECK(seg.type() == SegmentType::UNKNOWN);
  CHECK(seg.L() == Approx(0.0));
}

TEST_CASE("Segment kinematics follow constant-acceleration equations", "[utils][Segment]")
{
  const real L = 10.0;
  const real v0 = 5.0;
  const real a = 1.0;

  Segment seg(L, v0, /*k0=*/0.1, /*k1=*/0.2);
  seg.set_a(a);
  seg.set_type(SegmentType::FORWARD);

  const real expected_vf = std::sqrt(2.0 * a * L + v0 * v0);

  CHECK(seg.s0() == Approx(0.0));
  CHECK(seg.s1() == Approx(L));
  CHECK(seg.VF() == Approx(expected_vf));
  CHECK(seg.V(0.0) == Approx(v0));
  CHECK(seg.V(L) == Approx(expected_vf));
  CHECK(seg.type() == SegmentType::FORWARD);
}

TEST_CASE("Segment curvature interpolates linearly between k0 and k1", "[utils][Segment]")
{
  Segment seg(/*L=*/10.0, /*v0=*/5.0, /*k0=*/0.1, /*k1=*/0.3);

  CHECK(seg.K(0.0) == Approx(0.1));
  CHECK(seg.K(10.0) == Approx(0.3));
  CHECK(seg.K(5.0) == Approx(0.2));
}

TEST_CASE("Segment lateral acceleration is curvature times V^2", "[utils][Segment]")
{
  Segment seg(/*L=*/10.0, /*v0=*/5.0, /*k0=*/0.1, /*k1=*/0.1);
  seg.set_a(0.0); // constant speed -> V(s) == v0 everywhere

  CHECK(seg.AY(0.0) == Approx(0.1 * 5.0 * 5.0));
  CHECK(seg.AX(3.0) == Approx(0.0));
}

TEST_CASE("Segment::set_times derives T0/T1 from the segment's own duration", "[utils][Segment]")
{
  Segment seg(/*L=*/10.0, /*v0=*/5.0, /*k0=*/0.0, /*k1=*/0.0);
  seg.set_a(1.0);

  const real T0 = 100.0;
  seg.set_times(T0);

  CHECK(seg.getT0() == Approx(T0));
  CHECK(seg.getT1() == Approx(T0 + seg.getT()));
  CHECK(seg.getT() == Approx(seg.T()));
}
