// shapes.hpp
#pragma once
#include <array>
#include <cstddef>

/// Cubic–Hermite shape functions on [0,1] (with element length h)
template <class T>
struct CubicHermite
{
    /// nodal displacement shape
    
    static T H1(T xi) noexcept { return 1.0 - 3.0*xi*xi + 2.0*xi*xi*xi; }
    static T H3(T xi) noexcept { return 3.0*xi*xi - 2.0*xi*xi*xi; }

    /// nodal slope (rotation) shape
    static T H2(T xi, T h) noexcept { return h*(xi - 2.0*xi*xi + xi*xi*xi); }
    static T H4(T xi, T h) noexcept { return h*(xi*xi*xi - xi*xi); }

    /// first derivatives w.r.t. reference ξ
    static T dH1_dxi(T xi) noexcept { return -6.0*xi + 6.0*xi*xi; }
    static T dH2_dxi(T xi, T h) noexcept { return h*(1.0 - 4.0*xi + 3.0*xi*xi); }
    static T dH3_dxi(T xi) noexcept { return  6.0*xi - 6.0*xi*xi; }
    static T dH4_dxi(T xi, T h) noexcept { return h*(3.0*xi*xi - 2.0*xi); }

    /// second derivatives w.r.t. reference ξ
    static T d2H1_dxi2(T xi) noexcept { return -6.0 + 12.0*xi; }
    static T d2H2_dxi2(T xi, T h) noexcept { return h*(-4.0 + 6.0*xi); }
    static T d2H3_dxi2(T xi) noexcept { return  6.0 - 12.0*xi; }
    static T d2H4_dxi2(T xi, T h) noexcept { return h*(6.0*xi - 2.0); }

    /// values [H1,H2,H3,H4]
    static std::array<T,4> values(T xi, T h) noexcept {
        return { H1(xi), H2(xi,h), H3(xi), H4(xi,h) };
    }

    /// first derivatives w.r.t. physical x (i.e. d/dx = (1/h) d/dξ)
    static std::array<T,4> derivs(T xi, T h) noexcept {
        return { dH1_dxi(xi)/h,
                 dH2_dxi(xi,h)/h,
                 dH3_dxi(xi)/h,
                 dH4_dxi(xi,h)/h };
    }

    /// second derivatives w.r.t. physical x (i.e. d²/dx² = (1/h²) d²/dξ²)
    static std::array<T,4> second_derivs(T xi, T h) noexcept {
        return { d2H1_dxi2(xi)/(h*h),
                 d2H2_dxi2(xi,h)/(h*h),
                 d2H3_dxi2(xi)/(h*h),
                 d2H4_dxi2(xi,h)/(h*h) };
    }
};

/// Simple linear (hat) shape on [0,1]
template <class T>
struct LinearShape
{
    static T M1(T xi) noexcept { return 1.0 - xi; }
    static T M2(T xi) noexcept { return xi; }

    static std::array<T,2> values(T xi) noexcept {
        return { M1(xi), M2(xi) };
    }
};
