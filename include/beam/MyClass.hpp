// include/mylib/MyClass.hpp
#pragma once

namespace beam {

/**
 * @brief A simple demonstration class
 * 
 * MyClass provides a basic example of class documentation and functionality.
 * It stores a base value and can perform computations using this value.
 * 
 *   The distance between \f$(x_1,y_1)\f$ and \f$(x_2,y_2)\f$ is 
 * \f(\sqrt{(x_2-x_1)^2+(y_2-y_1)^2}\f)
 * 
 *   \f[
    |I_2|=\left| \int_{0}^T \psi(t) 
             \left\{ 
                u(a,t)-
                \int_{\gamma(t)}^a 
                \frac{d\theta}{k(\theta,t)}
                \int_{a}^\theta c(\xi)u_t(\xi,t)\,d\xi
             \right\} dt
          \right|
  \f].
 */
class MyClass {
public:
    /**
     * @brief Construct a new MyClass object
     * @param init_value The initial base value to store
     */
    explicit MyClass(int init_value);
    
    /**
     * @brief Destroy the MyClass object
     */
    ~MyClass();

    /**
     * @brief Compute a result based on input and stored base value
     * @param x The input value to process
     * @return int The computed result (base + x²)
     */
    int compute(int x) const;

private:
    int _base;  ///< The stored base value
};

} // namespace beam
