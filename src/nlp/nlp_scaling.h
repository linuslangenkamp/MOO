#ifndef OPT_NLP_SCALING_H
#define OPT_NLP_SCALING_H

#include "base/fixed_vector.h"
#include "base/util.h"

/**
 * @brief Interface for NLP variable and function scaling.
 *
 * The Scaling interface defines how to scale and unscale the components of a nonlinear program (NLP).
 * This is useful for improving numerical conditioning and solver performance.
 *
 * All scaling operations should be implemented in-place and must be invertible (i.e. unscale must exactly
 * reverse the corresponding scale operation).
 *
 * This interface is entirely user-defined — the NLP implementation itself does not automatically apply
 * any scaling or unscaling. It is up to the user or solver implementation to call the appropriate methods
 * at the correct points during the solve.
 */
class Scaling {
public:
    virtual void scale_x(FixedVector<f64>& x) const = 0;
    virtual void unscale_x(FixedVector<f64>& x) const = 0;

    virtual void scale_f(f64& f) const = 0;
    virtual void unscale_f(f64& f) const = 0;

    virtual void scale_g(FixedVector<f64>& g) const = 0;
    virtual void unscale_g(FixedVector<f64>& g) const = 0;

    virtual void scale_lambda(FixedVector<f64>& lambda) const = 0;
    virtual void unscale_lambda(FixedVector<f64>& lambda) const = 0;

    virtual void scale_sigma(f64&) const = 0;
    virtual void unscale_sigma(f64&) const = 0;

    virtual ~Scaling() = default;
};

/**
 * @brief No-op scaling: identity transformations (set as default)
 */
class NoScaling : public Scaling {
public:
    void scale_x(FixedVector<f64>&) const override {}
    void unscale_x(FixedVector<f64>&) const override {}

    void scale_f(f64&) const override {}
    void unscale_f(f64&) const override {}

    void scale_g(FixedVector<f64>&) const override {}
    void unscale_g(FixedVector<f64>&) const override {}

    void scale_lambda(FixedVector<f64>&) const override {}
    void unscale_lambda(FixedVector<f64>&) const override {}

    void scale_sigma(f64&) const override {}
    void unscale_sigma(f64&) const override {}
};

#endif // OPT_NLP_SCALING_H
