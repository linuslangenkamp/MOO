#ifndef OPT_UTIL_H
#define OPT_UTIL_H

#include <functional>
#include <vector>
#include <cmath>

/* simple typedef for the Number, for using f32 or something later */
typedef double f64;

const f64 PLUS_INFINITY = std::numeric_limits<f64>::infinity();
const f64 MINUS_INFINITY = -std::numeric_limits<f64>::infinity();

template <typename T>
inline int int_size(const std::vector<T>& vec) {
    return static_cast<int>(vec.size());
}

/* AutoFree manages a list of raw pointers and their corresponding free functions.
 * When an AutoFree object goes out of scope, it automatically calls each stored free function
 * on its associated pointer to release resources and avoid memory leaks.
 * You register pointers and their free functions using the attach() method.
 * Used for pure C-style mallocs / callocs in the OPT module */
class AutoFree {
public:
    ~AutoFree() {
        for (auto& [ptr, free_fn] : _to_free) {
            free_fn(ptr);
        }
    }

    template <typename T>
    void attach(T* ptr, void (*free_fn)(T*)) {
        _to_free.emplace_back(ptr, [free_fn](void* p) { free_fn(static_cast<T*>(p)); });
    }

private:
    std::vector<std::pair<void*, std::function<void(void*)>>> _to_free;
};

#endif // OPT_UTIL_H
