#ifndef OPT_FIXED_VECTOR_H
#define OPT_FIXED_VECTOR_H

#include <cstdint>
#include <memory>
#include <algorithm>
#include <functional>

/**
 * @brief Templated, zero-initialized vector with fixed size.
 *
 * This class represents a fixed-size vector that is zero-initialized upon construction.
 * It provides basic functionality for accessing and manipulating elements, as well as
 * support for copying, nulling and moving vectors.
 *
 * @tparam T The type of elements stored in the vector.
 */
template<typename T>
class FixedVector {
public:

    constexpr FixedVector() noexcept : _size{0}, _data{} {}

    explicit FixedVector(std::size_t size) : _size{size}, _data{std::make_unique<T[]>(size)} {}

    FixedVector(std::size_t size, const std::function<T()> & generator) : FixedVector(size) {
        std::generate(_data.get(), _data.get() + _size, generator);
    }

    constexpr FixedVector(const FixedVector &other) : FixedVector(other._size) {
        *this = other;
    }
    /*
    template<typename... Args>
    explicit FixedVector(std::size_t first, Args... rest) : FixedVector(first) {
        _init_variadic(rest...);
    }
    */
    FixedVector(FixedVector &&other) noexcept = default;

    // assign based on iterator begin() and end()
    template<typename It>
    constexpr FixedVector(It first, It last) : FixedVector(last - first) {
        assign(first, last);
    }

    // access vector at index 0 <= index < vec.size()
    constexpr T& operator[](std::size_t index) {
        assert(index < _size);

        return _data[index];
    }

    // access vector at index 0 <= index < vec.size()
    constexpr const T& operator[](std::size_t index) const {
        assert(index < _size);

        return _data[index];
    }
    
    // assign vector to other vector of equal size
    constexpr FixedVector& operator=(const FixedVector &other) {
        assert(_size == other._size);

        memcpy(_data.get(), other._data.get(), _size * sizeof(T));

        return *this;
    }

    FixedVector& operator=(FixedVector&& other) noexcept = default;

    // fill entire vector with 0
    constexpr void fill_zero() {
        memset(_data.get(), 0, _size * sizeof(T));
    }

    constexpr void assign(const T* data, std::size_t len, std::size_t offset = 0) {
        assert(data != nullptr);
        assert(len <= _size - offset);

        memcpy(_data.get() + offset, data, len * sizeof(T));
    }

    template<typename It>
    constexpr void assign(It first, It last, std::size_t offset = 0) {
        assert(last - first <= _size - offset);

        std::copy(first, last, _data.get() + offset);
    }

    constexpr std::size_t size() const {
        return _size;
    }

    constexpr T& back() {
        assert(_size != 0);

        return _data[_size - 1];
    }

    constexpr const T& back() const {
        assert(_size != 0);

        return _data[_size - 1];
    }

    constexpr T* raw() {
        return _data.get();
    }

    constexpr const T* raw() const {
        return _data.get();
    }

    constexpr T* begin() noexcept {
        return _data.get();
    }

    constexpr T* end() noexcept {
        return _data.get() + _size;
    }

    constexpr const T* begin() const noexcept {
        return _data.get();
    }

    constexpr const T* end() const noexcept {
        return _data.get() + _size;
    }

    constexpr const bool empty() const noexcept {
        return (_size == 0);
    }

private:
    std::size_t _size;
    std::unique_ptr<T[]> _data;

    template<typename... Args>
    void _init_variadic(std::size_t first, Args... rest) {
        for (std::size_t i = 0; i < _size; i++) {
            _data[i] = FixedVector<T>(first, rest...);
        }
    }

    void _init_variadic(std::size_t first) {
        for (std::size_t i = 0; i < _size; i++) {
            _data[i] = FixedVector<T>(first);
        }
    }
};

template<typename T, std::size_t Dim>
struct FixedMultiVectorRecursive {
    using type = FixedVector<typename FixedMultiVectorRecursive<T, Dim - 1>::type>;
};

template<typename T>
struct FixedMultiVectorRecursive<T, 1> {
    using type = FixedVector<T>;
};

template<typename T, std::size_t Dim>
using FixedMultiVector = typename FixedMultiVectorRecursive<T, Dim>::type;

#endif // OPT_FIXED_VECTOR_H
