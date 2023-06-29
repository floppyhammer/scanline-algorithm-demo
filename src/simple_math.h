#ifndef MATH_SIMPLE_CPP
#define MATH_SIMPLE_CPP

const float PI = 3.14159265358979323846;

template <typename T>
inline T sign(T val) {
    return (T(0) < val) - (val < T(0));
}

template <typename T>
inline void swap(T* a, T* b) {
    T temp = *a;
    *a = *b;
    *b = temp;
}

/// Perform selection sorting.
template <typename T>
inline void sort(T arr[], int n) {
    // One by one move boundary of unsorted subarray.
    for (int i = 0; i < n - 1; i++) {
        // Find the minimum element in unsorted array.
        int min_idx = i;
        for (int j = i + 1; j < n; j++) {
            if (arr[j] < arr[min_idx]) min_idx = j;
        }

        // Swap the found minimum element with the first element.
        swap(&arr[min_idx], &arr[i]);
    }
}

void cubic_roots(const float* polynomial, float* root);

void bezier_coeffs(float P0, float P1, float P2, float P3, float* Z);

/// Computes intersection between a cubic spline and a line segment.
void compute_intersections(const float* p, const float* lx, const float* ly, float* I);

#endif // MATH_SIMPLE_CPP
