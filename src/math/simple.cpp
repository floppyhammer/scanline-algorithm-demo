#include "simple.h"

#include <cmath>

float my_round(float var, int decimal) {
    // 37.66666 * 100 =3766.66
    // 3766.66 + .5 =3767.16    for rounding off value
    // then type cast to int so value is 3767
    // then divided by 100 so the value converted into 37.67
    float factor = pow(10.0, decimal);
    float value = (int)(var * factor + 0.5);
    return (float)value / factor;
}

void cubic_roots(const float* polynomial, float* root) {
    // cubic function: f(x)=ax^3+bx^2+cx+d (a ≠ 0)
    float a = polynomial[0];
    float b = polynomial[1];
    float c = polynomial[2];
    float d = polynomial[3];

    //printf("\nPolynomial %f %f %f %f", a, b, c, d);

    if (a == 0) {
        if (b == 0) {
            if (c == 0) {
                root[0] = -1.0f;
                root[1] = -1.0f;
                root[2] = -1.0f;
                return;
            }
            root[0] = -1.0f * (d / c);
            root[1] = -1.0f;
            root[2] = -1.0f;

            // discard out of spec roots
            if (root[0] == 1.0f)
                root[0] = -1.0f;

            // sort but place -1 at the end
            sort(root, 3);
            return;
        }

        // quadratic discriminant
        float DQ = pow(c, 2.0) - 4.0 * b * d;
        if (DQ >= 0) {
            DQ = sqrt(DQ);

            root[0] = -1.0f * ((DQ + c) / (2.0f * b));
            root[1] = ((DQ - c) / (2.0f * b));
            root[2] = -1.0f;

            // discard out of spec roots
            for (int i = 0; i < 2; i++) {
                if (root[i] == 1.0)
                    root[i] = -1.0;
            }
        }

        // sort but place -1 at the end
        sort(root, 3);

        return;
    }

    float A = b / a;
    float B = c / a;
    float C = d / a;

    float Q = (3.0f * B - pow(A, 2.0f)) / 9.0f;
    float R = (9.0f * A * B - 27.0f * C - 2.0f * pow(A, 3.0f)) / 54.0f;

    // polynomial discriminant
    float D = pow(Q, 3.0f) + pow(R, 2.0f);
    //printf("\nD %f", D);
    // complex or duplicate roots
    if (D >= 0.0f) {
        float S = sign(R + sqrt(D)) * pow(abs(R + sqrt(D)), (1.0f / 3.0f));
        float T = sign(R - sqrt(D)) * pow(abs(R - sqrt(D)), (1.0f / 3.0f));

        root[0] = -A / 3.0f + (S + T); // real root
        root[1] = -A / 3.0f - (S + T) / 2.0f; // real part of complex root
        root[2] = -A / 3.0f - (S + T) / 2.0f; // real part of complex root

        float imaginary = abs(sqrt(3.0f) * (S - T) / 2.0f); // complex part of root pair

        // discard complex roots
        if (imaginary != 0.0f) {
            root[1] = -1.0f;
            root[2] = -1.0f;
        }
    }
    // distinct real roots
    else {
        float th = acos(R / sqrt(-pow(Q, 3.0f)));

        root[0] = 2.0f * sqrt(-Q) * cos(th / 3.0f) - A / 3.0f;
        root[1] = 2.0f * sqrt(-Q) * cos((th + 2.0f * PI) / 3.0f) - A / 3.0f;
        root[2] = 2.0f * sqrt(-Q) * cos((th + 4.0f * PI) / 3.0f) - A / 3.0f;
    }

    // discard out of roots exceding the curve range (0, 1)
    for (int i = 0; i < 3; i++) {
        if (root[i] < 0.0f || root[i] > 1.0f) root[i] = -1.0f;
    }

    // sort but place -1 at the end
    sort(root, 3);
}

void bezier_coeffs(float P0, float P1, float P2, float P3, float* Z) {
    Z[0] = -P0 + 3.0f * P1 - 3.0f * P2 + P3;
    Z[1] = 3.0f * P0 - 6.0f * P1 + 3.0f * P2;
    Z[2] = -3.0f * P0 + 3.0f * P1;
    Z[3] = P0;
}

/**
 * Compute intersection points between the scanline and a Bezier curve.
 *
 * @param p Bezier points.
 * @param lx Scanline x.
 * @param ly Scanline y.
 * @param I Intersection points.
 * @return Intersection points.
 */
void compute_intersections(const float* p, const float* lx, const float* ly, float* I) {
    // A=y2-y1
    float A = ly[1] - ly[0];
    // B=x1-x2
    float B = lx[0] - lx[1];
    // C=x1*(y1-y2)+y1*(x2-x1)
    float C = lx[0] * (ly[0] - ly[1]) + ly[0] * (lx[1] - lx[0]);

    //printf("ABC %f %f %f", A, B, C);

    float bx[] = {0.0, 0.0, 0.0, 0.0};
    float by[] = {0.0, 0.0, 0.0, 0.0};

//    printf("\nP");
//    for (int i = 0; i < 8; i++) {
//        printf(" %f", p[i]);
//    }

    bezier_coeffs(p[0], p[2], p[4], p[6], bx);
    bezier_coeffs(p[1], p[3], p[5], p[7], by);

//    printf("\nbx");
//    for (int i = 0; i < 4; i++) {
//        printf(" %f", bx[i]);
//    }
//    printf("\nby");
//    for (int i = 0; i < 4; i++) {
//        printf(" %f", by[i]);
//    }
//
//    printf("\nZ1");
//    for (float i : bx) {
//        printf(" %f", i);
//    }
//    printf("\nZ2");
//    for (float i : by) {
//        printf(" %f", i);
//    }

    float P[] = {0.0, 0.0, 0.0, 0.0};
    P[0] = A * bx[0] + B * by[0]; // t^3
    P[1] = A * bx[1] + B * by[1]; // t^2
    P[2] = A * bx[2] + B * by[2]; // t
    P[3] = A * bx[3] + B * by[3] + C; // 1

    float r[] = {-1.0, -1.0, -1.0};
    cubic_roots(P, r);

//    for (int i = 0; i < 3; i++) {
//        //if (r[i] > -1) {
//            printf("\nroot %d %f", i, r[i]);
//        //}
//    }

    // Verify the roots are in bounds of the linear segment.
    for (int i = 0; i < 3; i++) {
        // Linear coordinate
        float t = r[i];

        float X[2];
        X[0] = bx[0] * pow(t, 3.0f) + bx[1] * pow(t, 2.0f) + bx[2] * t + bx[3];
        X[1] = by[0] * pow(t, 3.0f) + by[1] * pow(t, 2.0f) + by[2] * t + by[3];

        // above is intersection point assuming infinitely long line segment,
        // make sure we are also in bounds of the line
        float s;

        // if not vertical line
        if ((lx[1]-lx[0]) != 0) {
            s = (X[0]-lx[0])/(lx[1]-lx[0]);
        } else {
            s = (X[1]-ly[0])/(ly[1]-ly[0]);
        }

        // in bounds?
        if (t < 0.0f || t > 1.0f || s < 0.0f || s > 1.0f) {
            // move off screen
            X[0] = -10000.0f;
            X[1] = -10000.0f;
        }

        // move intersection point
        I[i * 2] = X[0];
        I[i * 2 + 1] = X[1];
    }
}
