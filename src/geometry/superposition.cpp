/**
 * @file superposition.cpp
 * @brief Horn (1987) quaternion least-squares superposition.
 */
#include <pairfinder/geometry/superposition.hpp>

#include <cmath>

namespace pairfinder::geometry {

namespace {

// Symmetric 4x4 Jacobi eigensolver. On return, eigenvectors are columns of V
// and eigenvalues in d. (Cyclic Jacobi; ample iterations for a 4x4.)
void jacobi4(double a[4][4], double d[4], double v[4][4]) {
    for (int i = 0; i < 4; ++i) {
        for (int j = 0; j < 4; ++j) v[i][j] = (i == j) ? 1.0 : 0.0;
    }
    for (int sweep = 0; sweep < 100; ++sweep) {
        double off = 0.0;
        for (int p = 0; p < 4; ++p)
            for (int q = p + 1; q < 4; ++q) off += a[p][q] * a[p][q];
        if (off < 1e-30) break;
        for (int p = 0; p < 4; ++p) {
            for (int q = p + 1; q < 4; ++q) {
                if (std::fabs(a[p][q]) < 1e-300) continue;
                const double theta = (a[q][q] - a[p][p]) / (2.0 * a[p][q]);
                double t = (theta >= 0 ? 1.0 : -1.0) /
                           (std::fabs(theta) + std::sqrt(theta * theta + 1.0));
                const double c = 1.0 / std::sqrt(t * t + 1.0);
                const double s = t * c;
                for (int i = 0; i < 4; ++i) {
                    const double aip = a[i][p], aiq = a[i][q];
                    a[i][p] = c * aip - s * aiq;
                    a[i][q] = s * aip + c * aiq;
                }
                for (int i = 0; i < 4; ++i) {
                    const double api = a[p][i], aqi = a[q][i];
                    a[p][i] = c * api - s * aqi;
                    a[q][i] = s * api + c * aqi;
                }
                for (int i = 0; i < 4; ++i) {
                    const double vip = v[i][p], viq = v[i][q];
                    v[i][p] = c * vip - s * viq;
                    v[i][q] = s * vip + c * viq;
                }
            }
        }
    }
    for (int i = 0; i < 4; ++i) d[i] = a[i][i];
}

}  // namespace

SuperpositionResult superpose(const std::vector<Vector3d>& from,
                              const std::vector<Vector3d>& to) {
    SuperpositionResult out;
    const std::size_t n = from.size();
    if (n < 3 || to.size() != n) return out;

    Vector3d cf, ct;
    for (std::size_t i = 0; i < n; ++i) { cf = cf + from[i]; ct = ct + to[i]; }
    cf = cf / static_cast<double>(n);
    ct = ct / static_cast<double>(n);

    // Cross-covariance S_ab = sum (from_centered)_a * (to_centered)_b
    double Sxx = 0, Sxy = 0, Sxz = 0, Syx = 0, Syy = 0, Syz = 0, Szx = 0, Szy = 0, Szz = 0;
    for (std::size_t i = 0; i < n; ++i) {
        const Vector3d f = from[i] - cf;
        const Vector3d t = to[i] - ct;
        Sxx += f.x * t.x; Sxy += f.x * t.y; Sxz += f.x * t.z;
        Syx += f.y * t.x; Syy += f.y * t.y; Syz += f.y * t.z;
        Szx += f.z * t.x; Szy += f.z * t.y; Szz += f.z * t.z;
    }

    // Horn's 4x4 key matrix N (symmetric).
    double N[4][4] = {
        {Sxx + Syy + Szz, Syz - Szy,       Szx - Sxz,        Sxy - Syx},
        {Syz - Szy,       Sxx - Syy - Szz, Sxy + Syx,        Szx + Sxz},
        {Szx - Sxz,       Sxy + Syx,       -Sxx + Syy - Szz, Syz + Szy},
        {Sxy - Syx,       Szx + Sxz,       Syz + Szy,        -Sxx - Syy + Szz},
    };

    double d[4], V[4][4];
    jacobi4(N, d, V);
    int best = 0;
    for (int i = 1; i < 4; ++i) if (d[i] > d[best]) best = i;
    double q0 = V[0][best], q1 = V[1][best], q2 = V[2][best], q3 = V[3][best];
    const double qn = std::sqrt(q0 * q0 + q1 * q1 + q2 * q2 + q3 * q3);
    if (qn < 1e-12) return out;
    q0 /= qn; q1 /= qn; q2 /= qn; q3 /= qn;

    // Rotation from unit quaternion (w=q0, x=q1, y=q2, z=q3); rotates from->to.
    auto& R = out.rotation;
    R[0][0] = 1 - 2 * (q2 * q2 + q3 * q3);
    R[0][1] = 2 * (q1 * q2 - q0 * q3);
    R[0][2] = 2 * (q1 * q3 + q0 * q2);
    R[1][0] = 2 * (q1 * q2 + q0 * q3);
    R[1][1] = 1 - 2 * (q1 * q1 + q3 * q3);
    R[1][2] = 2 * (q2 * q3 - q0 * q1);
    R[2][0] = 2 * (q1 * q3 - q0 * q2);
    R[2][1] = 2 * (q2 * q3 + q0 * q1);
    R[2][2] = 1 - 2 * (q1 * q1 + q2 * q2);

    auto apply = [&R](const Vector3d& p) {
        return Vector3d{R[0][0] * p.x + R[0][1] * p.y + R[0][2] * p.z,
                        R[1][0] * p.x + R[1][1] * p.y + R[1][2] * p.z,
                        R[2][0] * p.x + R[2][1] * p.y + R[2][2] * p.z};
    };
    out.translation = ct - apply(cf);

    double se = 0.0;
    for (std::size_t i = 0; i < n; ++i) {
        const Vector3d mapped = apply(from[i]) + out.translation;
        se += (mapped - to[i]).dot(mapped - to[i]);
    }
    out.rms = std::sqrt(se / static_cast<double>(n));
    out.valid = true;
    return out;
}

}  // namespace pairfinder::geometry
