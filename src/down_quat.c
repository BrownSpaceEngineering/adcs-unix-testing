#include "include/down_quat.h"
#include "include/laextension.h"
#include "include/quat.h"
#include "math.h"
#define EPS 1e-9
static void normalize3(double *v)
{
    double n = l2_norm(v, 3);
    if (n > EPS) {
        v[0] /= n;
        v[1] /= n;
        v[2] /= n;
    }
}

static double dot3(const double *a, const double *b)
{
    return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
}
// --- Main function ---

void down_quat(double* from, double* to, double* q_body_to_eci, double* goal_q)
{
    double nadir[3];

    // nadir = to - from: axis pointing from 'from' toward 'to'
    nadir[0] = to[0] - from[0];
    nadir[1] = to[1] - from[1];
    nadir[2] = to[2] - from[2];

    double nadir_norm = l2_norm(nadir, 3);
    if (nadir_norm < EPS) {
        // Degenerate case: return identity quaternion
        goal_q[0] = 1.0;
        goal_q[1] = 0.0;
        goal_q[2] = 0.0;
        goal_q[3] = 0.0;
        return;
    }
    normalize3(nadir);

    // Get current rotation matrix (Body → ECI)
    double R[9];
    quat2rotm(q_body_to_eci, R);

    // Extract body Y-axis in ECI (2nd column)
    double y[3] = { R[1], R[4], R[7] };

    // Project y onto plane perpendicular to nadir
    double proj = dot3(y, nadir);
    double new_y[3] = {
        y[0] - proj * nadir[0],
        y[1] - proj * nadir[1],
        y[2] - proj * nadir[2]
    };

    double new_y_norm = l2_norm(new_y, 3);

    // Handle degeneracy (y nearly parallel to nadir)
    if (new_y_norm < EPS) {
        // Pick arbitrary orthogonal vector
        if (fabs(nadir[0]) < 0.9) {
            new_y[0] = 0.0;
            new_y[1] = -nadir[2];
            new_y[2] = nadir[1];
        } else {
            new_y[0] = -nadir[1];
            new_y[1] = nadir[0];
            new_y[2] = 0.0;
        }
    }
    normalize3(new_y);

    // Compute x = nadir × new_y
    double x[3];
    cross(nadir, new_y, x);
    normalize3(x);

    // Recompute new_y = x × nadir (ensures orthogonality)
    cross(x, nadir, new_y);

    // Build rotation matrix (row-major)
    double new_R[9];

    new_R[0] = x[0];     new_R[1] = new_y[0];  new_R[2] = nadir[0];
    new_R[3] = x[1];     new_R[4] = new_y[1];  new_R[5] = nadir[1];
    new_R[6] = x[2];     new_R[7] = new_y[2];  new_R[8] = nadir[2];

    // Convert to quaternion
    rotm_to_quat(new_R, goal_q);
}
