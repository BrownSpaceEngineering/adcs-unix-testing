#ifndef TORQUE2MOMENTS
#define TORQUE2MOMENTS
// Minimum-norm dipole m with m x B = (component of torques perpendicular to B).
// Writes zeros if |B| ~ 0.
void torque_2_moments(const float* B, const float* torques, float* moments);
#endif
