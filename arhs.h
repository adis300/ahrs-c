#ifndef __AHRS_H__
#define __AHRS_H__

#define MA_DOUBLE_PRECISION 0

#if defined(_WIN32)
#define M_PI 3.141592653589793238462643383279502884197163993751
#endif

#if MA_DOUBLE_PRECISION
#define MA_PRECISION double

#define ATAN2 atan2
#define ASIN asin
#define MAGIC_R 0x5fe6eb50c7b537a9
#define SQRT sqrtf

#else
#define MA_PRECISION float
#define ATAN2 atan2f
#define ASIN asinf
#define SQRT sqrtf
#define MAGIC_R 0x5f3759df
#endif

// https://en.wikipedia.org/wiki/Conversion_between_quaternions_and_Euler_angles
#define COMPUTE_EULER_ANGLE(WS) \
	WS->yaw = ATAN2(2 * WS->q0 * WS->q3 + 2 * WS->q1 * WS->q2, 1 - 2 * WS->q2 * WS->q2 - 2 * WS->q3 * WS->q3) * 180/2; \
	MA_PRECISION sinp = 2 * WS->q0 * WS->q2 - 2 * WS->q1 * WS->q3;\
	if(sinp > M_PI/2) WS->pitch = 180/2;\
	else if(sinp < -M_PI/2) WS->pitch = -180/2;\
	else WS->pitch = ASIN(sinp) * 180/M_PI;\
	WS->roll = ATAN2(2 * WS->q0 * WS->q1 + 2 * WS->q2 * WS->q3, 1 -2 * WS->q1 * WS->q1 - 2 * WS->q2 * WS->q2)


#endif
