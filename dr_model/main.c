/* main.c
 *
 * Test the DR model
 *
 * This files runs the DR model
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include <time.h>
#include <memory.h>
#include <string.h>
#include "dr_model.h"

#ifdef LINUX
#define _isnan isnan
#endif

int main(int argc, char* argv[])
{
	double eta = 0.05;
	double phi = 0.8;
	int L = 10;
	double Pm = 50;
	int quiet = argc>1 && strcmp(argv[1],"-quiet");
	DRMODEL *drm = drm_create(L,eta,phi);
	int t, T = quiet?1000000:100;
	time_t now = time(NULL), start=now;

	if (drm==NULL)
	{
		printf("%s\n",drm_error());
		return 1;
	}

	if (!quiet)
	{
		printf("time, eta, load, price, forecast\n");
		printf("%4d, %+.2f, %.2f, %.2f, %.2f\n", 0, eta, phi,Pm,Pm);
	}
	for (t=1; t<T; t++)
	{
		double Q, P, Pf;

		/* make a copy of the DR model */
		DRMODEL *old = drm_copy(drm);
		if (old==NULL)
			break;

		/* give a random control signal */
		eta = (rand()/(double)RAND_MAX)*2-1;

		/* calculate new Q */
		Q = drm_update(drm,0,eta,phi);
		if (_isnan(Q))
			break;

		/* calculate price needed to reach new Q from old Q */
		P = drm_price(old,Q,phi,-0.1) + Pm;
		Pf = drm_forecast(drm,0,-0.1,0,eta,phi)+Pm;
		if (!quiet)
			printf("%4d, %+.2f, %.2f, %.2f, %.2f\n", t, eta, Q, P,Pf);
		else if (now!=time(NULL))
		{
			now = time(NULL);
			printf("%d done\n", t);
		}

		/* destroy the old copy */
		free(old);
	}
	if (t<T)
		printf("%4d, %.2f, %s\n", t, eta, drm_error());
	else if (quiet)
	{
		printf("%d timesteps tested ok in %d seconds\n", t, time(NULL)-start);
		printf("estimated performance is %.0f,000 timesteps/second\n", ((double)t/(time(NULL)-start)+500)/1000);
	}
	return 0;
}
